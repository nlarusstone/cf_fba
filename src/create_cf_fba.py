import Bio.PDB.Polypeptide as pdb
import cobra
import re
import utils
import cf_io
import argparse
import difflib
from constants import varner_to_ijo

parser = argparse.ArgumentParser(description='Create cell free model')
parser.add_argument('-c', '--conc-final', metavar='c', type=str, help='Path for final concentration file')
parser.add_argument('-t', '--txtl', dest='txtl', help='Correlation loss', default=False, action='store_true')
parser.add_argument('--no-txtl', dest='txtl', help='Correlation loss', action='store_false')
parser.add_argument('-d', '--dataset', type=str, default='nls')

def get_aa_metab(model, aa, cmpt='c'):
    return model.metabolites.query('{0}__._{1}'.format(aa, cmpt))

def replace_metab(mod, metab):
    new_id = re.sub(r'_.$', '_c', metab.id)
    try:
        cyt = mod.metabolites.get_by_id(new_id)
    except:
        cyt = metab
        cyt.id = new_id
        cyt.compartment = 'c'
    return cyt

def coalesce_cmpts(model):
    mod = model.copy()
    for rxn in mod.reactions:
        if 'p' in rxn.compartments or 'e' in rxn.compartments:
            #mod.remove_reactions(reactions=[rxn])
            for metab, amt in rxn.metabolites.items():
                cyt = replace_metab(mod, metab)
                rxn.add_metabolites({metab: -1 * amt})
                rxn.add_metabolites({cyt: amt})
            rxn.comparments = set('c')
            #mod.add_reaction(reaction=rxn)
    for m in mod.metabolites.query(r'.*_e$'):
        assert(len(m.reactions) == 0)
        m.remove_from_model(destructive=True)
    for m in mod.metabolites.query(r'.*_p$'):
        assert(len(m.reactions) == 0)
        m.remove_from_model(destructive=True)
    return mod

def strip_exchanges(mod, reactants):
    # Delete transmembrane transport reactions
    model = mod.copy()

    exs = set()
    for metab in reactants:
        if metab == 'trna':
            for trna in model.metabolites.query('trna'):
                exs = exs.union(trna.reactions.intersection(model.exchanges))
        elif metab.upper() in pdb.aa3:
            aas = get_aa_metab(model, metab.lower(), cmpt='c')
            for aa in aas:
                exs = exs.union(aa.reactions.intersection(model.exchanges))
        else:
            m = model.metabolites.get_by_id('{0}_c'.format(metab))
            exs = exs.union(m.reactions.intersection(model.exchanges))
    model.remove_reactions(exs)
    #['EX_glc_e', 'EX_pi_e', 'EX_mg2_e', 'EX_k_e', 'EX_nh4_e'])

    # As objective function, we selected the exchange reaction which corresponds to the target metabolite 
    # for which a pathway should be determined.   
    return model

def build_medium(model, cfps_conc):
    mod = model.copy()
    
    for metab, vals in cfps_conc.iterrows():
        flux = utils.conc_to_flux(vals['final_conc'])

        if metab == 'trna':
            ms = model.metabolites.query('trna')
        elif metab.upper() in pdb.aa3:
            ms = get_aa_metab(model, metab.lower(), cmpt='c')
        else:
            ms = mod.metabolites.query(r'^{0}_c'.format(metab))
        for m in ms:
            rxn_nm = 'EX_' + m.id
            mod.add_boundary(metabolite=m, type='exchange', lb=0, ub=flux)
            #mod.add_boundary(metabolite=m, type='cfps-medium', reaction_id=rxn_nm, lb=0, ub=flux) 
    return mod

def extract_txtl_rxns(model):
    aa_metabs = []
    for aa in pdb.aa3:
        aa_metabs += model.metabolites.query(aa.lower())
    aa_rxns = []
    for aa_metab in aa_metabs:
        aa_rxns += aa_metab.reactions
    mrna_rxns = model.reactions.query(re.compile('mrna', re.IGNORECASE))
    trna_rxns = model.reactions.query('tRNA_c')
    tx_rxns = model.reactions.query('transcription')
    tl_rxns = model.reactions.query('translation')
    prot_rxns = model.reactions.query('PROTEIN')
    #txtl_rxns = list(set(aa_rxns).union(tx_rxns).union(tl_rxns).union(prot_rxns).union(mrna_rxns))
    txtl_rxns = list(set(tx_rxns).union(tl_rxns).union(prot_rxns).union(mrna_rxns).union(trna_rxns))
    return txtl_rxns

def varner_to_cobra(model, metab, metab_ids, varner_to_ijo):
    if metab.id.startswith('M_'):
        metab_stem = metab.id.split('M_')[1].rsplit('_c', 1)[0]
        #print metab_stem
        if 'tRNA' in metab_stem:
            aa = metab_stem.split('_', 1)[0]
            metab_name = aa + 'trna'
        elif not metab_stem in metab_ids:
            #query = varner_to_ijo[metab_stem]
            #print metab_stem
            if metab_stem in varner_to_ijo:
                #print 'matched'
                metab_name = varner_to_ijo[metab_stem]
            elif '_L' in metab_stem or '_D' in metab_stem:
                #print difflib.get_close_matches(metab_stem, metab_ids, 1, 0.7)
                metab_name = difflib.get_close_matches(metab_stem, metab_ids, 1, 0.7)[0]
            else:
                print 'TODO: ', metab_stem
                raise Exception
        else:
            metab_name = metab_stem
    else:
        try:
            model.metabolites.get_by_id(metab_name)
        except:
            model.metabolites.add(metab)
    return model.metabolites.get_by_id(metab_name + '_c')

def add_txtl(model, txtl_rxns):
    mod = model.copy()
    metab_ids = [m.id.rsplit('_c', 1)[0] for m in mod.metabolites if m.compartment == 'c']
    for rxn in txtl_rxns:
        #print rxn
        for metab, amt in rxn.metabolites.items():
            if not metab.id.startswith('M_'):
                #print 'EXCEPT:', metab
                continue
            new_metab = varner_to_cobra(mod, metab, metab_ids, varner_to_ijo)
            rxn.add_metabolites({metab: -1 * amt})
            rxn.add_metabolites({new_metab: amt})
        mod.add_reaction(rxn)
    return mod

def add_but(model):
    alc_dehydr = cobra.Reaction(id='ALCDBUT', name='Alcohol dehydrogenase (butanal)', subsystem='c')
    model.add_reaction(alc_dehydr)
    alc_dehydr.add_metabolites({'btcoa_c': -1, 'h_c': -1, 'nadh_c': -1, 'nad_c': 1, 'btal_c': 1})

    butanol = cobra.Metabolite(id='btol_c', name='1-butanol', compartment='c', charge=0, formula='C4H9OH')
    but_synth = cobra.Reaction(id='BUTSYN', name='Butanol synthesis', subsystem='c')
    model.add_reaction(but_synth)
    but_synth.add_metabolites({'btal_c': -1, 'h_c': -1, 'nadh_c': -1, 'nad_c': 1, butanol: 1})
    return model

if __name__ == '__main__':
    args = parser.parse_args()
    print args
    print 'Generate final concentrations'
    cfps_conc = cf_io.get_conc(cfps_final=args.conc_final)
    print 'Read in GEM'
    model = cobra.io.read_sbml_model(filename='../models/iJO1366.xml')
    if args.txtl:
        print 'Adding TXTL reactions'
        varner = cobra.io.load_json_model('../models/varner.json')
        txtl_rxns = extract_txtl_rxns(varner)
        model = add_txtl(model, txtl_rxns)
    if args.dataset == 'karim':
        model = add_but(model)
        utils.change_obj(model, model.metabolites.btol_c)
    print 'Moving all reactions to same compartment'
    model_cyt = coalesce_cmpts(model)
    print 'Rebuilding medium'
    model_bare = strip_exchanges(model_cyt, cfps_conc.index[:-1])
    model_cf = build_medium(model_bare, cfps_conc)
    if args.txtl:
        utils.change_obj(model_cf, model_cf.metabolites.PROTEIN_RFP)
    print 'Writing out'
    cobra.io.write_sbml_model(filename='../models/{0}_ecoli_cf_base{1}.sbml'.format(args.dataset, '_txtl' if args.txtl else ''), cobra_model=model_cf)
    cfps_conc.to_csv(path_or_buf='../data/{0}_concs'.format(args.dataset))
