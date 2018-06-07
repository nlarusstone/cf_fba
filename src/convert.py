import argparse
import difflib
import os
import re
import Bio.PDB.Polypeptide as pdb
import cobra
import pandas as pd
import cf_io
from constants import varner_to_ijo
import fba_utils as futils
import utils

def coalesce_cmpts(model):
    mod = model.copy()
    for rxn in mod.reactions:
        if rxn.compartments != set('c'):
        #if 'p' in rxn.compartments or 'e' in rxn.compartments:
            #mod.remove_reactions(reactions=[rxn])
            for metab, amt in rxn.metabolites.items():
                cyt = futils.replace_metab(mod, metab)
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
    model = mod.copy()
    exs = set()
    for metab in reactants:
        if metab == 'trna':
            for trna in model.metabolites.query('trna'):
                exs = exs.union(trna.reactions.intersection(model.exchanges))
        elif metab.upper() in pdb.aa3:
            aas = futils.get_aa_metab(model, metab.lower(), cmpt='c')
            for aa in aas:
                exs = exs.union(aa.reactions.intersection(model.exchanges))
        else:
            try:
                m = model.metabolites.get_by_id('{0}_c'.format(metab))
                exs = exs.union(m.reactions.intersection(model.exchanges))
            except:
                continue
    model.remove_reactions(exs)
    return model

def build_medium(model, cfps_conc):
    mod = model.copy()
    for metab, vals in cfps_conc.iterrows():
        flux = utils.conc_to_flux(vals['final_conc'])

        if metab == 'trna':
            ms = model.metabolites.query('trna')
        elif metab.upper() in pdb.aa3:
            ms = futils.get_aa_metab(model, metab.lower(), cmpt='c')
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

def add_txtl_rxns(model, txtl_rxns):
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

def convert_to_cf_model(model_f, add_txtl, obj='BIOMASS_Ec_iJO1366_core_53p95M', cfps_sys='nls', conc_file=None):
    final_concs_fname = '../bio_models/{0}/final_concs.csv'.format(cfps_sys)
    if os.path.exists(final_concs_fname):
        print 'Read in final concentrations'
        cfps_conc = pd.read_csv(final_concs_fname, index=compound)
    else:
        print 'Generate final concentrations'
        cfps_conc = cf_io.get_conc(cfps_final=conc_file)
        print 'Writing out CFPS system "{0}" concentrations to {1}'.format(cfps_sys, final_concs_fname)
        cfps_conc.to_csv(path_or_buf=final_concs_fname)

    model_base = model_f.rsplit('/', 1)[1].rsplit('.', 1)[0]
    print 'Read in GEM: {0}'.format(model_base)
    if model_f.endswith('json'):
        model = cobra.io.load_json_model(model_f)
    else:
        model = cobra.io.read_sbml_model(filename=model_f)

    if add_txtl:
        print 'Adding TXTL reactions'
        varner = cobra.io.load_json_model('../bio_models/varner.json')
        txtl_rxns = extract_txtl_rxns(varner)
        model = add_txtl_rxns(model, txtl_rxns)

    if cfps_sys == 'karim':
        model = utils.add_but(model)

    print 'Moving all reactions to same compartment'
    model_cyt = coalesce_cmpts(model)
    print 'Rebuilding medium'
    model_bare = strip_exchanges(model_cyt, cfps_conc.index[:-1])
    model_cf = build_medium(model_bare, cfps_conc)

    print 'Updating objective to {0}'.format(obj)
    try:
        futils.change_obj(model, metab=obj)
    except:
        futils.change_obj(model, rxn=obj)

    cf_model_fname = '../bio_models/{0}/{1}_cf{2}.sbml'.format(cfps_sys, model_base, '_txtl' if add_txtl else '')
    print 'Writing out cf_model to {0}'.format(cf_model_fname)
    cobra.io.write_sbml_model(filename=cf_model_fname, cobra_model=model_cf)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create cell free model')
    parser.add_argument('-c', '--conc-file', metavar='c', type=str, help='Path for final concentration file', default=None)
    parser.add_argument('-t', '--txtl', dest='txtl', help='Toggle to add txtl reactions', default=False, action='store_true')
    parser.add_argument('--no-txtl', dest='txtl', help='Toggle to not add txtl reactions', action='store_false')
    parser.add_argument('-o', '--obj', type=str, default='BIOMASS_Ec_iJO1366_core_53p95M')
    parser.add_argument('-s', '--sys', type=str, default='nls')
    parser.add_argument('-m', '--model', type=str, help='Path to model file', default='../bio_models/iJO1366.xml')

    args = parser.parse_args()

    if not os.path.exists('../bio_models/{0}'.format(args.sys)):
        os.mkdir('../bio_models/{0}'.format(args.sys))
    convert_to_cf_model(args.model, args.txtl, args.obj, args.sys, args.conc_file)
