import Bio.PDB.Polypeptide as pdb
import cobra
import re
import utils
import cf_io

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
        flux = utils.conc_to_flux(vals['final_conc']) * 100

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

if __name__ == '__main__':
    model = cobra.io.read_sbml_model(filename='../models/iJO1366.xml')
    model_cyt = coalesce_cmpts(model)
    cfps_conc = cf_io.get_conc()
    model_bare = strip_exchanges(model_cyt, cfps_conc.index[:-1])
    model_cf = build_medium(model_bare, cfps_conc)
    cobra.io.write_sbml_model(filename='../models/ecoli_simp_cf_base.sbml', cobra_model=model_cf)
