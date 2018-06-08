import re
import Bio.PDB.Polypeptide as pdb
import cobra

# Gets the metabolite associated with an amino acid
def get_aa_metab(model, aa, cmpt='c'):
    return model.metabolites.query('{0}__._{1}'.format(aa, cmpt))

def replace_metab(model, metab):
    new_id = re.sub(r'_.$', '_c', metab.id)
    try:
        cyt = model.metabolites.get_by_id(new_id)
    except:
        cyt = metab
        cyt.id = new_id
        cyt.compartment = 'c'
    return cyt

# Changes objective function of model
# Works in place
def change_obj(model, metab=None, rxn=None):
    if rxn is not None:
        try:
            model.objective = model.reactions.get_by_id(rxn)
        except:
            raise Exception('Objective reaction does not exist')
    elif metab is not None:
        try:
            met = model.metabolites.get_by_id(metab)
        except:
            raise Exception('Metabolite does not exist')
        if (not 'DM_{0}'.format(met.id) in model.reactions):
            model.add_boundary(metabolite=met, type='demand')
            model.objective = model.reactions.get_by_id('DM_{0}'.format(met.id))
    else:
        print 'Need to specify either a metabolite or reaction as an objective'
        raise Exception

# Change the concentration of exchange reactions to a given cell-free setup
def change_conc(model, cfps_conc):
    mod = model.copy()
    
    for metab, vals in cfps_conc.iteritems():
        flux = conc_to_flux(vals)

        if metab == 'trna':
            ms = model.metabolites.query('trna')
        elif metab.upper() in pdb.aa3:
            ms = get_aa_metab(model, metab.lower(), cmpt='c')
        else:
            ms = mod.metabolites.query(r'^{0}_c'.format(metab))
        for m in ms:
            rxn_nm = 'EX_' + m.id
            rxn = mod.reactions.get_by_id(rxn_nm)
            rxn.lower_bound = -1 * flux
            rxn.upper_bound = flux
    return mod

# Based on MetaboTools
# flux = metConc/(cellConc*cellWeight*t*1000);
# cellConc: 10 mg/ml
# cellWeight: 500 * (10 ** -11)
# t = 24
# t: Experiment duration
def conc_to_flux(metab_conc, t=24):
    # Taken from MetaboTools, units are gDW/cell
    #cell_weight = 500 * (10 ** -11)
    # We have 10 mg/ml
    # need cells / ml
    cell_conc = 10 * (1/ 1000.0) #* (1 / cell_weight)
    flux = metab_conc / (cell_conc * t * 1000.0)
    return flux * 100
