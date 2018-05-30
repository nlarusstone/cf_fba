import cobra
import cobra.test
import numpy as np
from functools import partial

def convert_cmpts(model, metab, from_cmpt='c', to_cmpt='e'):
    if type(metab) is cobra.Metabolite:
        metab = metab.id.rsplit('_', 1)[0]
    metab_name = metab + '_{0}'.format(to_cmpt)
    try:
        metab_cnvrt = model.metabolites.get_by_id(metab_name)
    except KeyError:
        metab_cnvrt = model.metabolites.get_by_id(metab + '_{0}'.format(from_cmpt)).copy()
        metab_cnvrt.id, metab_cnvrt.compartment = metab_name, to_cmpt
    return metab_cnvrt

def gen_metab_dict(model, metab_nms, metab_vals=None, cnvt=True, from_cmpt='c', to_cmpt='e'):
    metabs = []
    for metab_nm in metab_nms:
        if cnvt:
            metab = convert_cmpts(model, metab_nm, from_cmpt, to_cmpt)
        else:
            metab = metab_nm
        metabs.append(metab)
    if not metab_vals:
        metab_vals = [0] * len(metabs)
    metab_conc = dict(zip(metabs, metab_vals))
    return metab_conc

def add_exchange(model, metab_dict, additive=False):
    rxn_dict = {}
    for med_met, conc in metab_dict.items():
        if conc < 0:
            ty = 'exchange'
            ty_name = 'EX'
	    lb, ub = conc, 0
        else:
            ty = 'demand'
            ty_name = 'DM'
	    lb, ub = conc, 1000
        rxn_nm = '{0}_'.format(ty_name) + med_met.id
        if rxn_nm in model.reactions:
            rxn = model.reactions.get_by_id(rxn_nm)
            if additive:
                rxn.lower_bound,rxn.upper_bound = lb + rxn.lower_bound, ub + rxn.upper_bound
            #    conc += rxn.upper_bound
            #    model.reactions.remove(rxn)
            #    model.add_boundary(metabolite=med_met, type='exchange', lb=conc)
            else:
                rxn.lower_bound, rxn.upper_bound = lb, 0
        else:
            model.add_boundary(metabolite=med_met, type='medium-diff', reaction_id=rxn_nm, lb=lb, ub=0)
        rxn_dict[model.reactions.get_by_id(rxn_nm)] = conc
        #print model.reactions.get_by_id(rxn_nm).upper_bound
    return rxn_dict

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

def build_medium(mod):
    model = mod.copy()
    # in mM
    aa_conc = -1 * conc_to_flux(2)
    trna_conc = -1 * conc_to_flux(0.17)
    metab_conc = map(lambda x: conc_to_flux(x) * -1, [1.2, 0.85, 0.85, 0.85, 0.33, 0.27, 1.50, 1.00, 130, 10, 12, 0.33])
    # TODO: map from conc to exchange fluxes
    
    aa3 = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
    aas = map(lambda x: x.lower(), aa3)
    aa_concs = {}
    for aa in aas:
        l_aa, d_aa = None, None
        try:
            l_aa = convert_cmpts(model, model.metabolites.get_by_id(aa + '__L_c'), from_cmpt='c', to_cmpt='e')
            d_aa = convert_cmpts(model, model.metabolites.get_by_id(aa + '__D_c'), from_cmpt='c', to_cmpt='e')
        except KeyError:
            if aa == 'gly':
                l_aa = convert_cmpts(model, model.metabolites.get_by_id(aa + '_c'), from_cmpt='c', to_cmpt='e')
        if d_aa:
            aa_concs[l_aa] = aa_conc / 2.0
            aa_concs[d_aa] = aa_conc / 2.0
        else:
            aa_concs[l_aa] = aa_conc
            
    trna_concs = {convert_cmpts(model, trna_metab, from_cmpt='c', to_cmpt='e') : trna_conc for trna_metab in model.metabolites.query('trna')}
    
    metab_nms = ['atp', 'gtp', 'utp', 'ctp', 'nad', 'coa', 'spmd', 'ptrc', 'k', 'nh4', 'mg2', 'pep']
    metab_concs = gen_metab_dict(model, metab_nms, metab_conc)
        
    # Not included: RNAP & Folinic Acid & cell extract
    base_medium = {}
    base_medium.update(add_exchange(model, metab_concs))
    base_medium.update(add_exchange(model, aa_concs))
    base_medium.update(add_exchange(model, trna_concs))
    
    #model.medium = base_medium
    return model

def add_addl_reactants(model, df):
    #mod = model.copy()
    addl_reagent_nms = ['mg2', 'nh4', 'k', 'glc__D', 'pi', 'nad', 'atp', 'coa']
    objs = []
    for i, row in df.iterrows():
	with model as mod:
	    metab_dict = gen_metab_dict(mod, addl_reagent_nms, map(lambda x: -1 * conc_to_flux(x), row[4:]))
	    rxn = add_exchange(mod, metab_dict, additive=True)
	    #mod.add_reactions(reaction_list=sol[0])
	    obj = mod.slim_optimize()
            if np.isnan(obj):
                obj = 0
	    objs.append(obj)
	    #print 'Obj: {0}'.format(obj.objective_value)
    return objs

def gen_fluxes_addl_reactants(model, df, parallel):
    #mod = model.copy()
    #objs, fluxes = [], []
    #for i, row in df.iterrows():
    #res = parallel(delayed(get_fluxes)(model, row) for i, row in df.iterrows())
    if parallel:
        res = parallel.map(partial(get_fluxes, model), [i[1] for i in df.iterrows()])
    else:
        res = []
        for i, row in df.iterrows():
            res.append(get_fluxes(model, row))
    return [i[0] for i in res], [i[1] for i in res]

def gen_fluxes(model):
    obj = 0
    with model as mod:
        sol = mod.optimize()
        obj = sol.objective_value
        fluxes = sol.fluxes
    return obj, fluxes

def get_fluxes(model, row):
    addl_reagent_nms = ['mg2', 'nh4', 'k', 'glc__D', 'pi', 'nad', 'atp', 'coa']
    with model as mod:
        metab_dict = gen_metab_dict(mod, addl_reagent_nms, map(lambda x: -1 * conc_to_flux(x), row[4:]))
        rxn = add_exchange(mod, metab_dict, additive=True)
        #mod.add_reactions(reaction_list=sol[0])
        sol = mod.optimize()
        if sol.status == 'infeasible' or np.isnan(sol.objective_value):
            obj =append(0)
        else:
            obj = sol.objective_value
        flux = sol.fluxes
    return obj, flux

def change_obj(model, metab):
    if (not 'DM_{0}'.format(metab.id) in model.reactions):
	model.add_boundary(metabolite=metab, type='demand')
    model.objective = model.reactions.get_by_id('DM_{0}'.format(metab.id))

def different_mediums(model1, model2):
    joint_keys = set(model1.medium.keys()).union(set(model2.medium.keys()))
    for key in joint_keys:
        if key in model2.medium and key in model1.medium:
            print 'both: ', key, model1.medium[key], model2.medium[key]
        elif key in model1.medium:
            print 'Model 1: ', key, model1.medium[key]
        elif key in model2.medium:
            print 'Model 2: ', key, model2.medium[key]

def add_reagents_to_model(model, row):
    mod = model.copy()
    addl_reagent_nms = ['mg2', 'nh4', 'k', 'glc__D', 'pi', 'nad', 'atp', 'coa']
    metab_dict = gen_metab_dict(mod, addl_reagent_nms, map(lambda x: -1 * conc_to_flux(x), row[4:]))
    rxn = add_exchange(mod, metab_dict, additive=True)
    return mod

def get_aa_metab(model, aa, cmpt='c'):
    return model.metabolites.query('{0}__._{1}'.format(aa, cmpt))
