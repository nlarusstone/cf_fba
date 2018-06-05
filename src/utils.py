from functools import partial
import Bio.SeqUtils as su
import cobra
import cobra.test
import numpy as np
from sklearn.preprocessing import maxabs_scale, minmax_scale, normalize, scale, robust_scale
from sklearn.model_selection import train_test_split

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

def add_but(model):
    alc_dehydr = cobra.Reaction(id='ALCDBUT', name='Alcohol dehydrogenase (butanal)', subsystem='c')
    model.add_reaction(alc_dehydr)
    alc_dehydr.add_metabolites({'btcoa_c': -1, 'h_c': -1, 'nadh_c': -1, 'nad_c': 1, 'btal_c': 1})

    butanol = cobra.Metabolite(id='btol_c', name='1-butanol', compartment='c', charge=0, formula='C4H9OH')
    but_synth = cobra.Reaction(id='BUTSYN', name='Butanol synthesis', subsystem='c')
    model.add_reaction(but_synth)
    but_synth.add_metabolites({'btal_c': -1, 'h_c': -1, 'nadh_c': -1, 'nad_c': 1, butanol: 1})
    return model

# amt in g, vol in mL, mw in g/mol
def calc_conc(amt, vol, mw=None, seq=None, seq_type=None):
    # seq can be DNA or protein or an amino acid
    if seq:
        mw = su.molecular_weight(seq, seq_type)
    elif not mw:
        raise Exception('Need a molecular weight for non-DNA')
    conc = (amt * 1000) / (vol * mw)
    # returns Molar concentrations
    return conc

def conc_dilution(start_conc, vol_add, tot_vol):
    return start_conc * (vol_add / tot_vol)

def scale_data(data, scale_type, in_place=True):
    if not in_place:
        scaled_data = data.copy()
    else:
        scaled_data = data
    if 'norm' in scale_type:
        scale_func = scale
    elif 'robust' in scale_type:
        scale_func = robust_scale
    elif 'maxabs' in scale_type:
        scale_func = maxabs_scale
    elif 'negone' in scale_type:
        scale_func = partial(minmax_scale, feature_range=(-1, 1))
    else:
        scale_func = minmax_scale
    if len(data.shape) < 3:
        print len(data.shape)
        scaled_data = scale_func(scaled_data, copy=(not in_place))
    else:
        if 'flux' in scale_type:
            for i in range(data.shape[2]):
                    scaled_data[:, :, i] = scale_func(scaled_data[:, :, i])
        elif 'exp' in scale_type:
            for i in range(data.shape[1]):
                    scaled_data[:, i, :] = scale_func(scaled_data[:, i, :])
        else:
            raise('No scale direction')
    
    if not in_place:
        return scaled_data

def resample(data, n_rows, n_experiments, n_rxns):
    resamp_data = np.empty((n_rows, n_experiments, n_rxns))
    for j in range(n_rows):
        inds = np.random.choice(range(data.shape[0]), size=n_experiments, replace=True)
        for i in range(n_experiments):
            ind = inds[i]
            flxs = data[ind, i, :]
            resamp_data[j][i][:] = flxs
        if j % (n_rows / 10) == 0:
            print j
    print 'done'
    return resamp_data

def biased_resample(sorted_samp_data, n_rows=2000):
    #sorted_samp_data = samp_data_scaled.copy()
    #for i in range(n_experiments):
    #    exp = samp_data_scaled[:, i, :]
    #    sorted_samp_data[:, i, :] = samp_data_scaled[exp[:,btol_col].argsort(), i, :]
    np.random.seed(42)
    resamp_data = np.empty((n_rows, n_experiments, max_sz[0]))
    for j in range(n_rows):
        #exps = []
        btol_val = 0
        for i in range(n_experiments):
            inds = (sorted_samp_data[:, i, btol_col] >= btol_val).nonzero()[0]
            if inds.any():
                ind = np.random.choice(inds[:5], size=1, replace=True)[0]
            else:
                ind = np.argmax(sorted_samp_data[:, i, btol_col])
            btol_val = sorted_samp_data[ind, i, btol_col]
            #exps.append(btol_val)
            flxs = sorted_samp_data[ind, i, :]
            resamp_data[j][i][:] = flxs
        #print scipy.stats.pearsonr(exps, df['AVG.1'])
        if j % (n_rows / 10) == 0:
            print j
    print 'done'
    return resamp_data

def gen_train_test(data, y=None):
    #train_ind = np.random.choice(data.shape[0], size=int(0.9 * data.shape[0]), replace=False)
    #test_ind = list(set(range(data.shape[0])) - set(train_ind))
    #y = np.array(range(41) * data.shape[0])
    if not y is None:
        return train_test_split(data, y, random_state=42)
    else:
        return train_test_split(data, random_state=42)

