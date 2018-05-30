import pandas as pd
import numpy as np
from cobra.flux_analysis import sample
import multiprocessing
import cobra
import utils
import Bio.PDB.Polypeptide as pdb
import Bio.SeqUtils as su

def get_aa_metab(model, aa, cmpt='c'):
    return model.metabolites.query('{0}__._{1}'.format(aa, cmpt))

def change_conc(model, cfps_conc):
    mod = model.copy()
    
    for metab, vals in cfps_conc.iteritems():
        flux = utils.conc_to_flux(vals) * 100

        if metab == 'trna':
            ms = model.metabolites.query('trna')
        elif metab.upper() in pdb.aa3:
            ms = get_aa_metab(model, metab.lower(), cmpt='c')
        else:
            ms = mod.metabolites.query(r'^{0}_c'.format(metab))
        for m in ms:
            rxn_nm = 'EX_' + m.id
            rxn = mod.reactions.get_by_id(rxn_nm)
            rxn.upper_bound = flux
            #mod.add_boundary(metabolite=m, type='exchange', lb=0, ub=flux)
            #mod.add_boundary(metabolite=m, type='cfps-medium', reaction_id=rxn_nm, lb=0, ub=flux) 
    return mod

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

def read_start_conds():
    aa_mix = pd.read_csv('../data/aa_mix.csv', index_col='AA')
    nrg_mix = pd.read_csv('../data/energy_mix.csv', index_col='compound')
    with open('../genes/rfp.txt', 'r') as f:
        seq = f.read()
    cfps_conc = pd.read_csv('../data/cfps_start.csv', index_col='compound')
    dna_conc = calc_conc(0.000750, 0.00496, seq=seq, seq_type='DNA')
    return aa_mix, nrg_mix, cfps_conc, dna_conc

def get_exp_data():
    df = pd.read_csv('../data/17_5_18_T7_mRFP_NLS.CSV', skiprows=6)
    print df.shape
    gain_diff = df.shape[0] / 5
    times = df["Unnamed: 1"]
    df.drop('Unnamed: 1', inplace=True, axis=1)
    # Bad data
    df.drop(['E09', 'F04', 'F05'], inplace=True, axis=1)
    df.head()
    gain2 = gain_diff
    outs = df[gain2:gain2+gain_diff - 1].mean(axis=0)
    conds = pd.read_csv('../data/17_5_18_exp_conditions.csv')
    conds.drop(conds.shape[0] - 1, inplace=True)
    conds_full = conds.reindex(np.repeat(conds.index.values, 2)).reset_index()
    conds_full = conds_full.drop(32).reset_index(drop=True)
    conds_full['OUT'] = outs.reset_index(drop=True)
    conds_avg = conds_full.groupby('index').mean()
    conds_norm = conds_avg / conds_avg.max()
    conds_norm.to_csv('../data/17_5_18_EXPERIMENT.csv')
    return conds_norm

if __name__ == '__main__':
    flux_sz = 10
    models = ['ecoli_cf_base.sbml', 'ecoli_cf_txtl_rfp_base.sbml', 'ecoli_cf_txtl_comb_base.sbml']

    df = get_exp_data() #pd.read_csv('../data/Karim_MetEng_2018_Figure2_Data.csv')
    print 'Read in data'
    #df.drop(columns=['Area_1', 'Area_2', 'Conc_1', 'Conc_2'], inplace=True)

    aa_mix, nrg_mix, cfps_conc, dna_conc = read_start_conds()
    rxn_amt = 5
    batch_size = 50 / rxn_amt
    cfps_conc['amt'] = cfps_conc['amt'] / batch_size
    aa_mix['start_conc'] = aa_mix.apply(lambda row: calc_conc(row['weight_add'], 1, 
                                                              seq=pdb.three_to_one(row.name.upper()), seq_type='protein'), axis=1)
    aa_mix['conc_add'] = conc_dilution(aa_mix['start_conc'], aa_mix['vol_add'], aa_mix['vol_add'].sum())
    pi_conc = calc_conc(0.15, 5, mw=611.77)
    nrg_mix['start_conc'] = nrg_mix.apply(lambda row: calc_conc(row['amt'], row['fill'], mw=row['mw']), axis=1)
    nrg_mix['conc_add'] = conc_dilution(nrg_mix['start_conc'], nrg_mix['vol_add'], nrg_mix['vol_add'].sum())

    for cmpnd, vals in nrg_mix.iterrows():
        cfps_conc.loc[cmpnd] = [vals['conc_add'], (5.0 / batch_size)]
    for aa, vals in aa_mix.iterrows():
        cfps_conc.loc[aa] = [vals['conc_add'], (10.0 / batch_size)]
    cfps_conc.loc['GENE'] = [dna_conc, (4.96 / batch_size)]

    print 'Time to generate model specific conditions'
    for model_f in models:
        model = cobra.io.read_sbml_model('../models/{0}'.format(model_f))
        print 'Model {0} read in'.format(model_f)
        for idx, row in df.iterrows():
            print idx
            cfps_conc_tmp = cfps_conc.copy()
            cfps_conc_tmp.loc['pi']['amt'] += row['pi']
            cfps_conc_tmp.loc['k']['amt'] += row['k']
            for nt in ['atp', 'gtp', 'ctp', 'utp']:
                cfps_conc_tmp.loc[nt]['amt'] += row['nts'] / 4
            for cmpnd, vals in nrg_mix.iterrows():
                cfps_conc_tmp.loc[cmpnd] = [vals['conc_add'], (5.0 + row['mdx'] / batch_size)]
            rxn_amt += 1
            exp_concs = conc_dilution(cfps_conc_tmp['start_conc'], cfps_conc_tmp['amt'], rxn_amt)
            model_i = change_conc(model, exp_concs)
            print model_i.medium
            samples_i = sample(model_i, flux_sz, processes=multiprocessing.cpu_count() - 1)
            samples_i.to_csv('../data/f{2}/{0}_fluxes_{1}'.format(model_f, idx, flux_sz))
