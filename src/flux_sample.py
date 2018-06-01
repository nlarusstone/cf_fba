import pandas as pd
import numpy as np
from cobra.flux_analysis import sample
import multiprocessing
import cobra
import utils
import Bio.PDB.Polypeptide as pdb
import Bio.SeqUtils as su
import argparse

parser = argparse.ArgumentParser(description='Sample fluxes from different CF models')
parser.add_argument('-s', '--samps', metavar='s', type=int, help='Number of samples', default=200)
parser.add_argument('-d', '--dataset', type=str, default='nls')
parser.add_argument('-f', '--froot', type=str, default='hand')
parser.add_argument('-r', '--rxn', type=int, default=5, help='Amount in uL for a reaction')
parser.add_argument('-a', '--addl', type=float, default=1, help='Additional amount of reactants added')
parser.add_argument('-b', '--batch_size', type=int, default=50, help='Total amount of liquid in uL in a batch')

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

def conc_dilution(start_conc, vol_add, tot_vol):
    return start_conc * (vol_add / tot_vol)

def get_exp_data(froot):
    # '../data/17_5_18_T7_mRFP_NLS.CSV'
    df = pd.read_csv('../data/{0}_data.CSV'.format(froot), skiprows=6)
    print df.shape
    gain_diff = df.shape[0] / 5
    times = df["Unnamed: 1"]
    df.drop('Unnamed: 1', inplace=True, axis=1)
    # Bad data
    if froot == 'hand':
        df.drop('E09', inplace=True, axis=1)
    # Remove negative control
    df.drop(df.columns[-2:], inplace=True, axis=1)
    
    gain2 = gain_diff
    outs = df[gain2:gain2+gain_diff - 1].mean(axis=0)
    # '../data/17_5_18_exp_conditions.csv'
    conds = pd.read_csv('../data/{0}_exp_conditions.csv'.format(froot))
    if froot == 'hand':
        conds.drop(conds.shape[0] - 1, inplace=True)
        conds_full = conds.reindex(np.repeat(conds.index.values, 2)).reset_index()
        conds_full = conds_full.drop(32).reset_index(drop=True)
    else:
        conds_full = conds.reindex(np.repeat(conds.index.values, 2)).reset_index()
    conds_full['OUT'] = outs.reset_index(drop=True)
    conds_avg = conds_full.groupby('index').mean()
    conds_norm = conds_avg / conds_avg.max()
    conds_norm.to_csv('../data/{0}_EXPERIMENT.csv'.format(froot))
    return conds_norm

def update_vals(cfps_conc_tmp, row, n_batches):
    for col in row.index[:-1]:
        if col == 'nts':
            for nt in ['atp', 'gtp', 'ctp', 'utp']:
                cfps_conc_tmp.loc[nt]['amt'] += row['nts'] / 4
        elif col == 'mdx':
            for cmpnd, vals in nrg_mix.iterrows():
                cfps_conc_tmp.loc[cmpnd] = [cfps_conc_tmp.loc[cmpnd]['start_conc'], (5.0 + row['mdx'] / n_batches)]
        elif col == 'aas':
            for aa in pdb.aa3:
                aa = aa.lower()
                cfps_conc_tmp.loc[aa] = [cfps_conc_tmp.loc[aa]['start_conc'], (10.0 + row['aas'] / n_batches)]
        else:
            cfps_conc_tmp.loc[col]['amt'] += row[col]
    return cfps_conc_tmp

if __name__ == '__main__':
    args = parser.parse_args()
    flux_sz = args.samps
    rxn_amt = args.rxn
    addl_amt = args.addl
    batch_size = args.batch_size
    n_batches = batch_size / rxn_amt
    models = ['{0}_ecoli_cf_base{1}.sbml'.format(args.dataset, txtl) for txtl in ['_txtl', '']]

    print 'Read in data'
    df = get_exp_data(args.froot) #pd.read_csv('../data/Karim_MetEng_2018_Figure2_Data.csv')
    #df.drop(columns=['Area_1', 'Area_2', 'Conc_1', 'Conc_2'], inplace=True)

    cfps_conc = pd.read_csv('../data/{0}_concs'.format(args.dataset), index_col='compound')
    #cfps_conc.index = cfps_conc.index.map(lambda x: x.lower())
    cfps_conc.drop('final_conc', inplace=True, axis=1)
    print cfps_conc
    nrg_mix = pd.read_csv('../data/energy_mix.csv', index_col='compound')

    print 'Time to generate model specific conditions'
    for model_f in models:
        model = cobra.io.read_sbml_model('../models/{0}'.format(model_f))
        print 'Model {0} read in'.format(model_f)
        for idx, row in df.iterrows():
            print idx
            cfps_conc_tmp = cfps_conc.copy()
            cfps_conc_tmp = update_vals(cfps_conc_tmp, row, n_batches)
            rxn_amt += addl_amt
            exp_concs = conc_dilution(cfps_conc_tmp['start_conc'], cfps_conc_tmp['amt'], rxn_amt)
            model_i = change_conc(model, exp_concs)
            print model_i.medium
            samples_i = sample(model_i, flux_sz, processes=multiprocessing.cpu_count() - 1)
            samples_i.to_csv('../data/f{2}/{3}_{0}_fluxes_{1}'.format(model_f, idx, flux_sz, args.froot))
