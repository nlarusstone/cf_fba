import argparse
import multiprocessing
import os
import Bio.PDB.Polypeptide as pdb
import cobra
from cobra.flux_analysis import sample
import numpy as np
import pandas as pd
import cf_io
import fba_utils as futils
import utils

parser = argparse.ArgumentParser(description='Sample fluxes from different CF models')
parser.add_argument('-s', '--samps', metavar='s', type=int, help='Number of samples, set to 0 if you want to generate models without sampling', default=0)
parser.add_argument('-d', '--dataset', type=str, default='nls')
parser.add_argument('-f', '--froot', type=str, default='manual')
parser.add_argument('-m', '--model', type=str, default='iJO1366')
parser.add_argument('-r', '--rxn', type=int, default=5, help='Amount in uL for a reaction')
parser.add_argument('-a', '--addl', type=float, default=1, help='Additional amount of reactants added')
parser.add_argument('-b', '--batch_size', type=int, default=50, help='Total amount of liquid in uL in a batch')
parser.add_argument('-t', '--txtl', dest='txtl', help='Toggle to add txtl reactions', default=False, action='store_true')
parser.add_argument('--no-txtl', dest='txtl', help='Toggle to not add txtl reactions', action='store_false')

# Helper function to change the exchange reactions of a model
def update_vals(cfps_conc_tmp, row, n_batches, nrg_mix=None, replace=False):
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
            if replace:
                cfps_conc_tmp.loc[col]['final_conc'] = row[col] / 1000.0
            else:
                cfps_conc_tmp.loc[col]['amt'] += row[col]
    return cfps_conc_tmp

if __name__ == '__main__':
    args = parser.parse_args()
    rxn_amt = args.rxn
    addl_amt = args.addl
    batch_size = args.batch_size
    n_batches = batch_size / rxn_amt

    print 'Read in data'
    df = cf_io.get_exp_data(args.froot)
    cfps_conc = pd.read_csv('../bio_models/{0}/final_concs.csv'.format(args.dataset), index_col='compound')
    if not args.froot == 'karim':
        cfps_conc.drop('final_conc', inplace=True, axis=1)
        nrg_mix = pd.read_csv('../data/energy_mix.csv', index_col='compound')
    else:
        nrg_mix = None

    print 'Generate model specific conditions'
    model_base = '{0}_cf{1}.sbml'.format(args.model, '_txtl' if args.txtl else '')
    model_f = '../bio_models/{0}/{1}'.format(args.dataset, model_base)
    model_i_base = '../bio_models/{0}/{1}/{2}'.format(args.dataset, args.froot, model_base)
    model = cobra.io.read_sbml_model(model_f)
    print 'Model {0} read in'.format(model_f)
    if not os.path.exists('../bio_models/{0}/{1}/'.format(args.dataset, args.froot)):
        os.makedirs('../bio_models/{0}/{1}/'.format(args.dataset, args.froot))
    if not os.path.exists('../data/f{0}/{1}/{2}/'.format(args.samps, args.dataset, args.froot)):
        os.makedirs('../data/f{0}/{1}/{2}/'.format(args.samps, args.dataset, args.froot))
    for idx, row in df.iterrows():
        print idx
        cfps_conc_tmp = cfps_conc.copy()
        cfps_conc_tmp = update_vals(cfps_conc_tmp, row, n_batches, nrg_mix=nrg_mix, replace=args.froot == 'karim')
        tot_rxn_amt = rxn_amt + addl_amt
        if args.froot == 'karim':
            exp_concs = cfps_conc_tmp['final_conc']
        else:
            exp_concs = utils.conc_dilution(cfps_conc_tmp['start_conc'], (rxn_amt / 5.0) * cfps_conc_tmp['amt'], tot_rxn_amt)
        model_i = futils.change_conc(model, exp_concs)
        model_i_prefix = model_i_base.rsplit('.', 1)[0]
        cobra.io.write_sbml_model(filename=model_i_prefix + '_' + str(idx) + '.sbml', cobra_model=model_i)
        print model_i.medium
        if args.samps > 0:
            samples_i = sample(model_i, args.samps, processes=multiprocessing.cpu_count() - 1)
            samples_i.to_csv('../data/f{0}/{1}/{2}/{3}_fluxes_{4}'.format(args.samps, args.dataset, args.froot, model_base, idx))
