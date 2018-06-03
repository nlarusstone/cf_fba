import pandas as pd
import numpy as np
from cobra.flux_analysis import sample
import multiprocessing
import cobra
import Bio.PDB.Polypeptide as pdb
import Bio.SeqUtils as su
import argparse
import utils
import flux_sample as fs

parser = argparse.ArgumentParser(description='Sample fluxes from different CF models')
parser.add_argument('-d', '--dataset', type=str, default='nls')
parser.add_argument('-f', '--froot', type=str, default='hand')
parser.add_argument('-r', '--rxn', type=int, default=5, help='Amount in uL for a reaction')
parser.add_argument('-a', '--addl', type=float, default=1, help='Additional amount of reactants added')
parser.add_argument('-b', '--batch_size', type=int, default=50, help='Total amount of liquid in uL in a batch')
parser.add_argument('-t', '--txtl', dest='txtl', help='Correlation loss', default=False, action='store_true')
parser.add_argument('--no-txtl', dest='txtl', help='Correlation loss', action='store_false')

if __name__ == '__main__':
    args = parser.parse_args()
    rxn_amt = args.rxn
    addl_amt = args.addl
    batch_size = args.batch_size
    n_batches = batch_size / rxn_amt
    model_f = '{0}_ecoli_cf_base{1}.sbml'.format(args.dataset, '_txtl' if args.txtl else '')

    print 'Read in data'
    df = fs.get_exp_data(args.froot)
    #df.drop(columns=['Area_1', 'Area_2', 'Conc_1', 'Conc_2'], inplace=True)

    cfps_conc = pd.read_csv('../data/{0}_concs'.format(args.dataset), index_col='compound')
    if not args.froot == 'karim':
        cfps_conc.drop('final_conc', inplace=True, axis=1)
        nrg_mix = pd.read_csv('../data/energy_mix.csv', index_col='compound')
    else:
        nrg_mix = None
    print cfps_conc

    print 'Time to generate model specific conditions'
    model = cobra.io.read_sbml_model('../models/{0}'.format(model_f))
    print 'Model {0} read in'.format(model_f)
    res = []
    for idx, row in df.iterrows():
        print idx
        cfps_conc_tmp = cfps_conc.copy()
        cfps_conc_tmp = fs.update_vals(cfps_conc_tmp, row, n_batches, nrg_mix=nrg_mix, replace=args.froot == 'karim')
        tot_rxn_amt = rxn_amt + addl_amt
        if args.froot == 'karim':
            exp_concs = cfps_conc_tmp['final_conc']
        else:
            exp_concs = fs.conc_dilution(cfps_conc_tmp['start_conc'], (rxn_amt / 5.0) * cfps_conc_tmp['amt'], tot_rxn_amt)
        model_i = fs.change_conc(model, exp_concs)
        #print model_i.reactions.EX_mg2_c.bounds
        #print model_i.reactions.EX_pi_c.bounds
        sol = model_i.optimize()
        res.append(sol.objective_value)
    res_df = pd.DataFrame(res)
    res_df['Actual'] = df['OUT']
    res_df.to_csv('../results/{0}_{1}_unoptimized_results{2}'.format(args.froot, args.dataset, '_txtl' if args.txtl else ''))
