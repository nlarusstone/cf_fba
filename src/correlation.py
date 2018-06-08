import argparse
import cobra
import numpy as np
import pandas as pd
import scipy.stats
import cf_io

parser = argparse.ArgumentParser(description='Sample fluxes from different CF models')
parser.add_argument('-d', '--dataset', type=str, default='nls')
parser.add_argument('-f', '--froot', type=str, default='manual')
parser.add_argument('-m', '--model', type=str, default='iJO1366')
parser.add_argument('-t', '--txtl', dest='txtl', help='Toggle to add txtl reactions', default=False, action='store_true')
parser.add_argument('--no-txtl', dest='txtl', help='Toggle to not add txtl reactions', action='store_false')
#parser.add_argument('-o', '--opt', type=str)

# Performs simple thresholding using the output of a VAE to reduce GEMs to cell-free models.
def optimized_corr(model_f, df, test_dec, cols, t_low=-2, t_high=-1, n_steps=20):
    thresh = np.logspace(start=t_low, stop=t_high, num=n_steps)
    bad_rxns_t = []
    res = []
    best = (0, None)
    model_f_base = model_f.rsplit('.', 1)
    for t in thresh:
        print t
        objs = []
        bad_rxns = []
        for idx, row in df.iterrows():
            model_f_i = model_f_base[0] + '_{0}.'.format(idx) + model_f_base[1]
            model_i = cobra.io.read_sbml_model(model_f_i)
            obj_name = str(model_i.objective.expression.args[0]).split('*')[1]
            dec_df = pd.DataFrame(data=test_dec[:, idx, :], columns=cols)
            bad_cols = cols[dec_df.mean() < t]
            bad_cols = bad_cols[bad_cols != obj_name]
            bad_rxns.append(bad_cols)
            model_i.remove_reactions(bad_cols)
            sol = model_i.optimize()
            objs.append(sol.objective_value)
        bad_rxns_t.append(bad_rxns)
        corr = scipy.stats.pearsonr(objs, df['OUT'])
        print corr
        if abs(corr[0]) > best[0]:
            best = (abs(corr[0]), bad_cols)
        res.append(objs)
    print 'Best: {0}'.format(best[0])
    return pd.DataFrame(res), pd.DataFrame(bad_rxns_t), pd.DataFrame(best[1])

# For base CFPS-like GEMs, calculates their correlation with experimental data
def unoptimized_corr(model_f, df):
    res = []
    for idx, row in df.iterrows():
        model_f_l = model_f.rsplit('.', 1)
        model_i = cobra.io.read_sbml_model(model_f_l[0] + '_' + str(idx) + '.' + model_f_l[1])
        sol = model_i.optimize()
        res.append(sol.objective_value)
    res_df = pd.DataFrame(res)
    res_df['Actual'] = df['OUT']
    return res_df

if __name__ == '__main__':
    args = parser.parse_args()
    model_f = '../bio_models/{0}/{1}/{2}_cf{3}.sbml'.format(args.dataset, args.froot, args.model, '_txtl' if args.txtl else '')

    df = cf_io.get_exp_data(args.froot)
    res_df = unoptimized_corr(model_f, df)
    res_df.to_csv('../results/{0}_{1}_{2}_unoptimized_results{3}'.format(args.froot, args.dataset, args.model, '_txtl' if args.txtl else ''))
