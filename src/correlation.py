import argparse
import cobra
import pandas as pd
import cf_io

parser = argparse.ArgumentParser(description='Sample fluxes from different CF models')
parser.add_argument('-d', '--dataset', type=str, default='nls')
parser.add_argument('-f', '--froot', type=str, default='manual')
parser.add_argument('-m', '--model', type=str, default='iJO1366')
parser.add_argument('-t', '--txtl', dest='txtl', help='Toggle to add txtl reactions', default=False, action='store_true')
parser.add_argument('--no-txtl', dest='txtl', help='Toggle to not add txtl reactions', action='store_false')

if __name__ == '__main__':
    args = parser.parse_args()
    model_f = '../bio_models/{0}/{1}/{2}_cf{3}.sbml'.format(args.dataset, args.froot, args.model, '_txtl' if args.txtl else '')

    df = cf_io.get_exp_data(args.froot)
    res = []
    for idx, row in df.iterrows():
        model_f_l = model_f.rsplit('.', 1)
        model_i = cobra.io.read_sbml_model(model_f_l[0] + '_' + str(idx) + '.' + model_f_l[1])
        #print model_i.reactions.EX_mg2_c.bounds
        #print model_i.reactions.EX_pi_c.bounds
        sol = model_i.optimize()
        res.append(sol.objective_value)
    res_df = pd.DataFrame(res)
    res_df['Actual'] = df['OUT']
    res_df.to_csv('../results/{0}_{1}_{2}_unoptimized_results{3}'.format(args.froot, args.dataset, args.model, '_txtl' if args.txtl else ''))
