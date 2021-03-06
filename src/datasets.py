import argparse
import h5py
import numpy as np
import pandas as pd
import cf_io
import utils

parser = argparse.ArgumentParser(description='Sample fluxes from different CF models')
parser.add_argument('-d', '--dataset', type=str, default='nls')
parser.add_argument('-f', '--froot', type=str, default='manual')
parser.add_argument('-m', '--model', type=str, default='iJO1366')
parser.add_argument('-r', '--resamp', dest='resamp', help='Resample', default=False, action='store_true')
parser.add_argument('--no-resamp', dest='resamp', help='Dont resample', action='store_false')
parser.add_argument('-t', '--txtl', dest='txtl', default=False, action='store_true')
parser.add_argument('--no-txtl', dest='txtl', action='store_false')

# Creates a h5py dataset of fluxes
def create_dataset(fname, X_train, X_test, obj_col, cols, y_vals, y_train=None, y_test=None):
    h5f = h5py.File(fname, 'w')
    h5f.create_dataset('X_train', data=X_train)
    h5f.create_dataset('X_test', data=X_test)
    if not y_train is None:
        h5f.create_dataset('y_train', data=y_train)
        h5f.create_dataset('y_test', data=y_test)
    h5f.create_dataset('obj_col', data=obj_col)
    h5f.create_dataset('cols', data=list(cols))
    h5f.create_dataset('y_vals', data=list(y_vals))
    h5f.close()

# Reads in a h5py dataset of fluxes
def get_dataset(fname):
    h5f = h5py.File(fname, 'r')
    X_train, X_test = h5f['X_train'], h5f['X_test']
    if 'y_train' in h5f:
        y_train, y_test = h5f['y_train'], h5f['y_test']
    else:
        y_train, y_test = None, None
    obj_col = h5f['obj_col']
    cols = h5f['cols']
    y_vals = h5f['y_vals']
    return X_train, y_train, X_test, y_test, obj_col, cols, y_vals

if __name__ == '__main__':
    args = parser.parse_args()
    df = cf_io.get_exp_data(args.froot)
    n_experiments = df.shape[0]

    if args.froot == 'karim':
        if args.txtl:
            raise Exception('No TXTL with this dataset')
        obj_name = 'DM_btol_c'
    elif args.txtl:
        obj_name = 'PROTEIN_export_RFP'
    else:
        obj_name = 'BIOMASS_Ec_iJO1366_core_53p95M'

    flux_root = '{0}_cf{1}.sbml_fluxes'.format(args.model, '_txtl' if args.txtl else '')
    X_train, y_train, X_test, y_test, obj_col, cols = cf_io.read_flux_data('../data/f2000/{0}/{1}'.format(args.dataset, args.froot), flux_root,
                                                                    n_experiments, obj_name, resamp=args.resamp, scale='flux_zero')
    y_vals = utils.scale_data(data=df['OUT'].values, scale_type='flux_zero', in_place=False)
    fname = '../data/{0}{1}_{2}_fluxes'.format(args.froot, '_txtl' if args.txtl else '', 'stacked' if args.resamp else 'flat')
    create_dataset(fname, X_train, X_test, obj_col, cols, y_vals, y_train, y_test)
