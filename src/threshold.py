import argparse
import numpy as np
import cf_io
import correlation as corr

parser = argparse.ArgumentParser(description='VAE on metabolic fluxes')
parser.add_argument('-d', '--dim', metavar='d', type=int, help='latent dimensionality for vae',
                    default=2)
parser.add_argument('--corr', dest='corr', help='Correlation loss', default=True, action='store_true')
parser.add_argument('--no-corr', dest='corr', help='Correlation loss', action='store_false')
parser.add_argument('-l','--layers', nargs='+', help='Layer sizes', type=int, default=[1024, 1024, 1024])
parser.add_argument('-r', '--resamp', dest='resamp', help='Resample', default=True, action='store_true')
parser.add_argument('--no-resamp', dest='resamp', help='Dont resample', action='store_false')
parser.add_argument('-t', '--txtl', dest='txtl', default=False, action='store_true')
parser.add_argument('--no-txtl', dest='txtl', action='store_false')
parser.add_argument('-f', '--froot', type=str, default='manual')
parser.add_argument('-s', '--dset', type=str, default='nls')
parser.add_argument('-m', '--model', type=str, default='iJO1366')
parser.add_argument('-n', '--epochs', type=int, default=200)


if __name__ == '__main__':
    args = parser.parse_args()
    froot = args.froot
    #dset = args.dset
    if froot == 'karim':
        dset = 'karim'
    else:
        dset = 'nls'
    txtl = args.txtl
    resamp = args.resamp
    latent_dim = args.dim
    layer_szs = args.layers
    use_corr = args.corr
    model = args.model

    encoder, generator, X_test, y_test, obj_col, cols, y_vals_d, test_enc, test_dec = cf_io.get_test_data(froot, txtl, resamp, latent_dim, layer_szs, use_corr=use_corr, n_epochs=args.epochs)
    cols = np.array(cols)
    df = cf_io.get_exp_data(froot)
    model_f = '../bio_models/{0}/{1}/{2}_cf{3}.sbml'.format(dset, froot, model, '_txtl' if txtl else '')
    res, bad_rxns, best = corr.optimized_corr(model_f, df, test_dec, cols, t_low=-2, t_high=-1, n_steps=20)
    result_f_base = '../results/{0}_{1}_{2}_reduced{3}_dim={4}_last_sz={5}'.format(dset, froot, model, '_txtl' if txtl else '', latent_dim, layer_szs[-1])
    res.to_csv(result_f_base + '_res_new')
    bad_rxns.to_csv(result_f_base + '_rxns_new')
    best.to_csv(result_f_base + '_best_new')
