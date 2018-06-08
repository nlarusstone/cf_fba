from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import cf_io

# Compare reconstructed standard deviations to the original standard deviations
def plt_stds(rxns, stds, recon):
    plt.scatter(rxns, stds)
    plt.ylabel('Std deviation of fluxes')
    plt.xlabel('Reaction ID')
    plt.title('Standard deviation of {0}fluxes'.format('reconstructed ' if recon else ''))
    plt.show()

# Plots the correlation between our generated fluxes and the experimental dataset as we increase noise
def plt_noise_corr(noise_data, orig_corr, ndim=2):
    orig_enc, noise_data = noise_data[-1], noise_data[:-1]
    plt.figure(figsize=(10, 8))
    plt.title('Latent dimension = {0}'.format(ndim))
    plt.axhline(y=-1 * orig_enc[1] if orig_enc[1] < 0 else orig_enc[1], label='Original data encoded')
    plt.axhline(y=orig_corr, label='Original data correlation', c='g')
    for noise, corr in noise_data:
        plt.scatter(x=noise, y=-1 * corr if corr < 0 else corr)
        plt.xlabel('Noise amount')
        plt.ylabel('Correlation')
    plt.xscale('log')
    plt.legend()
    plt.show()

# Code to plot the latent space of a VAE or PCA
def plt_latent_space(encoded_data, df, ax, y_test=None, flat=False, dim_1=0, dim_2=1, samp_range=None,
                     color_scheme='tab20', sz=60, legend=True, label='Dimension'):
    if ax is None:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
        legend = True
    df_sor = df.sort_values(by='OUT')
    n_experiments = df.shape[0]
    ind = df_sor.index
    cmap = cm.get_cmap(color_scheme, n_experiments)

    if samp_range is not None:
        #low, high = samp_range
        if type(samp_range) == int:
            samp_range = range(samp_range)
        encoded_data = encoded_data[samp_range]
        if y_test is not None:
            y_test = y_test[samp_range]
    if flat:
        colors = y_test
        d_1 = encoded_data[:, dim_1]
        d_2 = encoded_data[:, dim_2]
    else:
        colors = np.array([[k] * encoded_data.shape[0] for k in range(n_experiments)]).T
        d_1 = encoded_data[:, ind, dim_1]
        d_2 = encoded_data[:, ind, dim_2]

    xmin, xmax = np.amin(d_1), np.amax(d_1)
    ymin, ymax = np.amin(d_2), np.amax(d_2)
    x_diff = (xmax - xmin) / 10.0
    y_diff = (ymax - ymin) / 10.0    
    print 'Plot d_1 vs d_2'
    sc = ax.scatter(d_1, d_2, c=colors, cmap=cmap, s=sz)
    ax.set_xlim((xmin - x_diff, xmax + x_diff))
    ax.set_ylim((ymin - y_diff, ymax + y_diff))
    if legend:
        lp = lambda k: ax.plot([],color=sc.cmap(sc.norm(k)), ms=np.sqrt(sz), mec="none",
                                label=' '.join([col + ': {0}'.format(df_sor.iloc[k, col_i]) for col_i, col in enumerate(df_sor.columns[:-1])]),
#'Sugar: {0}, Phosphate: {1}, NTs: {2}, Potassium: {3}'.format(
#                                    df_sor.iloc[k, 0], df_sor.iloc[k, 1], df_sor.iloc[k, 2], df_sor.iloc[k, 3]), 
                                ls="", marker="o")[0]
        handles = [lp(k) for k in np.unique(range(n_experiments))]
        ax.legend(bbox_to_anchor=(1.1, 1.0))#handles=handles, 
    ax.set_xlabel('{0} {1}'.format(label, dim_1))
    ax.set_ylabel('{0} {1}'.format(label, dim_2))

# Compares a flat VAE to a stacked VAE
def compare_flat_stacked(df, test_enc, test_enc_f, y_test_f):
    fig, axarr = plt.subplots(nrows=1, ncols=2, figsize=(14, 8))
    plt_latent_space(test_enc_f, df, axarr[0], flat=True, y_test=y_test_f, samp_range=(0, 1000), legend=False)
    axarr[0].set_title('VAE latent space')
    plt_latent_space(test_enc, df, axarr[1], samp_range=(0, 1000), legend=True)
    axarr[1].set_title('Corr-VAE latent space')
    plt.show()
