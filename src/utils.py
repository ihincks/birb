import re
import sys
from random import shuffle
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import seaborn as sbn
import h5py
import pystan as ps
import os
import dill
import warnings
from six import StringIO
import copy
from tqdm import tqdm
from wurlitzer import pipes, STDOUT
import scipy.optimize as opt
import scipy.stats as st
from scipy.special import binom, gammaln, xlogy, expm1, log1p, logit, expit, psi

#===============================================================================
# FUNCTIONS
#===============================================================================

overwrite = True
def export_figure(fig, name, extensions=['.png', '.pdf']):
    for ext in extensions:
        do_save = False
        fname = '../fig/' + name + ext
        if not os.path.isfile(fname):
            print('Saving {}'.format(fname))
            do_save = True
        elif overwrite:
            print('Overwriting {}'.format(fname))
            do_save = True
        else:
            print('Skipping {}'.format(fname))
        if do_save:
            fig.savefig(fname)
            
def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = np.array(colorsys.rgb_to_hls(*mc.to_rgb(c)))
    return colorsys.hls_to_rgb(c[0],1-amount * (1-c[1]),c[2])

def advi(model, data, **kwargs):
    # Dumb pystan saves advi output to a csv file, so we have
    # write a bunch of code to go fetch and parse it
    for tries in range(5):
        try:
            output = model.vb(data=data, **kwargs)
            break
        except Exception as e:
            if tries == 4:
                raise e
    raw = np.genfromtxt(output['args']['sample_file'], skip_header=7, delimiter=',')
    with open(output['args']['sample_file']) as f:
        for i in range(5):
            names = f.readline()[:-1].split(',')
    #os.remove(output['args']['sample_file'])
    #os.rmdir(os.path.dirname(output['args']['sample_file']))
    def splitname(name):
        x = name.split('.')
        return [False, name, 0] if len(x) ==1 else [True, x[0], len(x)-1]
    scalar_names = []
    vector_names = []
    for is_vector, name, dim in map(splitname, names):
        if is_vector:
            vector_names += [(name, dim)]
        else:
            scalar_names += [name]
    vector_names = list(set(vector_names))
    output = {}
    for scalar in scalar_names:
        output[scalar] = raw[:,names.index(scalar)]
    for vector, dim in vector_names:
        idxs = []
        shapes = [1] * dim
        for idx, name in enumerate(names):
            if re.search(vector + r'\..+', name):
                idxs += [idx]
                for k, d in enumerate(map(int, name.split('.')[1:])):
                    shapes[k] = max(shapes[k], d)
        mat = raw[:,idxs]
        output[vector] = mat.reshape([mat.shape[0]] + shapes)
    return output

def _extract_flatname(fit, flatname, inc_warmup=False, join_chains=True):
    # parts = re.match(r'([^\[]+)(\[(.+)\])?', flatname).groups()
    # name = parts[0]
    # if len(parts) == 1:
    #     out =  fit.extract(name)[name]
    # else:
    #     dim_idxs = tuple(map(int, parts[-1].split(',')))
    #     out = fit.extract(name)[name][(slice(None),) + dim_idxs]
    # if out.ndim > 1:
    #     raise ValueError('Given name {} is not a flatname.'.format(flatname))
    # return out
    found = False
    for idx, name in enumerate(fit.flatnames):
        if name == flatname:
            found = True
            break
    if not found:
         raise ValueError('Given name {} is not a flatname.'.format(flatname))
    out = fit.extract(inc_warmup=inc_warmup)[...,idx]
    return out.flatten() if join_chains else out

def pairs(fit, flatnames=None, max_vars=10, shuffle_vars=True, panel_size=3.5, inc_warmup=True, diag_format='hist',label_loc=[0.05,0.9], true_vals={}, show_divergent_points=True):
    flatnames = fit.flatnames if flatnames is None else flatnames
    if max_vars is not None:
        if len(flatnames) > max_vars:
            if shuffle_vars:
                shuffle(flatnames)
            flatnames = flatnames[:max_vars]
    n = len(flatnames)

    if inc_warmup:
        accept_stat = np.concatenate([param['accept_stat__'] for param in fit.get_sampler_params()])
        div_mask = np.concatenate([param['divergent__'] for param in fit.get_sampler_params()])
    else:
        warmup = fit.sim['warmup']
        accept_stat = np.concatenate([param['accept_stat__'][warmup:] for param in fit.get_sampler_params()])
        div_mask = np.concatenate([param['divergent__'][warmup:] for param in fit.get_sampler_params()])
    med = np.median(accept_stat)
    upper_mask = accept_stat >= med
    lower_mask = accept_stat < med
    div_mask = div_mask.astype('bool')

    fig, axarr = plt.subplots(n, n, figsize=(panel_size * n,) * 2)
    for idx_row in range(n):
        for idx_col in range(n):
            if idx_row == idx_col:
                name = flatnames[idx_row]
                data = _extract_flatname(fit, name, inc_warmup=inc_warmup)
                ax = axarr[idx_row, idx_col]
                if diag_format == 'hist':
                    ax.hist(data, bins=40)
                elif diag_format == 'kde':
                    sbn.kdeplot(data, ax=ax)
                ax.text(label_loc[0], label_loc[1], '{}'.format(name),
                    transform = ax.transAxes, va='center', ha='left',
                    bbox={'facecolor':'white', 'alpha':0.5, 'pad':4, 'lw':1}
                )
                if name in true_vals:
                    ax.axvline(x=true_vals[name], color='r', lw=1.5)
            if idx_row < idx_col:
                namex, namey = flatnames[idx_row], flatnames[idx_col]
                datax = _extract_flatname(fit, namex, inc_warmup=inc_warmup)
                datay = _extract_flatname(fit, namey, inc_warmup=inc_warmup)
                for ax, mask in [[axarr[idx_row, idx_col], upper_mask],[axarr[idx_col, idx_row], lower_mask]]:
                    ax.plot(datax[mask], datay[mask], '*')
                    if show_divergent_points:
                        ax.plot(datax[div_mask], datay[div_mask], 'r*')
                    ax.text(label_loc[0], label_loc[1], '({},{})'.format(namex, namey),
                        transform = ax.transAxes, va='center', ha='left',
                        bbox={'facecolor':'white', 'alpha':0.5, 'pad':4, 'lw':1}
                    )
    return fig
    
def rb_wlsf(data, suffix='', guesses=np.array([1,1,0.5])):
    """
    Given the dictionary data containing Q, Nbin, and m, does
    a weighted least-squares fit to estimate (A-B) * p**m + B.
    The weights used for each sequence length are as prescribed 
    by Robin-Blume Kohout, where only the Gaussian case is implemented, and
    an error is returned if it is detected that we are not in the Gaussian case.
    """
    
    datavals = data['Q{}'.format(suffix)].astype('float')
    Nbin = data['Nbin']
    Nsamp = data['Nsamp']
    mlist = data['m']
        
    # give an error if we are in a bad noise regime
    if np.any(np.var(datavals, axis=1) < 10**-8):
        warnings.warn('Supplied data has a zero-variance sequence length and this method will not be reliable')
    if Nbin < 10:
        warnings.warn('It is unwise to do weighted least squares fitting with such a low Nbin')
    
    # compute a couple types of moments recommended by RBK
    f_hat = np.sum(datavals, axis=1) / (Nbin * Nsamp)
    delta_f_tight = f_hat * (f_hat - 1) / Nsamp
    delta_f_empirical = np.sum(datavals**2, axis=1) / (Nsamp * Nbin**2) - f_hat**2
    delta_f = np.sqrt(np.max(np.vstack([delta_f_tight, delta_f_empirical]), axis=0) / Nsamp)
    delta_f_tight
    delta_f_empirical
    
    # non-linear model with explanitory variable m
    def rb(m,p,A,B):
        return (A-B)*(p**m) + B
    
    # do non-linear least-squares fit
    pABopt, pABcov = curve_fit(
        rb, 
        mlist, 
        f_hat, 
        p0=guesses, 
        sigma=delta_f, 
        absolute_sigma=True
    )
    
    # return estimate, covariance matrix of estimate, and weights used at each sequence length
    return pABopt, pABcov, delta_f
    
def rb_lsf(data, suffix='', guesses=np.array([1,1,0.5])):
    """
    Given the dictionary data containing Q, Nbin, and m, does
    a least-squares fit to estimate (A-B) * p**m + B.
    """
    
    datavals = data['Q{}'.format(suffix)].astype('float')
    Nbin = data['Nbin']
    Nsamp = data['Nsamp']
    mlist = data['m']
    
    # compute mean at each sequence length
    f_hat = np.sum(datavals, axis=1) / (Nbin * Nsamp)
    
    # non-linear model with explanitory variable m
    def rb(m,p,A,B):
        return (A-B)*(p**m) + B
    
    # do non-linear least-squares fit
    pABopt, pABcov = curve_fit(
        rb, 
        mlist, 
        f_hat, 
        p0=guesses, 
        absolute_sigma=True
    )

    # return estimate, covariance matrix of estimate
    return pABopt, pABcov, None
    
def rb_naive_wlsf(data, suffix='', guesses=np.array([1,1,0.5])):
    """
    Given the dictionary data containing Q, Nbin, and m, does
    a weighted least-squares fit to estimate (A-B) * p**m + B.
    """
    
    datavals = data['Q{}'.format(suffix)].astype('float')
    Nbin = data['Nbin']
    Nsamp = data['Nsamp']
    mlist = data['m']
    
    # compute a couple types of moments recommended by RBK
    f_hat = np.sum(datavals, axis=1) / (Nbin * Nsamp)
    delta_f = np.sqrt(np.var(datavals / Nbin, axis=1))
    
    # give an error if we are in a bad noise regime
    if np.any(np.var(datavals, axis=1) < 10**-8):
        warnings.warn('Supplied data has a zero-variance sequence length and this method will not be reliable')
    if Nbin < 10:
        warnings.warn('It is unwise to do weighted least squares fitting with such a low Nbin')
    
    # non-linear model with explanitory variable m
    def rb(m,p,A,B):
        return (A-B)*(p**m) + B
    
    # do non-linear least-squares fit
    pABopt, pABcov = curve_fit(
        rb, 
        mlist, 
        f_hat, 
        p0=guesses, 
        sigma=delta_f, 
        absolute_sigma=True
    )
    
    # return estimate, covariance matrix of estimate, and weights used at each sequence length
    return pABopt, pABcov, delta_f
    
def rb_lsf_completion(data, suffix='', guesses=np.array([1,1,0.5]), method=rb_wlsf):
    """
    Given the dictionary data containing Q, Nbin, and m, does
    a least-squares fit to estimate (A-B) * p**m + B, adding
    p_est, A_est, B_est, p_std, A_std, B_std to a copy
    of the dictionary and returning it. The guess is generated by drawing a 
    sample from the WLF estimate's multivariate-normal distrubution, and 
    post-selecting on it being in-bounds.
    """
    
    pABopt, pABcov, aoeu = method(data, suffix=suffix, guesses=guesses)
    
    guess = pABopt
    while np.any(np.clip(guess, 0, 1) != guess):
        # if the point estimate is out of bounds, sample from the truncated normal
        guess = np.random.multivariate_normal(pABopt, pABcov)
    if aoeu is None:
        ts = 0.1 * np.ones(3)
    else:
        ts = np.diag(pABcov) / (guess * (1 - guess))**2
    
    estimates = data.copy()
    estimates.update({
        'p_est': guess[0], 'A{}_est'.format(suffix): guess[1], 'B{}_est'.format(suffix): guess[2], 'B_est': guess[2],
        'p_t': ts[0], 'A{}_t'.format(suffix): ts[1], 'B{}_t'.format(suffix): ts[2], 'B_t': ts[2]
    })
    
    return estimates
    
def save_stan_fit(name, fit):
    """
    Pickles the given StanFit4Model to file in ../data/fits/name.fit
    """
    folder = os.path.join('..', 'data', 'fits')
    if not os.path.exists(folder):
        os.makedirs(folder)
    storage_file = os.path.join(folder , name + '.fit')
    with open(storage_file, 'wb') as f: 
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            dill.dump(fit, f)
        
def load_stan_fit(name):
    """
    Tries to unpickle the given StanFit4Model from ../data/fits/name.fit
    """
    folder = os.path.join('..', 'data', 'fits')
    storage_file = os.path.join(folder , name + '.fit')
    try:
        with open(storage_file, 'rb') as f:
            fit = dill.load(f)
        return fit
    except IOError as e:
        # Downgrade to warning
        warnings.warn('{}'.format(e))
        return None
        
#===============================================================================
# MLE and bootstrapping
#===============================================================================

def psi_plus_k(z, k):
    """
    Computes the digamma function psi(z+k) returning both
    psi(z) and psi(z+k), where k is an array of integers.
    Exploits the relationship psi(z+1)=psi(z)+1/z.
    This should be much faster than psi(z+k) when max(k) is 
    not too big; do your own timing tests for your use case.
    """
    # begin by computing psi(z+i) for all i=0,...,max(k)
    # we are computing more than we need, but in our use case, k will
    # have many values in this range, and we also get the benefit of cumsum.
    max_k = np.max(k)
    harmonic_sum = np.cumsum(1 / (z[...,np.newaxis] + np.arange(0, max_k)), axis=-1)
    harmonic_sum = np.concatenate(
        [np.zeros(z.shape + (1,)), harmonic_sum], 
        axis=-1
    ).reshape(-1, max_k+1)
    
    # now it is a matter of correctly selecting the values we need
    psi_z = psi(z)
    harmonic_values = harmonic_sum[:,k.flatten().astype(np.int)]
    return psi_z, psi_z + harmonic_values.reshape(np.broadcast(z, k).shape)
    
def harmonic_sum(z, k):
    """
    Returns 1/z + 1/(z+1) + ... + 1/(z+i) for each i in k.
    Equivalent to psi(z+k) - psi(z).
    """
    # begin by computing psi(z+i) for all i=0,...,max(k)
    # we are computing more than we need, but in our use case, k will
    # have many values in this range, and we also get the benefit of cumsum.
    max_k = int(np.max(k))
    harmonic_sum = np.cumsum(1 / (z[...,np.newaxis] + np.arange(0, max_k)), axis=-1)
    harmonic_sum = np.concatenate(
        [np.zeros(z.shape + (1,)), harmonic_sum], 
        axis=-1
    ).reshape(-1, max_k+1)
    
    # now it is a matter of correctly selecting the values we need
    return harmonic_sum[np.arange(len(harmonic_sum))[:,np.newaxis], k.astype(np.int)]

def log_betabinom_r(k, a, b, n):
    return gammaln(n + 1) + gammaln(k + a) + gammaln(n - k + b) + gammaln(a + b) \
        - gammaln(k + 1) - gammaln(n - k + 1) - gammaln(n + a + b) - gammaln(a) - gammaln(b)
def grad_log_betabinom(k, a, b, n):
    s = a + b
    common = psi(s) - psi(s + n)
    da = np.sum(common + harmonic_sum(a, k), axis=1)
    db = np.sum(common + harmonic_sum(b, n-k), axis=1)
    # da = np.sum(common - psi(a) + psi(a + k), axis=1)
    # db = np.sum(common - psi(b) + psi(b + n - k), axis=1)
    return np.concatenate([da[...,np.newaxis], db[...,np.newaxis]], axis=1)
def jac_mur_ab(mu, v, transform):
    jac = np.empty(mu.shape + (2,2))
    if transform == 'r':
        jac[...,0,0] = 1 / (v * (1 - mu)**2) - 1
        jac[...,0,1] = 1 - 1/(v * mu**2)
        jac[...,1,0] = 1 / (v**2 * (mu - 1))
        jac[...,1,1] = - 1 / (v**2 * mu)
    elif transform == 't':
        jac[...,0,0] = 1 / v - 1
        jac[...,0,1] = -jac[...,0,0]
        jac[...,1,0] = -mu / v**2
        jac[...,1,1] = (mu - 1) / v**2
    return jac
def jac_ABp_mu(A, B, p, seq_lengths):
    pm = p ** seq_lengths
    jac = np.empty((seq_lengths.size, 3,))
    jac[:,0] = pm
    jac[:,1] = 1 - pm
    jac[:,2] = (seq_lengths * (A - B)) * p ** seq_lengths
    return jac
def log_likelihood(x, data, transform='t', compute_jac=True):
    ex = expit(x)
    A = ex[0]
    B = ex[1]
    p = ex[2]
    vs = ex[3:]
    
    n_samp = float(data['Nsamp'])
    n_bin = float(data['Nbin'])
    Q = data['Q'].astype('float64') # shape (n_seq_lengths, n_samp)
    seq_lengths = data['m']
    n_seq_lengths = seq_lengths.size
    
    # mus, rs, a, b will all have shape (n_seq_lengths,)
    mus = (A - B)* p ** seq_lengths[:] + B
    mus = (A - B)* p ** seq_lengths[:] + B
    if transform == 'r':
        a = 1 / (vs * (1 - mus)) - mus
        b = 1 / (vs * mus) - 1 + mus
    elif transform == 't':
        a = mus * (1 / vs - 1)
        b = (1 - mus) * (1 - vs) / vs
        
    log_lik = np.sum(
        log_betabinom_r(Q, a[:,np.newaxis], b[:,np.newaxis], n_bin)
    )
    
    if compute_jac:
        # we think of the objective function as four composed functions, whose
        # jacobians we must multiply. several of the functions have 
        # diagonal or block jacobians, and so we store more compact forms.
        
        # jacobian matrix of the first two compositions
        jac1_ABp = jac_ABp_mu(A, B, p, seq_lengths) * ex[np.newaxis,:3] * (1 - ex[np.newaxis,:3])
        jac1_r = ex[3:] * (1 - ex[3:])
        
        # jacobian matrix of the last two compusitions
        jac2 = np.matmul(
                jac_mur_ab(mus, vs, transform),
                grad_log_betabinom(Q, a[:,np.newaxis], b[:,np.newaxis], n_bin)[:,:,np.newaxis]
            )
        
        # compose the above two
        jac = np.empty((3 + n_seq_lengths, ))
        jac[:3] = np.sum(jac1_ABp * jac2[:,0,:], axis=0)
        jac[3:] = jac1_r * jac2[:,1,0]
        
        return log_lik, jac
    else:
        return log_lik
        
def grad_log_likelihood(x, data, transform='t', compute_jac=True):
    ex = expit(x)
    A = ex[0]
    B = ex[1]
    p = ex[2]
    vs = ex[3:]        
    
    n_samp = float(data['Nsamp'])
    n_bin = float(data['Nbin'])
    Q = data['Q'].astype('float64') # shape (n_seq_lengths, n_samp)
    seq_lengths = data['m']
    n_seq_lengths = seq_lengths.size
    
    # mus, rs, a, b will all have shape (n_seq_lengths,)
    mus = (A - B)* p ** seq_lengths[:] + B
    if transform == 'r':
        a = 1 / (vs * (1 - mus)) - mus
        b = 1 / (vs * mus) - 1 + mus
    elif transform == 't':
        a = mus * (1 / vs - 1)
        b = (1 - mus) * (1 - vs) / vs

    # we think of the objective function as four composed functions, whose
    # jacobians we must multiply. several of the functions have 
    # diagonal or block jacobians, and so we store more compact forms.
    
    # jacobian matrix of the first two compositions
    jac1_ABp = jac_ABp_mu(A, B, p, seq_lengths) * ex[np.newaxis,:3] * (1 - ex[np.newaxis,:3])
    jac1_r = ex[3:] * (1 - ex[3:])
    
    # jacobian matrix of the last two compusitions
    jac2 = np.matmul(
            jac_mur_ab(mus, vs, transform),
            grad_log_betabinom(Q, a[:,np.newaxis], b[:,np.newaxis], n_bin)[:,:,np.newaxis]
        )
    
    # compose the above two
    jac = np.empty((3 + n_seq_lengths, ))
    jac[:3] = np.sum(jac1_ABp * jac2[:,0,:], axis=0)
    jac[3:] = jac1_r * jac2[:,1,0]
    
    return -jac

def flatten_params(A, B, p, rs):
    return np.concatenate([
            np.atleast_1d(logit(A)),
            np.atleast_1d(logit(B)),
            np.atleast_1d(logit(p)),
            logit(rs).flatten()
    ])
def unflatten_params(x, data):
    n_seq_lengths = data['Nm']
    A = expit(x[0])
    B = expit(x[1])
    p = expit(x[2])
    rs = expit(x[3:].reshape(1, n_seq_lengths))
    return A, B, p, rs
def guess_params(A, B, p, data, random=True):
    n_seq_lengths = data['Nm']
    if random:
        rs = np.random.random((1, n_seq_lengths))
    else:
        rs = 0.5 * np.ones((1, n_seq_lengths))
    return flatten_params(A, B, p, rs)
def objective_function(x, data, transform='t', compute_jac=True):
    result = log_likelihood(x, data, transform=transform, compute_jac=compute_jac)
    if compute_jac:
        log_lik, jac = result
        return -log_lik, -jac
    else:
        return -result
        
def mle(data, guess=None, transform='t'):
    """
    Computes the MLE of the RB parameter p 
    with the given model, using a likelihood function
    which is beta-binomial distributed at every
    RB sequence. The parameters of this liklihood function
    are A,B,p and r_1,...,r_M for each sequence length.
    The r's model the variance of each beta distribution.
    """
    if guess is None:
        guess = guess_params(0.9,0.5,0.99,data)
    res = opt.minimize(
        objective_function, 
        guess,
        args=(data,transform,True),
        method='BFGS',
        jac=True
    )
    return unflatten_params(res.x, data)
    
def bootstrap_nonparam_sample(data, precomputed_mle=None):
    datavals = data['Q']
    n_samp = data['Nsamp']
    n_seq_lengths = data['Nm']
    choices = np.random.randint(0, n_samp, size=(n_seq_lengths,n_samp))
    new_datavals = np.empty(datavals.shape)
    for idx, choice_row in enumerate(choices):
        new_datavals[idx, :] = datavals[idx, choice_row]
    new_data = data.copy()
    new_data['Q'] = new_datavals
    return new_data
def bootstrap_param_sample(data, precomputed_mle=None):
    A,B,p,rs = mle(data) if precomputed_mle is None else precomputed_mle
    n_samp = data['Nsamp']
    n_bin = data['Nbin']
    seq_lengths = data['m']
    mus = (A - B)* p ** seq_lengths[:] + B
    a = 1 / (rs * (1 - mus)) - mus
    b = 1 / (rs * mus) - 1 + mus
    a = np.repeat(a, n_samp)
    b = np.repeat(b, n_samp)
    qs = st.beta.rvs(a=a,b=b)
    if not np.all(np.isfinite(qs)):
        # in case something weird happens like p=1
        return bootstrap_param_sample(data)
    new_data = data.copy()
    new_data['Q'] = st.binom.rvs(n=n_bin, p=qs).reshape((-1,n_samp))
    return new_data
def bootstrap_ABp_mle(data, n, sampler=bootstrap_param_sample, transform='t'):
    ps = np.empty((n, 3))
    x0 = flatten_params(*mle(data, transform=transform))
    for idx in range(n):
        result = mle(sampler(data), guess=x0, transform=transform)
        A, B, p = result[:3]
        ps[idx,:] = [A,B,p]
    return ps
def bootstrap_ABpv_mle(data, n, sampler=bootstrap_param_sample, transform='t'):
    out = np.empty((n, 3 + data['Nm']))
    data_mle = mle(data, transform=transform)
    x0 = flatten_params(*data_mle)
    for idx in range(n):
        A,B,p,vs = mle(sampler(data, precomputed_mle=data_mle), guess=x0, transform=transform)
        out[idx,:] = np.concatenate([np.array([A,B,p]), vs.flatten()])
    return out

#===============================================================================
# CLASSES
#===============================================================================

class RBSimData(object):
    def __init__(self, filename):
        self._filename = filename
        self.reload()

    def reload(self):
        h5f = h5py.File(self._filename, 'r')
        
        self.protocol_name = h5f['ProtocolName'].value
        self.gateset_name = h5f['GateSetName'].value
        
        self.sequence_lengths = h5f['SequenceLengths'].value
        self.experiment_types = h5f['ExperimentTypes'].value
        
        self.dimension = h5f['Dimension'].value
        self.size = h5f['Size'].value
        
        self.survival_data = h5f['SurvivalData'].value
        
        h5f.close()

        
class StanModelFactory(object):
    """
    Class to construct instances of pystan.StanModel, which first checks
    if the model has already been compiled and saved to disk, loading it 
    if it has.
    
    :param str filename: filename of stan code to load
    :param storage_folder: list of strings specifying storage folder
    """
    STORAGE_FOLDER = ['..', 'data']
    def __init__(self, filename, storage_folder=STORAGE_FOLDER):
        self._filename = filename
        storage_filename = os.path.splitext(os.path.basename(filename))[0] + '.pkl'
        self._storage_name = os.path.join(os.path.join(*storage_folder), storage_filename)
        self._model = None
        
    def _load_model_from_disk(self):
        """
        Tries to load a pickled StanModel object from _storage_name. Returns 
        this, or None if it fails.
        """
        try:
            with open(self._storage_name, 'rb') as f:
                model = dill.load(f)
        except IOError:
            model = None
        return model
        
        
    def _save_model_to_disk(self, model):
        """
        Pickles the given model and saves it to file.
        
        :param model: StanModel object to pickle and save to _storage_name
        """
        with open(self._storage_name, 'wb') as f: 
            dill.dump(model, f)
                
    def _get_model_code(self):
        """
        Reads _filename and returns its contents as a string.
        """
        with open(self._filename, 'r') as f:
            model_code = "".join(f.read())
        return model_code
        
    def _up_to_date(self):
        """
        Decides if _model is up-to-date. Returns True if _model exists and 
        has model_code equal to the current contents of _filename, False 
        otherwise.
        """
        if self._model is None:
            return False
        else:
            return self._model.model_code == self._get_model_code()
                
    def _get_model(self):
        """
        Loads and unpickles the StanModel from disk if it exists, and Returns
        it if it is up-to-date. Otherwise, compiles a new StanModel.
        """
        model = self._load_model_from_disk()
        if model is not None and model.model_code == self._get_model_code():
            return model

        model = ps.StanModel(self._filename)
        self._save_model_to_disk(model)
        return model
        
    @property
    def model(self):
        """
        A StanModel instance of the model code located at the filename given
        at construction, but up to date with the current contents of the file.
        """
        if not self._up_to_date():
            self._model = self._get_model()
        return self._model

class TqdmUpTo(tqdm):
    """Provides `update_to(n)` which uses `tqdm.update(delta_n)`."""
    def update_to(self, b=1, bsize=1, tsize=None):
        """
        b  : int, optional
            Number of blocks transferred so far [default: 1].
        bsize  : int, optional
            Size of each block (in tqdm units) [default: 1].
        tsize  : int, optional
            Total size (in tqdm units). If [default: None] remains unchanged.
        """
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)  # will also set self.n = b * bsize
        
class CallbackStringIO(StringIO):
    def __init__(self, callback, initial_value='', newline='\n'):
        self._callback=callback
        super(CallbackStringIO, self).__init__(
            initial_value=initial_value, 
            newline=newline
        )
    def write(self, s):
        self._callback(s)
        super(CallbackStringIO, self).write(s)
        
class StanModelWrapper(ps.StanModel):
    def __init__(self, model):
        self._model = model

    def get_cpp_code(self):
        return self._model.get_cpp_code()
    def get_cxxflags(self):
        return self._model.get_cxxflags()

    def sampling(self, *args, **kwargs):
        def print_stdout(s):
            print s
        out = CallbackStringIO(print_stdout)
        with pipes(stdout=out, stderr=STDOUT):
            result = self._model.sampling(*args, **kwargs)
        print out.getvalue()
        return result
        
