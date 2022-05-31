# -*- coding: utf-8 -*-
"""
@author: Sanjay Kottapalli (svk4001@med.cornell.edu)
@date: 05/31/2022
Landau Lab, WCM

This script executes functions for the genotyping of cells, 
starting from WT and MUT read counts.
"""
#imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
#import sklearn
import timeit

from sklearn.neighbors import KernelDensity
from numpy.polynomial.polynomial import Polynomial
from numpy.polynomial.polynomial import polyfit
from scipy.optimize import minimize
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
from scipy.integrate import quad
from collections import Counter
from sklearn.metrics.pairwise import euclidean_distances
from scipy.stats import gaussian_kde
from sklearn.metrics import adjusted_rand_score
from sklearn.cluster import SpectralClustering
from sklearn.semi_supervised import SelfTrainingClassifier
from sklearn.neighbors import KNeighborsClassifier

import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 500
sns.set(font_scale=1.2)
sns.set_style("ticks")
tab20=plt.get_cmap('tab20')
gen_cmap = {'MUT':tab20.colors[0],'WT':tab20.colors[2],
            'HET':tab20.colors[4],'NA':tab20.colors[6],'-1':tab20.colors[8]}

plt.clf()

def GotchaLabeling(path="", infile="", gene_id="", sample_id=""):
    time1 = timeit.default_timer()
    print("Reading in file.")
    typing = read_data(path+infile, gene_id, sample_id)
    sample_dir = path+sample_id+'/'
    try:
        os.mkdir(sample_dir)
    except:
        pass
    
    for i in ["WT", "MUT"]:
        print("Noise correcting {} read counts.".format(i))
        if i=="WT":
            typing, wt_min = noise_correct(typing,i,sample_dir)
        else:
            typing, mut_min = noise_correct(typing,i,sample_dir)
    
    print("Performing quadrant genotyping.")
    typing = quadrant_genotype(typing, wt_min, mut_min)
    
    print("Fitting KDE on whole dataset.")
    typing, data_smooth, bw = KDE_2D(typing, wt_min, mut_min)
    
    print("Choosing optimal number of clusters.")
    typing, n_clusters, opt_clusters, affinity = \
        cluster_number(typing, data_smooth, bw, wt_min, mut_min)
        
    print("Computing KDE mixture model-based clusters.")
    typing = KDE_mixture_model(typing, bw, wt_min, mut_min, 
                          data_smooth, affinity, 
                          opt_clusters, sample_dir)
    
    typing.to_csv(sample_dir+sample_id+'_genotype_labels.csv')
    print("All analysis complete!")
    time2 = timeit.default_timer()
    print("Total time to execute: {}".format(time2-time1))
    
    return typing

def read_data(infile="", gene_id="", sample_id=""):
    '''
    This function reads in the sequencing reads from only real cells
    (based on ATAC calling). Cells without GoTChA reads are dropped.
    '''
    cell_line = pd.read_csv(infile, index_col=0)
    genotyping = pd.DataFrame(index=cell_line.index)
    genotyping['WTcount'] = cell_line[gene_id+'_WTcount']
    genotyping['MUTcount'] = cell_line[gene_id+'_MUTcount']
    index = genotyping[['WTcount', 'MUTcount']].dropna().index
    genotyping = genotyping.loc[index,:]

    return genotyping

def noise_correct(typing, feature="", sample_dir=""):
    #add noise
    np.random.seed(0)
    pseudocount = 1
    X = typing[[feature+'count']]+pseudocount
    logged_counts = np.log(X) #change to log10 later, arctanh?
    
    noise = np.random.normal(0,0.3,typing.shape[0])
    logged_counts = logged_counts + noise.reshape(-1,1)
    typing['transf_{}'.format(feature)] = logged_counts
    
    plt.show()
    plt.hist(np.log(X), density=True, bins=50)
    plt.title(feature+' counts')
    plt.ylabel("Probability")
    plt.xlabel("Log(counts+1)")
    plt.savefig(sample_dir+"log_{}_counts.pdf".format(feature), 
                dpi=500, bbox_inches = "tight")
    plt.show()
    plt.clf()
        
    #bandwidth optimization
    print("Computing optimal bandwidth...")
    time1 = timeit.default_timer()
    data = typing[['transf_{}'.format(feature)]].values
    
    bw_list = []
    results_list = []
    for i in np.arange(0.1,1.0,0.025):
        bw_list.append(i)
        res = optimize_bandwidth([i], data)
        results_list.append(res)
    
    fit = polyfit(bw_list, results_list, 10)
    poly = Polynomial(fit)
    
    poly_x = np.linspace(0.1,1.0,1000)
    poly_y = poly(poly_x)
    plt.scatter(bw_list, results_list)#bw_list, results_list)
    plt.plot(poly_x, poly_y)
    plt.ylabel("Loss")
    plt.xlabel("Bandwidth")
    plt.savefig(sample_dir+"{}_bandwidth.pdf".format(feature), 
                dpi=500, bbox_inches = "tight")
    plt.show()
    plt.clf()
    
    bw_guess = bw_list[np.argmin(results_list)]
    result = minimize(poly, x0=np.array([bw_guess]), 
                      options={'disp':True}, bounds=((0.1, 0.975),))
    bw = result.x[0]
    
    time2 = timeit.default_timer()
    print("Optimal bandwidth computed in {} seconds.".format(
        time2-time1))
    
    #plot initial KDE
    kde = KernelDensity(bandwidth=bw)#kernel='exponential')
    kde.fit(typing['transf_{}'.format(feature)].values.reshape(-1, 1))
    x_bin = np.histogram(typing['transf_{}'.format(feature)], bins=50)[1]
    kde_x = np.linspace(min(x_bin),max(x_bin),1000)
    kde_smooth = np.exp(kde.score_samples(kde_x.reshape(-1, 1)))
    
    plt.hist(typing['transf_{}'.format(feature)], density=True, bins=50)
    plt.plot(kde_x, kde_smooth, color='red')
    plt.title('Initial KDE')
    plt.ylabel("Probability")
    plt.xlabel("Log(counts+1)")
    plt.savefig(sample_dir+"{}_kde_initial.pdf".format(feature), 
                dpi=500, bbox_inches = "tight")
    plt.show()
    plt.clf()
    
    #plot smoothed KDE
    kde_test = np.exp(kde.score_samples(kde_x.reshape(-1, 1)))
    kde_smooth = savgol_filter(kde_test, 333, 3)
    
    plt.hist(typing['transf_{}'.format(feature)], density=True, bins=50)
    plt.plot(kde_x, kde_smooth, color='red')
    plt.title('Smoothed KDE')
    plt.ylabel("Probability")
    plt.xlabel("Log(counts+1)")
    plt.savefig(sample_dir+"{}_kde_smoothed.pdf".format(feature), 
                dpi=500, bbox_inches = "tight")
    plt.show()
    plt.clf()
    
    #find minima
    '''
    CHANGE THIS PART TO ONLY DEPEND ON MINIMUM (linear combination)
    '''
    peaks, properties = find_peaks(kde_smooth)
    #top_2 = np.argsort(kde_smooth[peaks])[-2:]
    min_pdf, properties = find_peaks(-1*kde_smooth)
    min_pdf = [i for i in kde_x[min_pdf] if i>kde_x[peaks][0] 
               and i<kde_x[peaks][-1]]
    min_kde = np.exp(kde.score_samples(np.array(min_pdf).reshape(-1, 1)))
    arg_min = np.argmin(min_kde)
    min_pdf = min_pdf[arg_min]
    center_list = np.sort(kde_x[peaks][kde_x[peaks]>min_pdf])
    center = center_list[-1]
    #width = center-min_pdf
    
    kde_signal = KernelDensity(bandwidth=bw)
    signal_values = logged_counts
    signal_values = signal_values[signal_values>min_pdf]
    signal_values = signal_values.dropna()
    kde_signal.fit(signal_values.values.reshape(-1, 1))
    kde_mean2 = np.exp(kde_signal.score_samples(kde_x.reshape(-1, 1)))
    
    #plot mixture KDE
    plt.hist(typing['transf_{}'.format(feature)], density=True, bins=50)
    plt.plot(kde_x, kde_test, color='red')
    plt.plot(kde_x, 
             (kde_mean2/max(kde_mean2))*np.exp(kde.score(center.reshape(1,-1))), color='pink')
    plt.plot(kde_x, 
             kde_test-(kde_mean2/max(kde_mean2))*np.exp(kde.score(center.reshape(1,-1))), color='yellow')
    plt.title('Noise vs. Signal')
    plt.ylabel("Probability")
    plt.xlabel("Log(counts+1)")
    plt.xticks(np.arange(0,11,1.0))
    plt.savefig(sample_dir+"{}_kde_mixture.pdf".format(feature), 
                dpi=500, bbox_inches = "tight")
    plt.show()
    plt.clf()
    
    def subtract(x):
        total = np.exp(kde.score_samples(np.array(x).reshape(1, -1)))[0]
        signal = np.exp(kde_signal.score_samples(np.array(x).reshape(1, -1)))[0]
        signal = (signal/max(kde_mean2))*np.exp(kde.score(center.reshape(1,-1)))
        diff = total-signal
        if diff >= 0:
            return diff
        else:
            return 0.0
    
    factor = quad(func=subtract, a=-np.inf, b=np.inf)[0]
        
    def diff_pdf(x):
        return subtract(x)/factor
    
    def mean_func(x):
        return x*diff_pdf(x)
    
    def var_func(x):
        return x*x*diff_pdf(x)

    noise_mean = quad(func=mean_func, a=-np.inf, b=np.inf, limit=100)[0]
    noise_var = quad(func=var_func, a=-np.inf, 
                     b=np.inf, limit=100)[0] - noise_mean**2
    print("Noise mean, noise var:")
    print(noise_mean, noise_var)
    
    noise_std = np.sqrt(noise_var)
    typing['transf_{}'.format(feature)] = (logged_counts-noise_mean)/noise_std
    
    #added_mean = 0
    added_var = 0.3**2
    new_var = noise_var - added_var
    new_std = np.sqrt(new_var)
    
    typing['transf_{}'.format(feature)] = (np.log(X)-noise_mean)/new_std
    
    plt.hist(typing['transf_{}'.format(feature)], bins=50, density=True)
    plt.title('{} Z-scores'.format(feature))
    plt.ylabel("Probability")
    plt.xlabel("Z-score")
    plt.savefig(sample_dir+"{}_zscores.pdf".format(feature), 
                dpi=500, bbox_inches = "tight")
    plt.show()
    plt.clf()
    
    feat_min = (min_pdf-noise_mean)/new_std
    if feat_min<2.0:
        feat_min=2.0
    print("PDF minimum: "+str(feat_min))
    
    return typing, feat_min
        
def optimize_bandwidth(x, data_input):
    #print(x)
    kde = KernelDensity(bandwidth=x[0])
    
    n = data_input.shape[0]
    n_splits=10
    score_list = []
    
    for i in range(10):
        test_index = np.random.choice(np.arange(n), 
                                      size=round(n/n_splits))
        data_temp = data_input[test_index,:]
        #print(data_temp.shape)
        #print(data_temp)
        kde.fit(data_temp)
        loglike = kde.score(data_input)
        score_list.append(loglike)

    score = -1*np.mean(score_list)

    #print('----')
    #print(score)
    #print('****')
    return score

def quadrant_genotype(typing, wt_min, mut_min):
    typing['quadrant_class'] = 'NA'
    for i in typing.index:
        if typing.loc[i,'transf_WT'] < wt_min:
            if typing.loc[i,'transf_MUT'] < mut_min:
                typing.loc[i,'quadrant_class'] = 'NA'
            else:
                typing.loc[i,'quadrant_class'] = 'MUT'
        else:
            if typing.loc[i,'transf_MUT'] < mut_min:
                typing.loc[i,'quadrant_class'] = 'WT'
            else:
                typing.loc[i,'quadrant_class'] = 'HET'
    
    return typing
    
def KDE_2D(typing, wt_min, mut_min, sample_dir=""):
    np.random.seed(0)
    noise = np.random.multivariate_normal(np.array([0,0]), 
                                          np.array([[0.09,0],[0,0.09]]), 
                                          typing.shape[0])
    data = typing[['transf_WT', 'transf_MUT']].values
    data_smooth = data + noise
    
    print("Computing optimal bandwidth...")
    time1 = timeit.default_timer()
    
    bw_list = []
    results_list = []
    for i in np.arange(0.1,1.0,0.025):
        bw_list.append(i)
        res = optimize_bandwidth([i], data_smooth)
        results_list.append(res)
    
    fit = polyfit(bw_list, results_list, 10)
    poly = Polynomial(fit)
    
    poly_x = np.linspace(0.1,1.0,1000)
    poly_y = poly(poly_x)
    plt.scatter(bw_list, results_list)
    plt.plot(poly_x, poly_y)
    plt.ylabel("Loss")
    plt.xlabel("Bandwidth")
    plt.savefig(sample_dir+"total_bandwidth.pdf", 
                dpi=500, bbox_inches = "tight")
    plt.show()
    plt.clf()
    
    bw_guess = bw_list[np.argmin(results_list)]
    result = minimize(poly, x0=np.array([bw_guess]), 
                      options={'disp':True}, bounds=((0.1, 0.975),))
    bw = result.x[0]
    
    time2 = timeit.default_timer()
    print("Optimal bandwidth computed in {} seconds.".format(
        time2-time1))
    
    x = data_smooth[:,0]
    y = data_smooth[:,1]
    xy = np.vstack([x,y])
    z = gaussian_kde(xy, bw_method=bw)(xy)#, bw_method=bw
    
    fig, ax = plt.subplots()
    ax.scatter(x, y, c=z, s=7, cmap='plasma')#, edgecolor='')
    plt.axhline(y=mut_min, color='r', linestyle='-')
    plt.axvline(x=wt_min, color='r', linestyle='-')
    #ax.colorbar()
    plt.title('Density')
    plt.xlabel('WT Z-scores')
    plt.ylabel('MUT Z-scores')
    plt.savefig(sample_dir+"total_smoothed_density.pdf", 
                dpi=500, bbox_inches = "tight")
    plt.show()
    plt.clf()

    return typing, data_smooth, bw

def cluster_number(typing, data_smooth, bw, wt_min, mut_min, sample_dir=""):
    X = euclidean_distances(data_smooth, squared=True)
    
    def gauss_func(x):
        coef1 = 1/(bw*np.sqrt(2*np.pi))
        coef2 = -1/(2*bw**2)
        return coef1*np.exp(coef2*x)
    
    affinity = gauss_func(X)
    
    for label in np.unique(typing['quadrant_class']):
        try:
            plt.scatter(data_smooth[np.where(typing['quadrant_class']==label)[0],0], 
                        data_smooth[np.where(typing['quadrant_class']==label)[0],1], label=label, 
                        s=7, color=gen_cmap[label])
        except:
            pass
    
    plt.legend()
    plt.title('Quadrant Genotype')
    plt.xlabel('WT Z-scores')
    plt.ylabel('MUT Z-scores')
    plt.axhline(y=mut_min, color='r', linestyle='-')
    plt.axvline(x=wt_min, color='r', linestyle='-')
    plt.savefig(sample_dir+"smoothed_quadrant.pdf", 
                dpi=500, bbox_inches = "tight")
    plt.show()
    plt.clf()

    rand_list = []
    clusters_list = []
    for k in [2,3,4]:
        spectral = SpectralClustering(n_clusters=k, affinity='precomputed')
        clusters_k = spectral.fit_predict(affinity)
        
        x=data_smooth[:,0]
        y=data_smooth[:,1]
        
        plt.scatter(x,y,c=clusters_k, s=7, cmap='Dark2')
        plt.title('Spectral Clustering (k={})'.format(k))
        plt.xlabel('WT Z-scores')
        plt.ylabel('MUT Z-scores')
        
        plt.axhline(y=mut_min, color='r', linestyle='-')
        plt.axvline(x=wt_min, color='r', linestyle='-')
        
        plt.savefig(sample_dir+"spectral_{}.pdf".format(k), 
                    dpi=500, bbox_inches = "tight")
        plt.show()
        plt.clf()
    
        rand_list.append(
            adjusted_rand_score(typing['quadrant_class'], 
                                clusters_k))
        clusters_list.append(clusters_k)
        
    rand_min = np.argmax(rand_list)
    n_clusters = rand_min + 2
    opt_clusters = clusters_list[rand_min]
    
    print("Number of optimal clusters: "+str(n_clusters))
        
    return typing, n_clusters, opt_clusters, affinity

def KDE_mixture_model(typing, bw, wt_min, mut_min, 
                      data_smooth, affinity, 
                      opt_clusters, sample_dir=""):
    
    data = typing[['transf_WT', 'transf_MUT']].values
    cluster_kde = {}
    cluster_weight_init = {}
    n = data_smooth.shape[0]
    mixture_prob = pd.DataFrame(index=typing.index,
                                columns=np.unique(opt_clusters))
    
    for c in np.unique(opt_clusters):
        cluster_index = np.where(opt_clusters==c)[0]
        kde = KernelDensity(bandwidth=bw)
        kde.fit(data_smooth[cluster_index,:])
        
        cluster_kde[c] = kde
        cluster_weight_init[c] = len(cluster_index)/n
        mixture_prob.loc[:,c] = np.exp(kde.score_samples(data))

    adjust_prob = mixture_prob.multiply(cluster_weight_init, axis=1)
    new_weights = adjust_prob.divide(adjust_prob.sum(axis=1), axis=0)
    new_weights = new_weights.sum(axis=0)
    new_weights = new_weights/typing.shape[0]
    
    adjust_prob = mixture_prob.multiply(new_weights, axis=1)
    adjust_prob = adjust_prob.divide(adjust_prob.sum(axis=1), axis=0)
    
    typing['clusters'] = adjust_prob.idxmax(axis=1)
    
    centroids = {}
    cluster_types = {}
    for i in np.unique(typing['clusters']):
        centroids[i] = \
            np.median(typing.loc[typing['clusters']==i, 
                                 ['transf_WT', 'transf_MUT']], axis=0)
        if centroids[i][0] < wt_min:
            if centroids[i][1] < mut_min:
                cluster_types[i] = 'NA'
            else:
                cluster_types[i] = 'MUT'
        else:
            if centroids[i][1] < mut_min:
                cluster_types[i]  = 'WT'
            else:
                cluster_types[i] = 'HET'
    
    adjust_prob = adjust_prob.rename(cluster_types, axis=1)
    adjust_prob = adjust_prob.groupby(adjust_prob.columns, axis=1).sum()
    
    x=data[:,0]#typing.loc[:,['transf_wt']]
    y=data[:,1]#typing.loc[:,['transf_mut']]
    
    for c in adjust_prob.columns:
        plt.scatter(x,y,c=adjust_prob.loc[:,c], s=7, cmap='plasma')
        plt.title('Cluster {} Probability'.format(c))
        plt.xlabel('WT Z-scores')
        plt.ylabel('MUT Z-scores')
        plt.colorbar()
        plt.savefig(sample_dir+"cluster_{}_probability.pdf".format(c), 
                    dpi=500, bbox_inches = "tight")
        plt.show()
        plt.clf()
    
    print("Remove uncertain labels and re-label with self-training.")
    
    typing['clusters'] = -1
    for i in typing.index:
        for j in adjust_prob.columns:
            if adjust_prob.loc[i,j] >= 0.99:
                typing.loc[i, 'clusters'] = j
                
    for label in np.unique(typing['clusters'].astype(str)):
        try:
            plt.scatter(data[np.where(typing['clusters'].astype(str)==label)[0],0], 
                        data[np.where(typing['clusters'].astype(str)==label)[0],1], label=label, 
                        s=7, color=gen_cmap[label])
        except:
            pass
    plt.legend()
    plt.title('Confident Calls')
    plt.xlabel('WT Z-scores')
    plt.ylabel('MUT Z-scores')
    plt.savefig(sample_dir+"confident_calls.pdf", 
                dpi=500, bbox_inches = "tight")
    plt.show()
    plt.clf()    

    n = data.shape[0]
    n_neighbors = round(np.sqrt(n))
    
    def dist_weight(distances):
        gamma = 1/(data.shape[1]*data.var()) #rule of thumb
        return np.exp(-1*gamma*distances)
    
    time1 = timeit.default_timer()
    STC = SelfTrainingClassifier(KNeighborsClassifier(n_neighbors=n_neighbors, weights=dist_weight),
                            criterion='k_best', k_best=1, max_iter=None)
    STC.fit(data, typing['clusters'])
    time2 = timeit.default_timer()
    print("Label propagation completed in {} seconds.".format(
        time2-time1))
    
    cluster_pred = STC.transduction_
    typing['clusters'] = cluster_pred
    
    typing['genotype_pred'] = typing['clusters']
    del typing['clusters']
    
    ## plot final figures
    x=data[:,0]#typing.loc[:,['transf_wt']]
    y=data[:,1]#typing.loc[:,['transf_mut']]
    
    for label in np.unique(typing['genotype_pred'].astype(str)):
        try:
            plt.scatter(data[np.where(typing['genotype_pred'].astype(str)==label)[0],0], 
                        data[np.where(typing['genotype_pred'].astype(str)==label)[0],1], label=label,
                       s=7, color=gen_cmap[label])
        except:
            pass
    plt.legend()
    plt.title('Cluster Genotype')
    plt.xlabel('WT Z-scores')
    plt.ylabel('MUT Z-scores')
    
    plt.axhline(y=mut_min, color='r', linestyle='-')
    plt.axvline(x=wt_min, color='r', linestyle='-')
    
    #plt.colorbar()
    plt.savefig(sample_dir+"cluster_genotype.pdf", 
                dpi=500, bbox_inches = "tight")
    plt.show()
    plt.clf() 
    
    for label in np.unique(typing['quadrant_class']):
        try:
            plt.scatter(data[np.where(typing['quadrant_class']==label)[0],0], 
                        data[np.where(typing['quadrant_class']==label)[0],1], label=label,
                       s=7, color=gen_cmap[label])
        except:
            pass
    plt.title('Quadrant Genotype')
    plt.xlabel('WT Z-scores')
    plt.ylabel('MUT Z-scores')
    
    plt.axhline(y=mut_min, color='r', linestyle='-')
    plt.axvline(x=wt_min, color='r', linestyle='-')
    
    plt.legend()
    plt.savefig(sample_dir+"quadrant_genotype.pdf", 
                dpi=500, bbox_inches = "tight")
    plt.show()
    plt.clf() 
    
    print('Percent genotyped: (cluster method): ' + 
          str(sum(typing['genotype_pred']!='NA')/typing.shape[0]))
    
    print('Percent genotyped: (quadrant method): ' + 
          str(sum(typing['quadrant_class']!='NA')/typing.shape[0]))
    
    print("Quadrant labels:")
    print(Counter(typing['quadrant_class']))
    
    print("Cluster labels:")
    print(Counter(typing['genotype_pred']))
    
    return typing





    
    
    
    
    
