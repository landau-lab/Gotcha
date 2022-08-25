# -*- coding: utf-8 -*-
"""
@author: Sanjay Kottapalli (svk4001@med.cornell.edu)
@date: 08/11/2022
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
#from scipy.signal import savgol_filter
#from scipy.signal import find_peaks
from scipy.integrate import quad
from collections import Counter
#from sklearn.metrics.pairwise import euclidean_distances
from scipy.stats import gaussian_kde
#from sklearn.metrics import homogeneity_score
#from sklearn.cluster import SpectralClustering
from sklearn.semi_supervised import SelfTrainingClassifier
from sklearn.neighbors import KNeighborsClassifier

#from numpy.random import multivariate_normal
from sklearn.metrics import confusion_matrix
from scipy.stats import gmean

from joblib import Parallel, delayed

import numpy.matlib

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
            typing, wt_min = noise_correct(typing,i,sample_dir,sample_id)
        else:
            typing, mut_min = noise_correct(typing,i,sample_dir,sample_id)
    
    print("Performing quadrant genotyping.")
    typing = quadrant_genotype(typing, wt_min, mut_min)
        
    print("Computing KNN-based clusters.")
    typing = KNN_cluster(typing, wt_min, mut_min, 
                      1.0, sample_dir)
    
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
    cell_line = pd.read_csv(infile, index_col=0, sep=",")
    #print(cell_line.head())
    if sample_id in np.unique(cell_line['Sample'].values):
        cell_line = cell_line.loc[cell_line['Sample']==sample_id, :]

    genotyping = pd.DataFrame(index=cell_line.index)
    genotyping['WTcount'] = cell_line[gene_id+'_WTcount']
    genotyping['MUTcount'] = cell_line[gene_id+'_MUTcount']
    index = genotyping[['WTcount', 'MUTcount']].dropna().index
    genotyping = genotyping.loc[index,:]
    print("Number of cells: "+str(genotyping.shape[0]))
    
    return genotyping

def noise_correct(typing, feature="", sample_dir="", sample_id=""):
    np.random.seed(0)
    pseudocount = 1
    X = typing[[feature+'count']]+pseudocount
    logged_counts = np.log(X) #change to log10 later, arcsinh? not log2

    #noise = np.random.normal(0,0.0,typing.shape[0]) #no noise
    #logged_counts = logged_counts + noise.reshape(-1,1)
    typing['transf_{}'.format(feature)] = logged_counts

    plt.hist(logged_counts, density=True, bins=50)
    plt.title(feature+' counts')
    plt.ylabel("Probability")
    plt.xlabel("Log(counts+1)")
    #plt.savefig(sample_dir+"log_{}_counts.pdf".format(feature), 
    #            dpi=500, bbox_inches = "tight")
    plt.show()
    plt.clf()
        
    bw = 0.1
    
    kde = KernelDensity(bandwidth=bw)#kernel='exponential')
    kde.fit(typing['transf_{}'.format(feature)].values.reshape(-1, 1))
    x_bin = np.histogram(typing['transf_{}'.format(feature)], bins=50)[1]
    kde_x = np.linspace(min(x_bin)-0.5,max(x_bin)+0.5,10000)
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

    kde_test = kde_smooth
    
    def kde_func(x):
        return np.exp(kde.score_samples(x.reshape(1, -1))[0])
    
    plt.plot(kde_x, np.log(kde_test), color='red')
    plt.title("{} (all drops)".format(sample_id))
    plt.ylabel("Log-Probability")
    plt.xlabel("Log(WTcounts)")
    plt.show()
    plt.clf()
    
    kde_smooth2 = kde_test#np.log(kde_test)
    #kde_smooth2 = savgol_filter(kde_smooth2, 1001, 3)#, mode='constant', cval=-np.inf)
    plt.plot(kde_x, kde_smooth2, color='red')
    plt.show()
    plt.clf()
    
    fit = polyfit(kde_x, kde_smooth2, 12)
    poly = Polynomial(fit)
    poly_y = poly(kde_x)
    plt.plot(kde_x, poly_y, color='red')
    plt.show()
    plt.clf()
    
    mid = (max(kde_x)-min(kde_x))/2
    print(mid)
    result = minimize(poly, x0=np.array([mid]), 
                        options={'disp':False}, bounds=((min(kde_x), max(kde_x)),),
                     method='Nelder-Mead')
    new_min = result.x[0]
    print(new_min)
    
    result = minimize(kde_func, x0=np.array([new_min]), 
                        options={'disp':True}, bounds=((min(kde_x), max(kde_x)),),
                     method='Nelder-Mead')
    new_min = result.x[0]
    print(new_min)
    
    if feature=="WT":
        alt_feature="MUT"
    else:
        alt_feature="WT"

    if new_min==max(kde_x):
        zero_counts = typing['{}count'.format(alt_feature)]
        zero_counts = zero_counts[zero_counts==0].index

        new_min = np.percentile(typing.loc[zero_counts,'transf_{}'.format(feature)], 99.99)

    noise_values = logged_counts
    noise_values = noise_values[noise_values<new_min]
    noise_values = noise_values.dropna()
    prop_noise = len(noise_values)/len(logged_counts)

    kde_noise = KernelDensity(bandwidth=bw)#kernel='exponential')
    kde_noise.fit(noise_values.values.reshape(-1, 1))
    x_bin = np.histogram(typing['transf_{}'.format(feature)], bins=50)[1]
    #kde_x = np.linspace(min(x_bin),max(x_bin),1000)
    noise_smooth = np.exp(kde_noise.score_samples(kde_x.reshape(-1, 1)))
    
    signal_values = logged_counts
    signal_values = signal_values[signal_values>=new_min]
    signal_values = signal_values.dropna()
    prop_signal = 1-prop_noise

    kde_signal = KernelDensity(bandwidth=bw)#kernel='exponential')
    kde_signal.fit(signal_values.values.reshape(-1, 1))
    #x_bin = np.histogram(typing['transf_{}'.format(feature)], bins=50)[1]
    #kde_x = np.linspace(min(x_bin),max(x_bin),1000)
    signal_smooth = np.exp(kde_signal.score_samples(kde_x.reshape(-1, 1)))
    
    plt.hist(typing['transf_{}'.format(feature)], density=True, bins=50)
    plt.plot(kde_x, kde_test, color='red')
    plt.plot(kde_x, prop_noise*noise_smooth, color='pink')
    plt.plot(kde_x, prop_signal*signal_smooth, color='yellow')
    plt.title('Noise vs. Signal')
    plt.ylabel("Probability")
    plt.xlabel("Log(counts+1)")
    #plt.xticks(np.arange(0,6,1.0))
    plt.savefig(sample_dir+"{}_kde_mixture.pdf".format(feature), 
                dpi=500, bbox_inches = "tight")
    plt.show()
    plt.clf()
    
    def noise_func(x):
        return np.exp(kde_noise.score_samples(np.array(x).reshape(1, -1)))

    def mean_func(x):
        return x*noise_func(x)

    def var_func(x):
        return x*x*noise_func(x)
    
    noise_mean = quad(func=mean_func, a=-np.inf, b=np.inf, limit=100)[0]
    noise_var = quad(func=var_func, a=-np.inf, 
                     b=np.inf, limit=100)[0] - noise_mean**2
    print("Noise mean, noise var:")
    print(noise_mean, noise_var)
    
    noise_std = np.sqrt(noise_var)
    #typing['transf_{}'.format(feature)] = (logged_counts-noise_mean)/noise_std

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
    
    feat_min = (new_min-noise_mean)/new_std
    print(feat_min)
    
    return typing, feat_min

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

def KNN_cluster(typing, wt_min, mut_min, 
                      knn_window=1.0, sample_dir=""):
    
    typing['clusters'] = typing['quadrant_class']
    
    data = typing[['transf_WT', 'transf_MUT']].values
    indices1 = set(np.where(data[:,0]>wt_min-knn_window)[0])
    indices1 = indices1.intersection(set(np.where(data[:,0]<=wt_min+knn_window)[0]))
    indices2 = set(np.where(data[:,1]>mut_min-knn_window)[0])
    indices2 = indices2.intersection(set(np.where(data[:,1]<=mut_min+knn_window)[0]))
    indices = indices1.union(indices2)
    
    indices = np.array(list(indices))
    indices = typing.index[indices]
    
    typing.loc[indices,'clusters'] = -1
    
    n = list(dict(Counter(typing['quadrant_class'].values)).values())
    print(n)
    n_neighbors = round(np.sqrt(gmean(n)))
    print("Nearest neighbors: "+str(n_neighbors))
    
    print("Remove uncertain labels and re-label with self-training.")
    
    for label in np.unique(typing['clusters'].astype(str)):
        try:
            plt.scatter(data[np.where(typing['clusters'].astype(str)==label)[0],0], 
                        data[np.where(typing['clusters'].astype(str)==label)[0],1], label=label, 
                        s=7, color=gen_cmap[label])
        except:
            pass
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    #plt.legend()
    plt.title('Confident Calls')
    plt.xlabel('WT Z-scores')
    plt.ylabel('MUT Z-scores')
    plt.savefig(sample_dir+"confident_calls.pdf", 
                dpi=500, bbox_inches = "tight")
    plt.show()
    plt.clf()
    
    time1 = timeit.default_timer()
    STC = SelfTrainingClassifier(KNeighborsClassifier(n_neighbors=n_neighbors, weights='distance'),#gauss_func),
                            criterion='k_best', k_best=1, max_iter=None)
    STC.fit(data, typing['clusters'])
    time2 = timeit.default_timer()
    print("Label propagation completed in {} seconds.".format(
        time2-time1))
    
    cluster_pred = STC.transduction_
    typing['clusters'] = cluster_pred
    
    typing['genotype_pred'] = typing['clusters']
    del typing['clusters']
    
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
    
    print('Proportion genotyped: (cluster method): ' + 
          str(sum(typing['genotype_pred']!='NA')/typing.shape[0]))
    
    print('Proportion genotyped: (quadrant method): ' + 
          str(sum(typing['quadrant_class']!='NA')/typing.shape[0]))
    
    print("Quadrant labels:")
    print(Counter(typing['quadrant_class']))
    
    print("Cluster labels:")
    print(Counter(typing['genotype_pred']))
    
    return typing





    
    
    
    
    
