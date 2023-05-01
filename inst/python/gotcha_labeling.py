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
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
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
gen_cmap = {'MUT':'#CB2026','WT':'#38459C',
            'HET':'#F89622','NA':'#D3D3D3','-1':tab20.colors[9]}

plt.clf()

def GotchaLabeling(path="", infile="", gene_id="", sample_id="", sample_column="", saturation=False):
    time1 = timeit.default_timer()
    print("Reading in file.")
    typing = read_data(path+infile, gene_id, sample_id, sample_column)
    sample_dir = path+sample_id+'/'
    try:
        os.mkdir(sample_dir)
    except:
        pass
    
    for i in ["WT", "MUT"]:
        print("Noise correcting {} read counts.".format(i))
        if i=="WT":
            typing, wt_min, wt_thresh, wt_auto = noise_correct(typing,i,sample_dir,sample_id, 'auto', 0.15, 0.25)
        else:
            typing, mut_min, mut_thresh, mut_auto = noise_correct(typing,i,sample_dir,sample_id, 'auto', 0.15, 0.25)
    
    pairs = [(a,b) for a in wt_thresh for b in mut_thresh]
    pairs = list(set(pairs))
    pairs.remove((wt_auto, mut_auto))
    print(pairs)
    
    print("Performing quadrant genotyping.")
    typing = quadrant_genotype(typing, wt_min, mut_min)
        
    print("Computing KNN-based clusters.")
    typing = KNN_cluster(typing, wt_min, mut_min, 
                      5, sample_dir)
    
    typing.to_csv(sample_dir+sample_id+'_genotype_labels.csv')
    
    typing_orig = typing.copy()
    
    if saturation:
        print("Performing saturation analysis.")
        saturation_analysis(typing, sample_dir)
        
    ########
    for i in range(len(pairs)):
        new_dir = sample_dir+'Output_{}/'.format(i+1)
        try:
            os.mkdir(new_dir)
        except:
            pass
        
        print("Analyzing case {}.".format(i+1))
        
        for j in ["WT", "MUT"]:
            print("Noise correcting {} read counts.".format(i))
            if j=="WT":
                typing, wt_min, wt_thresh, wt_auto = noise_correct(typing,j,new_dir,sample_id, pairs[i], 0.15, 0.25)
            else:
                typing, mut_min, mut_thresh, mut_auto = noise_correct(typing,j,new_dir,sample_id, pairs[i], 0.15, 0.25)
        print("Performing quadrant genotyping.")
        typing = quadrant_genotype(typing, wt_min, mut_min)
        
        print("Computing KNN-based clusters.")
        typing = KNN_cluster(typing, wt_min, mut_min, 5, new_dir)
    
        typing.to_csv(new_dir+sample_id+'_genotype_labels.csv')
    
    ########
    
    print("All analysis complete!")
    time2 = timeit.default_timer()
    print("Total time to execute: {}".format(time2-time1))
    
    return typing_orig

def read_data(infile="", gene_id="", sample_id="", sample_column=""):
    '''
    This function reads in the sequencing reads from only real cells
    (based on ATAC calling). Cells without GoTChA reads are dropped.
    '''
    cell_line = pd.read_csv(infile, index_col=0, sep=",")
    
    cell_line.fillna(0.0, inplace=True)
    #print(cell_line.head())
    
    try:
        genotyping = pd.DataFrame(index=cell_line.index)
        cell_line[sample_column] = cell_line[sample_column].astype(str)
        if sample_id in np.unique(cell_line[sample_column].values):
            cell_line = cell_line.loc[cell_line[sample_column]==sample_id, :]
        genotyping['WTcount'] = cell_line[gene_id+'_WTcount']
        genotyping['MUTcount'] = cell_line[gene_id+'_MUTcount']
    except:
        cell_line = cell_line.set_index("WhiteListMatch")
        genotyping = pd.DataFrame(index=cell_line.index)
        genotyping['WTcount'] = cell_line['WTcount']
        genotyping['MUTcount'] = cell_line['MUTcount']
    
    #print(genotyping.head())
    index = genotyping[['WTcount', 'MUTcount']].dropna().index
    genotyping = genotyping.loc[index,:]
    print("Number of cells: "+str(genotyping.shape[0]))
    
    return genotyping

'''
def optimize_bandwidth(x, data_input):
    print(x)
    kde = KernelDensity(bandwidth=x[0])
    
    n = data_input.shape[0]
    n_splits=10
    score_list = []
    
    for i in range(10):
        train_index = np.random.choice(
            np.arange(n), size=round(n/n_splits), replace=False)
        data_train = data_input[train_index,:]
        #print(data_temp.shape)
        #print(data_temp)
        kde.fit(data_train)
        #print(data_train.shape)
        
        #data_test = data_input.copy()
        mask = np.ones(n, dtype=bool)
        mask[train_index] = False
        data_test = data_input[mask,:].copy()
        #print(data_test.shape)
        
        loglike = kde.score(data_test)
        score_list.append(loglike)
        #print(loglike)

    score = -1*np.mean(score_list)

    print('----')
    print(score)
    print('****')
    return score
'''

def noise_correct(typing, feature="", sample_dir="", sample_id="", quadrants='auto', bw=0.3, rel_height=0.25):
    np.random.seed(0)
    pseudocount = 1
    X = typing[[feature+'count']]+pseudocount
    logged_counts = np.log(X) #change to log10 later, arcsinh? not log2

    factor = 1.15
    logged_counts = logged_counts ** factor
    #noise = np.random.normal(0,0.1,typing.shape[0]) #no noise
    #logged_counts = logged_counts + noise.reshape(-1,1)
    
    typing['transf_{}'.format(feature)] = logged_counts
    
    
    ######################
    if not os.path.isfile(sample_dir+"initial_scatter.pdf"):
        fig, ax = plt.subplots()

        ax.scatter(np.log(typing[['WTcount']]+pseudocount)**factor, np.log(typing[['MUTcount']]+pseudocount)**factor, s=5)
        plt.title("Initial scatter plot")
        plt.xlabel("Log WT counts")
        plt.ylabel("Log MUT counts")

        start, end = ax.get_xlim()
        ax.xaxis.set_ticks(np.arange(0, end, 1.0))
        start, end = ax.get_ylim()
        ax.yaxis.set_ticks(np.arange(0, end, 1.0))

        plt.savefig(sample_dir+"initial_scatter.pdf", 
                    dpi=500, bbox_inches = "tight")
        plt.show()
        plt.clf()
    ################

    plt.hist(logged_counts, density=True, bins=50)
    plt.title(feature+' counts')
    plt.ylabel("Probability")
    plt.xlabel("Log(counts+1)")
    #plt.savefig(sample_dir+"log_{}_counts.pdf".format(feature), 
    #            dpi=500, bbox_inches = "tight")
    plt.show()
    plt.clf()
    
    #bw=0.3
          
    '''
    #bandwidth optimization
    print("Computing optimal bandwidth...")
    time1 = timeit.default_timer()
    data = typing[['transf_{}'.format(feature)]].values
    
    
    def optimize_bandwidth(x):
        print(x)
        kde = KernelDensity(bandwidth=x)

        data_input = logged_counts.values.copy()

        n = data_input.shape[0]
        n_splits=10
        score_list = []

        for i in range(10):
            train_index = np.random.choice(
                np.arange(n), size=round(n/n_splits), replace=False)
            data_train = data_input[train_index,:]
            #print(data_temp.shape)
            #print(data_temp)
            kde.fit(data_train)
            #print(data_train.shape)

            #data_test = data_input.copy()
            mask = np.ones(n, dtype=bool)
            mask[train_index] = False
            data_test = data_input[mask,:].copy()
            #print(data_test.shape)

            loglike = kde.score(data_test)
            score_list.append(loglike)
            #print(loglike)

        score = -1*np.mean(score_list)

        #print('----')
        #print(score)
        #print('****')
        return x, score
    
    
    bw_list = np.arange(0.05,0.3,0.01)
    r = Parallel(n_jobs=-1, verbose=10, max_nbytes=None)(delayed(optimize_bandwidth)(bw) 
                                        for bw in bw_list)
    bw_list, results_list = zip(*r)
    
    
    
    fit = polyfit(bw_list, results_list, 10)
    poly = Polynomial(fit)
    
    poly_x = np.linspace(0.05,0.29,1000)
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
                      options={'disp':True}, bounds=((0.05,0.3-0.01),))
    bw = result.x[0]
    
    time2 = timeit.default_timer()
    print("Optimal bandwidth computed in {} seconds.".format(
        time2-time1))
    print("Bandwidth: {}".format(bw))
    '''
    
    bw = bw
    print("Bandwidth try: {}".format(bw))
    
    
    kde = KernelDensity(bandwidth=bw)#kernel='exponential')
    kde.fit(typing['transf_{}'.format(feature)].values.reshape(-1, 1))
    x_bin = np.histogram(typing['transf_{}'.format(feature)], bins=50)[1]
    kde_x = np.linspace(min(x_bin)-0.50,max(x_bin)+0.50,10000)
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
    
    kde_smooth2 = kde_test#np.log(kde_test)
    
    '''
    if quadrants=='auto':
        print("Peaks:")
        peaks, properties = find_peaks(kde_test, width=(None,None), height=0.005, rel_height=rel_height)
        print(kde_x[peaks])
        print(properties)

        prop_array = properties['widths']/(properties['peak_heights']+1)
        print("Prop array: ")
        print(prop_array)

        modes = np.argsort(prop_array)[-2:]
        modes = peaks[modes]

        modes = np.sort(modes)
        print(modes)
        print(kde_x[modes])

        new_range = kde_test[modes[0]:modes[1]]
        #print(new_range)
        new_min = np.argmin(new_range) + modes[0]
        new_min = kde_x[new_min]
        print(new_min)
    if quadrants != 'auto':
        if feature=='WT':
            new_min = quadrants[0]
        else:
            new_min = quadrants[1]
    '''
    ##############
    
    print("Peaks:")
    peaks, properties = find_peaks(kde_test, width=(None,None), height=0.005, rel_height=rel_height)
    print(kde_x[peaks])
    print(properties)

    prop_array = properties['widths']#/(properties['peak_heights']+1)
    print("Prop array: ")
    print(prop_array)
    
    thresh_list = []

    if quadrants=='auto':
        modes = np.argsort(prop_array)[-3:]
        modes = peaks[modes]
        modes = np.sort(modes)
        print(modes)
        print(kde_x[modes])

        pairs = [np.sort([a, b]) for idx, a in enumerate(modes) for b in modes[idx + 1:]]
        thresh_list = []
        
        '''
        print("--------------------------------------")
        fit = polyfit(kde_x, kde_test, 20)
        poly = Polynomial(fit)
        poly_y = poly(kde_x)
        #poly_y = savgol_filter(kde_test, 3333, 3, mode='constant', cval=0.0)
        plt.scatter(kde_x, poly_y, s=7)
        plt.title("smoothed fit")
        plt.ylim([0.0,0.9]
        plt.show()
        plt.gcf()
        
        
        print('--------------------------------------')
        print("second derivative: ")
        dy=np.diff(poly_y,1)
        dx=np.diff(kde_x,1)
        yfirst=dy/dx
        xfirst=0.5*(kde_x[:-1]+kde_x[1:])
        
        dyfirst=np.diff(yfirst,1)
        dxfirst=np.diff(xfirst,1)
        ysecond=dyfirst/dxfirst

        xsecond=0.5*(xfirst[:-1]+xfirst[1:])
        
        plt.scatter(xsecond, ysecond, s=7)
        plt.title("second derivative")
        plt.show()
        plt.gcf()
        '''
        
        for i in range(len(pairs)):
            new_range = kde_test[pairs[i][0]:pairs[i][1]]
            new_min = np.argmin(new_range) + pairs[i][0]
            new_min = kde_x[new_min]
            thresh_list.append(new_min)
            
        modes = np.argsort(prop_array)[-2:]
        modes = peaks[modes]
        modes = np.sort(modes)
        
        new_range = kde_test[modes[0]:modes[1]]
        #print(new_range)
        new_min = np.argmin(new_range) + modes[0]
        new_min = kde_x[new_min]
        print(new_min)
        
    if quadrants != 'auto':
        if feature=='WT':
            new_min = quadrants[0]
        else:
            new_min = quadrants[1]
    
    ##############
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

    added_var = 0.1**2
    new_var = noise_var - added_var
    new_std = np.sqrt(new_var)

    typing['transf_{}'.format(feature)] = (np.log(X)**factor-noise_mean)/new_std

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
    
    return typing, feat_min, thresh_list, new_min

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
                      knn_window=5, sample_dir=""):
    
    typing['clusters'] = typing['quadrant_class']
    
    data = typing[['transf_WT', 'transf_MUT']].values
    
    for label in np.unique(typing['quadrant_class']):
        try:
            plt.scatter(data[np.where(typing['quadrant_class']==label)[0],0], 
                        data[np.where(typing['quadrant_class']==label)[0],1], label=label,
                       s=5, color=gen_cmap[label])
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
    
    '''
    range_wt = max(data[:,0])-min(data[:,0])
    range_mut = max(data[:,0])-min(data[:,0])
    indices1 = set(np.where(data[:,0]>wt_min-knn_window*range_wt)[0])
    indices1 = indices1.intersection(set(np.where(data[:,0]<=wt_min+knn_window*range_wt)[0]))
    indices2 = set(np.where(data[:,1]>mut_min-knn_window*range_mut)[0])
    indices2 = indices2.intersection(set(np.where(data[:,1]<=mut_min+knn_window*range_mut)[0]))
    indices = indices1.union(indices2)
    
    indices = np.array(list(indices))
    indices = typing.index[indices]
    '''
    
    #knn_window = 5
    
    print("2D Density plot")
    x = data[:,0]
    y = data[:,1]
    
    xy = np.vstack([x,y])
    z = np.log(gaussian_kde(xy, bw_method=0.1)(xy))#, bw_method=bw
    
    fig, ax = plt.subplots()
    ax.scatter(x, y, c=z, s=5, cmap='plasma')#, edgecolor='')
    plt.axhline(y=mut_min, color='r', linestyle='-')
    plt.axvline(x=wt_min, color='r', linestyle='-')
    #ax.colorbar()
    plt.title('Density')
    plt.xlabel('WT Z-scores')
    plt.ylabel('MUT Z-scores')
    plt.savefig(sample_dir+"total_density.pdf", 
                dpi=500, bbox_inches = "tight")
    plt.show()
    
    #print(knn_window)
    
    wt_percentile = (len(x[x<wt_min])/len(x))*100
    print("WT quadrant percentile: "+str(wt_percentile))
    #print(np.percentile(x, wt_percentile-knn_window))
    #print(np.percentile(x, wt_percentile+knn_window))
    
    mut_percentile = (len(y[y<mut_min])/len(y))*100
    print("MUT quadrant percentile: "+str(mut_percentile))
    #print(np.percentile(y, mut_percentile-knn_window))
    #print(np.percentile(y, mut_percentile+knn_window))
    
    if wt_percentile<99-knn_window and wt_percentile>=knn_window+1:
        indices1 = set(np.where(x>np.percentile(x, wt_percentile-knn_window))[0])
        indices1 = indices1.intersection(set(np.where(x<np.percentile(x, wt_percentile+knn_window))[0]))
    elif wt_percentile<knn_window+1:
        indices1 = set(np.where(x>np.percentile(x, wt_percentile))[0])
        indices1 = indices1.intersection(set(np.where(x<np.percentile(x, wt_percentile+knn_window))[0]))
    else:
        indices1 = set(np.where(x>np.percentile(x, wt_percentile-knn_window))[0])
        indices1 = indices1.intersection(set(np.where(x<np.percentile(x, wt_percentile))[0]))
        
    if mut_percentile<99-knn_window and mut_percentile>=knn_window+1:
        indices2 = set(np.where(y>np.percentile(y, mut_percentile-knn_window))[0])
        indices2 = indices2.intersection(set(np.where(y<np.percentile(y, mut_percentile+knn_window))[0]))
    elif mut_percentile<knn_window+1:
        indices2 = set(np.where(y>np.percentile(y, mut_percentile))[0])
        indices2 = indices2.intersection(set(np.where(y<np.percentile(y, mut_percentile+knn_window))[0]))
    else:
        indices2 = set(np.where(y>np.percentile(y, mut_percentile-knn_window))[0])
        indices2 = indices2.intersection(set(np.where(y<np.percentile(y, mut_percentile))[0]))
    
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
                        s=5, color=gen_cmap[label])
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
                       s=5, color=gen_cmap[label])
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
    
    print('Proportion genotyped: (cluster method): ' + 
          str(sum(typing['genotype_pred']!='NA')/typing.shape[0]))
    
    print('Proportion genotyped: (quadrant method): ' + 
          str(sum(typing['quadrant_class']!='NA')/typing.shape[0]))
    
    print("Quadrant labels:")
    print(Counter(typing['quadrant_class']))
    
    print("Cluster labels:")
    print(Counter(typing['genotype_pred']))
    
    return typing
    
def saturation_analysis(typing, sample_dir):
    n_cells = typing.shape[0]
    typing['WhiteListMatch']=typing.index
    typing['total_reads']=typing['MUTcount']+typing['WTcount']
    typing['total_reads']=typing['total_reads'].astype(int)
    typing['index_start']=1

    typing['index'] = typing.apply(lambda row: list(range(row['index_start'], row['total_reads']+1)), axis=1)
    typing = typing.explode('index')
    typing['read'] = typing['WhiteListMatch']+typing['index'].astype(str)

    saturations = {np.round(i,2): len(typing.sample(int(len(typing)*i))['WhiteListMatch'].unique()) for i in np.arange(0.0,1.1,0.1)}

    saturation_df = pd.DataFrame.from_dict(data=[saturations]).transpose()
    
    saturation_df = saturation_df/n_cells
    
    plt.plot(saturation_df.index,saturation_df[0])
    plt.savefig(sample_dir+"saturation.pdf", 
                dpi=500, bbox_inches = "tight")
    
    return None




    
    
    
    
    
