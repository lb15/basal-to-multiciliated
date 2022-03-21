#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 12:26:54 2019

@author: lb
"""

import scanpy as sc
import pandas as pd
import numpy as np

path='/Users/lb/Documents/Reiter_Seq/10X_042519/NS_analysis/'
results_file = '/Users/lb/Documents/Reiter_Seq/10X_042519/v16/NS_v16_scanpy.h5ad'  # the file that will store the analysis results

sc.settings.set_figure_params(dpi=80)

adata = sc.read_10x_mtx(path + 'filtered_feature_bc_matrix',  # the directory with the `.mtx` file
    var_names='gene_symbols',                  # use gene symbols for the variable names (variables-axis index)
    cache=True)                                # write a cache file for faster subsequent reading

adata

sc.pl.highest_expr_genes(adata, n_top=20)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

mito_genes = adata.var_names.str.startswith('mt-')
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1


sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],jitter=0.4, multi_panel=True)

adata = adata[adata.obs['n_genes'] > 3500, :]
adata = adata[adata.obs['percent_mito'] < 0.10, :]
adata = adata[adata.obs['n_counts'] < 120791, :]
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],jitter=0.4, multi_panel=True)


sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)

adata.raw = adata

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)

adata = adata[:, adata.var['highly_variable']]


# error with scipy.misc 
#sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)

adata.write(results_file)

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)

sc.tl.umap(adata)

sc.pl.umap(adata, color=['Foxj1', 'Krt5', 'Scgb3a1'])

sc.tl.louvain(adata)
sc.pl.umap(adata, color=['louvain'])

adata.write(results_file)


####### doublets ##########

import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os

counts_matrix = adata.raw.X
print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=.035)
doublet_scores, predicted_doublets = scrub.scrub_doublets()

predicted_doublets = scrub.call_doublets(threshold=0.16)

scrub.plot_histogram();

print('Running UMAP...')
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True);
#plt.savefig('/Users/lb/Documents/Reiter_Seq/10X_070219/v13/WtAd3_scrublet_12perc_thresh28_umap.png')

sc.pl.umap(adata, color=['Mycl', 'Krt5', 'Scgb3a1'])

cellbarc = adata.obs.index

df=pd.DataFrame({
        'cell_barcodes':cellbarc,
        'doublet_score':scrub.doublet_scores_obs_,
        'predicted_doublet':scrub.predicted_doublets_
        })
df.to_csv(path+'v16/NS_scrublet_3.5perc_thresh0.16_output_table.csv',index=False)
