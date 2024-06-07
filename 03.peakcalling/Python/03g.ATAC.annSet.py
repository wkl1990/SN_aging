import os
import sys

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

import snapatac2 as sa2
import anndata

import pyreadr

proj_dir = "/SN_aging"

meta_file = os.path.join(proj_dir, "rds/RNA/cluster/integration/meta/SN.allen.integration.meta.final.barcode.csv")
metadata = pd.read_csv(meta_file)
metadata['barcode2'] = [f"{i}.{j}" for i, j in zip(metadata['orig.ident'].to_list(), metadata['barcode'].to_list())]
metadata.set_index("barcode2", inplace = True)
samples = metadata["orig.ident"].unique().tolist()

pmat_dir = os.path.join(proj_dir, "rds/ATAC/after_integra/pmat/")

# * get AnnDataSet and AnnData for union peak
suffix = ".pmat.h5ad"
sample2files = [(sample, os.path.join(pmat_dir, f"{sample}{suffix}")) for sample in samples]
out_pmat = os.path.join(proj_dir, "rds/ATAC/after_integra/pmat/", "ATAC.combined.pmat.h5ad")
#if not os.path.exists(os.path.dirname(out_pmat)):
#    os.makedirs(os.path.dirname(os.path.dirname(out_pmat)), exist_ok = True)
    
conset = sa2.AnnDataSet(adatas = sample2files, filename = out_pmat, add_key = 'sample')
conset.close()

# load data
conset = sa2.read_dataset(filename = out_pmat, mode = 'r')
out_adata = os.path.join(proj_dir, "rds/ATAC/after_integra/", "combine/ATAC.combined.pmat.adata.h5ad")
adata = conset.to_adata(copy_x = True, file = out_adata)
conset.close()

def modify_obs_name(sds, obs_key = "sample"):
    obs_names = [f"{i}.{j}"
        for i, j in zip(
            sds.obs[obs_key].to_list(), sds.obs_names)]
    return obs_names

new_obs_names = modify_obs_name(adata, obs_key = "sample")
adata.obs_names = new_obs_names
adata.close()

# for shuffle peaks
# * get AnnDataSet and AnnData for union peak
shuffle_suffix = ".shuffle.pmat.h5ad"
shuffle_sample2files = [(sample, os.path.join(pmat_dir, f"{sample}{shuffle_suffix}")) for sample in samples]
shuffle_out_pmat = os.path.join(proj_dir, "rds/ATAC/after_integra/pmat/", "ATAC.shuffle.combined.pmat.h5ad")
#if not os.path.exists(os.path.dirname(shuffle_out_pmat)):
#    os.makedirs(os.path.dirname(os.path.dirname(shuffle_out_pmat)), exist_ok = True)
    
shuffle_conset = sa2.AnnDataSet(adatas = shuffle_sample2files, filename = shuffle_out_pmat, add_key = 'sample')
shuffle_conset.close()

# load data
shuffle_conset = sa2.read_dataset(filename = shuffle_out_pmat, mode = 'r')
shuffle_out_adata = os.path.join(proj_dir, "rds/ATAC/after_integra/", "combine/ATAC.shuffle.combined.pmat.adata.h5ad")
shuffle_adata = shuffle_conset.to_adata(copy_x = True, file = shuffle_out_adata)
shuffle_conset.close()

new_obs_names = modify_obs_name(shuffle_adata, obs_key = "sample")
shuffle_adata.obs_names = new_obs_names
shuffle_adata.close()


