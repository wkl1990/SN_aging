import os
import sys

#import logging
from pathlib import Path
from typing import Dict, List, Tuple
from multiprocessing import Pool
import pickle

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

import snapatac2 as sa2
import anndata
# import scanpy as sc

#import pyprojroot
import pyreadr


def modify_obs_name(sds: sa2.AnnData | sa2.AnnDataSet,
                    obs_key = "sample") -> List[str]:
    obs_names: List[str] = [f"{i}.{j}"
              for i, j in zip(
                      sds.obs[obs_key].to_list(), sds.obs_names)]
    return obs_names

def rename_allensubclass(sc: str) -> str:
    return sc.replace(" ", '_').replace("/", "-")

def read_bed4(fnm) -> pd.DataFrame:
    bed: pd.DataFrame = pd.read_csv(fnm, sep = "\t",
                      header = None,
                      index_col = False)
    bed.columns = ['chrom', 'start', 'end', 'name']
    bed.set_index('name', drop = False, inplace = True)
    return bed

proj_dir='/SN_aging'
# * meta
## lastest metadata version: v9.6
atacmetafnm = os.path.join(proj_dir,
                           "rds/RNA/cluster/integration/meta/SN.allen.integration.meta.final.barcode.csv")
#atacmeta: pd.DataFrame = pyreadr.read_r(atacmetafnm)[None]
atacmeta: pd.DataFrame = pd.read_csv(atacmetafnm, sep=",", header=0, names=None)
atacmeta['barcode2'] = [f"{i}.{j}" for i, j in zip(atacmeta['orig.ident'].to_list(), atacmeta['barcode'].to_list())]
atacmeta.set_index("barcode2", inplace = True)

samples: List[str] = atacmeta["orig.ident"].unique().tolist()

# * get fraction of cells for peaks
# * rnd peak
adata_rnd = sa2.read(os.path.join(proj_dir, 
                 "rds/ATAC/after_integra/combine/ATAC.shuffle.combined.pmat.adata.h5ad"),
                 backed = "r")
adata_rnd = adata_rnd.to_memory()
groupby_rnd = atacmeta.loc[adata_rnd.obs_names, "subclass_label_id"]
subclasses = groupby_rnd.unique().tolist()

def get_peakfrac_rnd(group: str):
    t = adata_rnd.X[groupby_rnd == group, :]
    r = (t != 0).sum(axis = 0) / t.shape[0]
    return r
# fast, but each will cost 20G
with Pool(10) as p:
    peakfrac_rnd = p.map(get_peakfrac_rnd, subclasses)
peakfrac_rnd_mat = np.concatenate(peakfrac_rnd, axis = 0)
# to pandas
peakfrac_rnd_df = pd.DataFrame(peakfrac_rnd_mat, index = subclasses)
# save to pickle
peakfrac_rnd_df.to_pickle(
    os.path.join(proj_dir, "peak_calling/after_integra/scfilter", "peakfrac_rnd.pkl"))

# * union peak
adata_union = sa2.read(os.path.join(proj_dir, 
                                    "rds/ATAC/after_integra/combine/ATAC.combined.pmat.adata.h5ad"),
                       backed = "r")
adata_union = adata_union.to_memory()

groupby_union = atacmeta.loc[adata_union.obs_names, "subclass_label_id"]
subclasses = groupby_union.unique().tolist()

def get_peakfrac_union(group: str):
    t = adata_union.X[groupby_union == group, :]
    r = (t != 0).sum(axis = 0) / t.shape[0]
    return r
# 2 * each cost 150G
# with additional one with 150G
with Pool(2) as p:
    peakfrac_union = p.map(get_peakfrac_union, subclasses)
peakfrac_union_mat = np.concatenate(peakfrac_union, axis = 0)
# to pandas
peakfrac_union_df = pd.DataFrame(peakfrac_union_mat, index = subclasses)
# set colnames
peaknms = adata_union.var_names.tolist()
peakfrac_union_df.columns = peaknms
# save to pickle
peakfrac_union_df.to_pickle(
    os.path.join(proj_dir, "peak_calling/after_integra/scfilter", "peakfrac_union.pkl"))


# * single-cell level pmat
# when pmat without any description, it means count
barcode2subclass: pd.Series = atacmeta[
    'subclass_label_id'].apply(rename_allensubclass)

# load finalized peaks
final_peak_bed_fnm = os.path.join(
    proj_dir, "peak_calling/after_integra/final/",
    "SN_integra.final.peak.srt.bed")
final_peaks: pd.DataFrame = read_bed4(final_peak_bed_fnm)


# * load AnnData with pmat of union peak
adata_union = sa2.read(
    os.path.join(proj_dir, "rds/ATAC/after_integra/combine/",
                 "ATAC.combined.pmat.adata.h5ad"),
    backed = 'r'
)

# once adata loaded into memory,
# it becomes adata.AnnData, not SnapATAC2 AnnData
# so we can use functions from adata.AnnData
# snap = adata_union.to_memory()
# check if all the final peaks are there.
adata_var_mask = pd.Series(adata_union.var_names).isin(final_peaks.name)
adata_obs_group = barcode2subclass[adata_union.obs_names]
nds = 1000
seed = 0
def myds(data):
    if len(data) <= nds:
        return data
    return(data.sample(n = nds, replace = False, random_state = seed))
adata_obs_group_ds = adata_obs_group.groupby(adata_obs_group).apply(myds)
adata_obs_ds = adata_obs_group_ds.index.get_level_values('barcode2')
out_ds = os.path.join(proj_dir, "rds/ATAC/after_integra/scfilter/SN_integra.final.peak.srt.ds1000.h5ad")
adata_ds = adata_union.subset(obs_indices = adata_obs_ds, var_indices = final_peaks.name, out = out_ds)
adata_ds.obs['subclass_label_id'] = adata_obs_group_ds

# to anndata's AnnData object
ann_ds: anndata.AnnData = adata_ds.to_memory()
adata_ds.close()
out_ann2ds = os.path.join(proj_dir, "rds/ATAC/after_integra/scfilter/SN_integra.final.peak.srt.ds1000.ann.h5ad")
ann_ds.write_h5ad(out_ann2ds)

out_union = os.path.join(proj_dir, "rds/ATAC/after_integra/scfilter/SN_integra.final.peak.srt.pmat.h5ad")
adata_finalpeak = adata_union.subset(var_indices = final_peaks.name, out = out_union)
adata_finalpeak.obs['subclass_label_id'] = adata_obs_group

ann_union = adata_finalpeak.to_memory()
adata_union.close()
adata_finalpeak.close()
out_ann = os.path.join(proj_dir, "rds/ATAC/after_integra/scfilter/SN_integra.final.peak.srt.pmat.ann.h5ad")
#ann_union.obs['celltype'] = atacmeta['subclass_label_id'][ann_union.obs_names]
ann_union.write_h5ad(out_ann)

# * get sample and subclass-level CPM
## ds1000 data
annds: anndata.AnnData = anndata.read_h5ad(os.path.join(proj_dir, "rds/ATAC/after_integra/scfilter/SN_integra.final.peak.srt.ds1000.ann.h5ad"))
npeak = annds.shape[1]
ann_meta = annds.obs
ann_meta['subclass_sample'] = [f"{i}.{j}" for i, j in zip(ann_meta['subclass_label_id'].to_list(), ann_meta['sample'].to_list())] 
#grouped = ann_meta.groupby("subclass_label_id")
grouped = ann_meta.groupby("subclass_sample")
nsc = grouped.ngroups
cnt_pbysc = pd.DataFrame(
    np.zeros( (npeak, nsc) , dtype = np.float64),
    columns = list(grouped.groups.keys()),
    index = annds.var_names)

for group, idx in grouped.indices.items():
    cnt_pbysc[group] = np.ravel(
        annds.X[idx, ].sum(axis = 0, dtype = np.float64))

cnt_pbysc.to_csv(
    os.path.join(proj_dir, "rds/ATAC/after_integra/scfilter", "count_peakBySampleSubclass.ds1000.csv"),
    index = True,
    header = True
)

sum_pofsc = cnt_pbysc.sum(axis = 0)

values = cnt_pbysc.values
factors = sum_pofsc.values

cpm_pbysc = (values / factors ) * 10e6
cpm_pbysc_df = pd.DataFrame(cpm_pbysc, columns = cnt_pbysc.columns)
cpm_pbysc_df.index = cnt_pbysc.index
cpm_pbysc_df.to_csv(
    os.path.join(proj_dir, "rds/ATAC/after_integra/scfilter", "cpm_peakBySampleSubclass.ds1000.csv"),
    index = True,
    header = True
)

# subclass
annds: anndata.AnnData = anndata.read_h5ad(os.path.join(proj_dir, "rds/ATAC/after_integra/scfilter/SN_integra.final.peak.srt.ds1000.ann.h5ad"))
npeak = annds.shape[1]
ann_meta = annds.obs
grouped = ann_meta.groupby("subclass_label_id")
nsc = grouped.ngroups
cnt_pbysc = pd.DataFrame(
    np.zeros( (npeak, nsc) , dtype = np.float64),
    columns = list(grouped.groups.keys()),
    index = annds.var_names)

for group, idx in grouped.indices.items():
    cnt_pbysc[group] = np.ravel(
        annds.X[idx, ].sum(axis = 0, dtype = np.float64))

cnt_pbysc.to_csv(
    os.path.join(proj_dir, "rds/ATAC/after_integra/scfilter", "count_peakBySubclass.ds1000.csv"),
    index = True,
    header = True
)

sum_pofsc = cnt_pbysc.sum(axis = 0)

values = cnt_pbysc.values
factors = sum_pofsc.values

cpm_pbysc = (values / factors ) * 10e6
cpm_pbysc_df = pd.DataFrame(cpm_pbysc, columns = cnt_pbysc.columns)
cpm_pbysc_df.index = cnt_pbysc.index
cpm_pbysc_df.to_csv(
    os.path.join(proj_dir, "rds/ATAC/after_integra/scfilter", "cpm_peakBySubclass.ds1000.csv"),
    index = True,
    header = True
)

# all data
annds: anndata.AnnData = anndata.read_h5ad(os.path.join(proj_dir, "rds/ATAC/after_integra/scfilter/SN_integra.final.peak.srt.pmat.ann.h5ad"))
npeak = annds.shape[1]
ann_meta = annds.obs
ann_meta['subclass_sample'] = [f"{i}.{j}" for i, j in zip(ann_meta['subclass_label_id'].to_list(), ann_meta['sample'].to_list())] 
#grouped = ann_meta.groupby("subclass_label_id")
grouped = ann_meta.groupby("subclass_sample")
nsc = grouped.ngroups
cnt_pbysc = pd.DataFrame(
    np.zeros( (npeak, nsc) , dtype = np.float64),
    columns = list(grouped.groups.keys()),
    index = annds.var_names)

for group, idx in grouped.indices.items():
    cnt_pbysc[group] = np.ravel(
        annds.X[idx, ].sum(axis = 0, dtype = np.float64))

cnt_pbysc.to_csv(
    os.path.join(proj_dir, "rds/ATAC/after_integra/scfilter", "count_peakBySampleSubclass.csv"),
    index = True,
    header = True
)

sum_pofsc = cnt_pbysc.sum(axis = 0)

values = cnt_pbysc.values
factors = sum_pofsc.values

cpm_pbysc = (values / factors ) * 10e6
cpm_pbysc_df = pd.DataFrame(cpm_pbysc, columns = cnt_pbysc.columns)
cpm_pbysc_df.index = cnt_pbysc.index
cpm_pbysc_df.to_csv(
    os.path.join(proj_dir, "rds/ATAC/after_integra/scfilter", "cpm_peakBySampleSubclass.csv"),
    index = True,
    header = True
)

# subclass
annds: anndata.AnnData = anndata.read_h5ad(os.path.join(proj_dir, "rds/ATAC/after_integra/scfilter/SN_integra.final.peak.srt.pmat.ann.h5ad"))
npeak = annds.shape[1]
ann_meta = annds.obs
grouped = ann_meta.groupby("subclass_label_id")
nsc = grouped.ngroups
cnt_pbysc = pd.DataFrame(
    np.zeros( (npeak, nsc) , dtype = np.float64),
    columns = list(grouped.groups.keys()),
    index = annds.var_names)

for group, idx in grouped.indices.items():
    cnt_pbysc[group] = np.ravel(
        annds.X[idx, ].sum(axis = 0, dtype = np.float64))

cnt_pbysc.to_csv(
    os.path.join(proj_dir, "rds/ATAC/after_integra/scfilter", "count_peakBySubclass.csv"),
    index = True,
    header = True
)

sum_pofsc = cnt_pbysc.sum(axis = 0)

values = cnt_pbysc.values
factors = sum_pofsc.values

cpm_pbysc = (values / factors ) * 10e6
cpm_pbysc_df = pd.DataFrame(cpm_pbysc, columns = cnt_pbysc.columns)
cpm_pbysc_df.index = cnt_pbysc.index
cpm_pbysc_df.to_csv(
    os.path.join(proj_dir, "rds/ATAC/after_integra/scfilter", "cpm_peakBySubclass.csv"),
    index = True,
    header = True
)


