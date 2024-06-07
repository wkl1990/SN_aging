import os
import sys
from typing import Tuple, List, Dict
import pandas as pd
import snapatac2 as sa2

# load all the h5ads
prefix = "data/final"
suffix = "final_cells.h5ad"

sample: List[str] = [
    "2m_rep1",
    "2m_rep2",
    "6m_rep1",
    "6m_rep2",
    "12m_rep1",
    "12m_rep2",
    "18m_rep1",
    "18m_rep2"
]

sample2fnm: List[Tuple[str, str]] = [
    (i, f"{prefix}/{i}.{suffix}") for i in sample]

# * add tile to ann
sample2ann: Dict[str, sa2.AnnData] = {
    k: sa2.read(v, backed = 'r+') for k, v in sample2fnm
}

# add X for each h5ad
for _, v in sample2ann.items():
    print(f"add tile to {v}")
    sa2.pp.add_tile_matrix(adata = v, bin_size = 500,
                           inplace = True,
                           exclude_chroms=['chrM', 'chrY', 'M', 'Y'],
                           counting_strategy= 'insertion')
# save the data
for _, v in sample2ann.items():
    v.close()

# * merge all the h5ads
all_ann = sa2.AnnDataSet(
    adatas=sample2fnm,
    filename = "data/all_sa2_dataset.h5ad",
    add_key = 'sample'
)

all_ann = all_ann.to_adata(copy_x = True,
                           file = "data/all_sa2_ann.h5ad")
all_ann.close()

# * load cell meta data
cellmeta: pd.DataFrame = pd.read_csv(
    "data/SN.allen.integration.meta.final.barcode.csv",
    sep = ",",
    header=0
)
barcodes = [f'{r["sampleID"]}_{r["barcode"]}'
            for _, r in cellmeta.iterrows()]
cellmeta.insert(0, 'ubarcode', barcodes)
cellmeta.set_index("ubarcode", drop = False, inplace = True)

# * embedding
all_ann = sa2.read(filename = "data/all_sa2_ann.h5ad", backed = 'r+')
sa2.pp.select_features(adata = all_ann, n_features = 500000,
                       filter_lower_quantile=0.005,
                       filter_upper_quantile=0.005,
                       inplace = True)
sa2.tl.spectral(adata = all_ann, n_comps = 30,
                features = 'selected',
                inplace = True,
                weighted_by_sd = True)
sa2.tl.umap(adata = all_ann, n_comps=2,
            use_rep = 'X_spectral',
            key_added = "umap_no_batchcrctn",
            inplace=True)

# * add sample and cluster information
old_obs_names = all_ann.obs_names
uniq_obs_names = [f"{s}_{b}" for s, b in zip(all_ann.obs['sample'],
                                                all_ann.obs_names)]
all_ann.obs_names = uniq_obs_names
clusters = cellmeta.loc[all_ann.obs_names].seurat_clusters
all_ann.obs['cluster'] = [f"cl-{i}" for i in clusters]
all_ann.obs['class_label'] = cellmeta.loc[all_ann.obs_names].class_label
all_ann.obs['subclass_label'] = cellmeta.loc[all_ann.obs_names].subclass_label


# * umap plot
sa2.pl.umap(adata = all_ann, use_rep = "X_umap_no_batchcrctn",
            color = "cluster",
            out_file = "out/umap_no_batchcrctn.pdf")

sa2.pl.umap(adata = all_ann, use_rep = "X_umap_no_batchcrctn",
            color = "sample",
            out_file = "out/sample_umap_no_batchcrctn.pdf")
sa2.pl.umap(adata = all_ann, use_rep = "X_umap_no_batchcrctn",
            color = "class_label",
            out_file = "out/class_umap_no_batchcrctn.pdf")
sa2.pl.umap(adata = all_ann, use_rep = "X_umap_no_batchcrctn",
            color = "subclass_label",
            out_file = "out/subclass_umap_no_batchcrctn.pdf")

# * perform harmony for batch correction
sa2.pp.harmony(adata = all_ann, batch = "sample",
               use_rep = "X_spectral",
               inplace = True)
sa2.tl.umap(adata = all_ann, n_comps = 2,
            use_rep = "X_spectral_harmony",
            key_added = "umap_batchcrctn",
            inplace = True)

# * plot after batch correction
# sa2.pl.umap(adata = all_ann, use_rep = "X_umap_no_batchcrctn",
#             color = "cluster",
#             out_file = "out/umap_no_batchcrctn.pdf")

sa2.pl.umap(adata = all_ann, use_rep = "X_umap_batchcrctn",
            color = "sample",
            out_file = "out/sample_umap_batchcrctn.pdf")
sa2.pl.umap(adata = all_ann, use_rep = "X_umap_batchcrctn",
            color = "class_label",
            out_file = "out/class_umap_batchcrctn.pdf")
sa2.pl.umap(adata = all_ann, use_rep = "X_umap_batchcrctn",
            color = "subclass_label",
            out_file = "out/subclass_umap_batchcrctn.pdf")

all_ann.close()

# * save data to csv and prepare to plot using R
all_ann = sa2.read(filename = "data/all_sa2_ann.h5ad", backed = 'r')
cell_umap = pd.DataFrame(all_ann.obsm['X_umap_batchcrctn'],
                         index = all_ann.obs_names,
                         columns = ["UMAP1", "UMAP2"])
cell_umap.insert(0, 'barcode', cell_umap.index)
cell_umap.insert(1, 'class_id', cellmeta.loc[cell_umap.index].class_id)
cell_umap.insert(2, 'subclass_id', cellmeta.loc[cell_umap.index].subclass_id)
cell_umap.insert(3, 'class_color', cellmeta.loc[cell_umap.index].class_color)
cell_umap.insert(4, 'subclass_color', cellmeta.loc[cell_umap.index].subclass_color)
cell_umap.to_csv("out/cell_umap_after_batch_correction.csv",
                 sep = ",", header = True, index = False)
all_ann.close()
