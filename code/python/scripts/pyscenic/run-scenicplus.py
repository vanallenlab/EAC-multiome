import dill
import scanpy as sc
import pandas as pd
import numpy as np
import os
import warnings

warnings.filterwarnings("ignore")

import pyranges as pr

# Set stderr to null to avoid strange messages from ray
import sys

_stderr = sys.stderr
null = open(os.devnull, "wb")
import pathlib as pl
import pickle

from scenicplus.scenicplus_class import create_SCENICPLUS_object
from scenicplus.preprocessing.filtering import *
from scenicplus.cistromes import *
from scenicplus.enhancer_to_gene import (
    get_search_space,
    calculate_regions_to_genes_relationships,
    GBM_KWARGS,
)
from scenicplus.TF_to_gene import *
from scenicplus.grn_builder.gsea_approach import build_grn

from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser()
    parser.add_argument("-sn", "--samplename", type=str, required=True)
    parser.add_argument("-n", "--ncpus", type=int, required=True)
    return parser.parse_args()


def main():
    args = get_args()
    print(f"Running on sample {args.samplename} with n_cpus = {args.ncpus}")

    sample_name = args.samplename
    ncpus = args.ncpus

    basedir = pl.Path("/add/path/here")
    work_dir = pl.Path("/add/path/here")
    tmp_dir = "/add/path/here"
    annot_dir = pl.Path("/add/path/here")

    refined_annotations = pd.read_csv(
        basedir / "refined_wCNMF_programs_and_sampleid.csv", index_col=0
    )

    datadir = basedir / "full_scrna_data"
    adata = sc.read_h5ad(datadir / "full_cohort.h5ad")

    adata.obs = pd.concat(
        [
            adata.obs,
            refined_annotations[
                ["refined_annotation", "refined_wcancer", "highlevel_annotation"]
            ],
        ],
        axis=1,
    )

    treatment_mapping = {
        "Neoadjuvant CROSS": "Neoadj. chemo",
        "Neoadjuvent carboplatin": "Neoadj. chemo",
    }

    clinical = pd.read_csv(basedir / "EAC_clinical_info.csv", index_col=0)

    metastatic = (clinical["Tumor?"] == "Yes ") & (
        clinical["Site"].str.contains("metastasis")
    )
    metastatic.name = "Metastatic?"
    clinical["Metastatic?"] = metastatic

    clinical["Location"] = clinical["Site"].replace(
        {"GEJ": "Esophagus/GEJ", "Esophagus": "Esophagus/GEJ"}
    )
    clinical["Location"][clinical["Location"].str.contains("Liver")] = "Liver"
    clinical["Location"][clinical["Location"].str.contains("Adrenal")] = "Adrenal gland"
    clinical["Location"][clinical["Location"].str.contains("Peritoneal")] = "Peritoneum"

    clinical["Stage"] = clinical["Grade/stage"].replace(
        {
            "Stage IV ": "IV",
            "Stage IV": "IV",
            "Moderately differentiated; ypT1aN0": "I",
            "Moderately differentiated; pT1aN0": "I",
            "Poorly differentiated; ypT2N0": "II",
            "Presented with stage III became stage IV during esophagectomy when pleural metastases were identified": "III/IV",
        }
    )

    clinical["Treatment"] = [
        "Neoadj. chemo",
        "None",
        "Neoadj. chemo + ICI + RT",
        "None",
        "None",
        "Chemo + HER2 targeted + ICI",
        "Neoadj. chemo + HER2 targeted",
        "Neoadj. chemo + ICI",
        "None",
        "Neoadj. chemo + VEGFR2i",
    ]

    clinical["HER2 status"] = clinical["HER2"].replace({"HER 2 1+": "1+/equivocal"})

    clinical = clinical.sort_values(by=["Tumor?", "Metastatic?", "Location"])

    cistopic_obj = dill.load(open(work_dir / sample_name / "cistopic_obj.pkl", "rb"))
    menr = dill.load(open(work_dir / sample_name / "motifs" / "menr.pkl", "rb"))

    middle_annotations = cistopic_obj.cell_data.copy()

    middle_annotations.highlevel_annotation[
        middle_annotations.highlevel_annotation == "Carcinoma"
    ] = np.nan

    middle_annotations = middle_annotations.highlevel_annotation.fillna(
        middle_annotations["refined_wcancer"]
    )

    smallgroups = middle_annotations.value_counts()[
        (middle_annotations.value_counts()) < 20
    ].index
    smallgroup_mapping = {k: "Other" for k in smallgroups}

    middle_annotations = middle_annotations.replace(smallgroup_mapping)

    cistopic_obj.cell_data["highlevel_wcancer"] = middle_annotations.str.replace(
        "/", "_"
    )

    adata.X = adata.layers["counts"].copy()

    sampleadata = adata[adata.obs.sample_id == sample_name].copy()

    idx = sampleadata.obs.index.str.split("-").str[:2]
    sampleadata.obs.index = ["-".join(ix) for ix in idx]

    scplus_obj = create_SCENICPLUS_object(
        GEX_anndata=sampleadata,
        cisTopic_obj=cistopic_obj,
        menr=menr,
        bc_transform_func=lambda x: f"{x}-{sample_name}",  # function to convert scATAC-seq barcodes to scRNA-seq ones
    )
    scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())

    # only keep the first two columns of the PCA embedding in order to be able to visualize this in SCope
    scplus_obj.dr_cell["GEX_X_pca"] = scplus_obj.dr_cell["GEX_X_pca"].iloc[:, 0:2]
    scplus_obj.dr_cell["GEX_X_umap"] = scplus_obj.dr_cell["GEX_X_umap"].iloc[:, 0:2]

    filter_genes(scplus_obj, min_pct=0.1)
    filter_regions(scplus_obj, min_pct=0.1)

    merge_cistromes(scplus_obj)

    get_search_space(
        scplus_obj,
        biomart_host=None,
        species=None,
        assembly=None,
        pr_annot=pr.PyRanges(pd.read_csv(annot_dir / "annot_ensembl.csv")),
        pr_chromsizes=pr.PyRanges(pd.read_csv(annot_dir / "chromsizes_ensembl.csv")),
        upstream=[1000, 150000],
        downstream=[1000, 150000],
    )
    tmp_dir = "/add/path/here"
    calculate_regions_to_genes_relationships(
        scplus_obj,
        ray_n_cpu=ncpus,
        _temp_dir=tmp_dir,
        importance_scoring_method="GBM",
        importance_scoring_kwargs=GBM_KWARGS,
    )

    print("Saving scplus_obj")
    with open(work_dir / sample_name / "scplus_obj.pkl", "wb") as f:
        pickle.dump(scplus_obj, f)

    infile = open(work_dir / sample_name / "scplus_obj.pkl", "rb")
    scplus_obj = pickle.load(infile)
    infile.close()

    tmp_dir = "/add/path/here"
    tf_file = "/add/path/here/utoronto_human_tfs_v_1.01.txt"
    calculate_TFs_to_genes_relationships(
        scplus_obj,
        tf_file=tf_file,
        ray_n_cpu=ncpus,
        method="GBM",
        _temp_dir=tmp_dir,
        key="TF2G_adj",
    )

    print("Saving scplus_obj")
    with open(work_dir / sample_name / "scplus_obj.pkl", "wb") as f:
        pickle.dump(scplus_obj, f)

    infile = open(work_dir / sample_name / "scplus_obj.pkl", "rb")
    scplus_obj = pickle.load(infile)
    infile.close()

    tmp_dir = "/add/path/here"
    build_grn(
        scplus_obj,
        min_target_genes=10,
        adj_pval_thr=1,
        min_regions_per_gene=0,
        quantiles=(0.85, 0.90, 0.95),
        top_n_regionTogenes_per_gene=(5, 10, 15),
        top_n_regionTogenes_per_region=(),
        binarize_using_basc=True,
        rho_dichotomize_tf2g=True,
        rho_dichotomize_r2g=True,
        rho_dichotomize_eregulon=True,
        rho_threshold=0.05,
        keep_extended_motif_annot=True,
        merge_eRegulons=True,
        order_regions_to_genes_by="importance",
        order_TFs_to_genes_by="importance",
        key_added="eRegulons_importance",
        cistromes_key="Unfiltered",
        disable_tqdm=True,  # If running in notebook, set to True
        ray_n_cpu=ncpus,
        _temp_dir=tmp_dir,
    )

    print("Saving scplus_obj")
    with open(work_dir / sample_name / "scplus_obj.pkl", "wb") as f:
        dill.dump(scplus_obj, f)


if __name__ == "__main__":
    main()
