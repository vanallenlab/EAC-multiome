import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)
import sys
import os

_stderr = sys.stderr
null = open(os.devnull, "wb")

import pathlib as pl
import scanpy as sc
import pandas as pd
import numpy as np
from tqdm.notebook import tqdm
import pycisTopic

from typing import List, Optional, Union

import matplotlib
import matplotlib.pyplot as plt
import scipy
import scipy.sparse as sparse
import sklearn
from scipy.stats import ranksums

from pycisTopic.cistopic_class import *
from pycisTopic.utils import *
import pickle
from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import *
from pycisTopic.clust_vis import *
from pycisTopic.topic_binarization import *
from pycisTopic.diff_features import *

import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
import pybiomart as pbm

from scenicplus.wrappers.run_pycistarget import run_pycistarget

from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser()
    parser.add_argument("-sn", "--samplename", type=str, required=True)
    parser.add_argument("-n", "--ncpus", type=int, required=True)
    return parser.parse_args()


def main():

    args = get_args()
    print(f"Running on sample {args.samplename} with n_cpus = {args.ncpus}")
    work_dir = pl.Path("/add/path/here/")
    os.makedirs(work_dir, exist_ok=True)

    atacdir = pl.Path("/add/path/here/")
    macsdir = pl.Path("/add/path/here/")

    print("Downloading scRNA-seq data...")
    adata = sc.read_h5ad("/add/path/here/full_cohort.h5ad")

    sample_name = args.samplename

    print(sample_name)

    os.makedirs(work_dir / sample_name, exist_ok=True)

    sampledir = atacdir / sample_name

    cell_data = adata[adata.obs.sample_id == sample_name].obs.copy()
    cell_data["highlevel_celltype"] = cell_data["highlevel_celltype"].astype(
        str
    )  # set data type of the celltype column to str, otherwise the export_pseudobulk function will complain.

    idx = cell_data.index.str.split("-").str[:2]
    cell_data.index = ["-".join(ix) for ix in idx]

    fragments_dict = {sample_name: os.path.join(sampledir, "atac_fragments.tsv.gz")}
    path_to_regions = {
        sample_name: os.path.join(macsdir, sample_name, f"{sample_name}.bed")
    }
    path_to_blacklist = "/add/path/here/hg38-blacklist.v2.bed"
    metadata_bc = {
        sample_name: pd.read_csv(sampledir / "per_barcode_metrics.csv", index_col=0)
    }

    atac_cell_bc = metadata_bc[sample_name][metadata_bc[sample_name].is_cell == 1].index

    key = sample_name
    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments=fragments_dict[key],
        path_to_regions=path_to_regions[key],
        path_to_blacklist=path_to_blacklist,
        metrics=metadata_bc[key],
        valid_bc=list(atac_cell_bc.intersection(cell_data.index)),
        n_cpu=1,
        project=key,
        split_pattern="-",
    )
    cistopic_obj.add_cell_data(cell_data, split_pattern="-")
    print(cistopic_obj)

    pickle.dump(cistopic_obj, open(work_dir / sample_name / "cistopic_obj.pkl", "wb"))

    tmp_dir = "/add/path/here/"
    cistopic_obj = pickle.load(open(work_dir / sample_name / "cistopic_obj.pkl", "rb"))

    models = run_cgs_models(
        cistopic_obj,
        n_topics=[2, 4, 10, 16],
        n_cpu=args.ncpus,
        n_iter=500,
        random_state=555,
        alpha=50,
        alpha_by_topic=True,
        eta=0.1,
        eta_by_topic=False,
        save_path=None,
        _temp_dir=os.path.join(tmp_dir + "ray_spill"),
    )

    os.makedirs(work_dir / sample_name / "models", exist_ok=True)
    pickle.dump(
        models,
        open(work_dir / sample_name / "models" / "models_500_iter_LDA.pkl", "wb"),
    )

    models = pickle.load(
        open(work_dir / sample_name / "models" / "models_500_iter_LDA.pkl", "rb")
    )
    cistopic_obj = pickle.load(open(work_dir / sample_name / "cistopic_obj.pkl", "rb"))
    print("creating topic models...")
    model = evaluate_models(
        models,
        return_model=True,
        metrics=["Arun_2010", "Cao_Juan_2009", "Minmo_2011", "loglikelihood"],
        plot_metrics=False,
    )

    cistopic_obj.add_LDA_model(model)
    pickle.dump(cistopic_obj, open(work_dir / sample_name / "cistopic_obj.pkl", "wb"))

    region_bin_topics_otsu = binarize_topics(cistopic_obj, method="otsu")
    region_bin_topics_top3k = binarize_topics(cistopic_obj, method="ntop", ntop=3000)

    imputed_acc_obj = impute_accessibility(
        cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6
    )
    normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
    variable_regions = find_highly_variable_features(
        normalized_imputed_acc_obj, plot=False
    )
    #

    cistopic_obj.cell_data["highlevel_celltype"] = cistopic_obj.cell_data[
        "highlevel_celltype"
    ].str.replace("/", "_")

    markers_dict = find_diff_features(
        cistopic_obj,
        imputed_acc_obj,
        variable="highlevel_celltype",
        var_features=variable_regions,
        split_pattern="-",
    )

    os.makedirs(work_dir / sample_name / "candidate_enhancers", exist_ok=True)
    pickle.dump(
        region_bin_topics_otsu,
        open(
            work_dir
            / sample_name
            / "candidate_enhancers"
            / "region_bin_topics_otsu.pkl",
            "wb",
        ),
    )
    pickle.dump(
        region_bin_topics_top3k,
        open(
            work_dir
            / sample_name
            / "candidate_enhancers"
            / "region_bin_topics_top3k.pkl",
            "wb",
        ),
    )
    pickle.dump(
        markers_dict,
        open(work_dir / sample_name / "candidate_enhancers" / "markers_dict.pkl", "wb"),
    )

    region_bin_topics_otsu = pickle.load(
        open(
            work_dir
            / sample_name
            / "candidate_enhancers"
            / "region_bin_topics_otsu.pkl",
            "rb",
        )
    )
    region_bin_topics_top3k = pickle.load(
        open(
            work_dir
            / sample_name
            / "candidate_enhancers"
            / "region_bin_topics_top3k.pkl",
            "rb",
        )
    )
    markers_dict = pickle.load(
        open(work_dir / sample_name / "candidate_enhancers" / "markers_dict.pkl", "rb")
    )

    region_sets = {}
    region_sets["topics_otsu"] = {}
    region_sets["topics_top_3"] = {}
    region_sets["DARs"] = {}
    for topic in region_bin_topics_otsu.keys():
        regions = region_bin_topics_otsu[topic].index[
            region_bin_topics_otsu[topic].index.str.startswith("chr")
        ]  # only keep regions on known chromosomes
        if len(regions) > 0:
            region_sets["topics_otsu"][topic] = pr.PyRanges(
                region_names_to_coordinates(regions)
            )
    for topic in region_bin_topics_top3k.keys():
        regions = region_bin_topics_top3k[topic].index[
            region_bin_topics_top3k[topic].index.str.startswith("chr")
        ]  # only keep regions on known chromosomes
        if len(regions) > 0:
            region_sets["topics_top_3"][topic] = pr.PyRanges(
                region_names_to_coordinates(regions)
            )
    for DAR in markers_dict.keys():
        regions = markers_dict[DAR].index[
            markers_dict[DAR].index.str.startswith("chr")
        ]  # only keep regions on known chromosomes
        if len(regions) > 0:
            region_sets["DARs"][DAR] = pr.PyRanges(region_names_to_coordinates(regions))

    # in case there are no DARs, e.g.
    todel = []
    for k, v in region_sets.items():
        if len(v) == 0:
            todel.append(k)
    for k in todel:
        del region_sets[k]

    auxiliarypath = pl.Path("/add/path/here/")
    # corresponds to Screen v10 region-based databases, SCENIC+, #1
    rankings_db = (
        auxiliarypath / "hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
    )
    # corresponds to Screen v10 region-based databases, SCENIC+, #2
    scores_db = auxiliarypath / "hg38_screen_v10_clust.regions_vs_motifs.scores.feather"
    # corresponds to Motif v10 annotation, SCENIC+
    motif_annotation = auxiliarypath / "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

    os.makedirs(work_dir / sample_name / "motifs", exist_ok=True)
    # running with custom species because no connection to internet through cluster
    # however could just be run with species="homo_sapiens" in which case the custom_annot
    # argument doesn't have to be populated with annot_dem.

    # corresponds to Annotation for local pycisTarget run
    annot_dem = pd.read_csv(auxiliarypath / "annot_ensembl.csv", index_col=0)

    run_pycistarget(
        region_sets=region_sets,
        species="custom",
        custom_annot=annot_dem,
        save_path=os.path.join(work_dir, f"{sample_name}/motifs"),
        ctx_db_path=rankings_db,
        dem_db_path=scores_db,
        path_to_motif_annotations=motif_annotation,
        run_without_promoters=False,
        n_cpu=args.ncpus,
        save_partial=True,
        _temp_dir=os.path.join(tmp_dir, "ray_spill"),
        annotation_version="v10nr_clust",
    )


if __name__ == "__main__":

    main()
