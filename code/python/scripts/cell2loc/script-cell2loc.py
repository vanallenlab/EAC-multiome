import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns
import pathlib as pl

import cell2location

from matplotlib import rcParams

from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser()
    parser.add_argument("-p", "--patient", type=str, required=True)
    return parser.parse_args()


def get_preprocessed_sample(
    sample_path: pl.Path, min_counts: int, pct_mt: int, min_cells: int
) -> sc.AnnData:

    adata = sc.read_visium(path=sample_path)

    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    adata.obsm["spatial"] = adata.obsm["spatial"].astype(int)

    sc.pp.filter_cells(adata, min_counts=min_counts)
    adata = adata[adata.obs["pct_counts_mt"] < pct_mt]
    print(f"#cells after MT filter: {adata.n_obs}")
    sc.pp.filter_genes(adata, min_cells=min_cells)

    return adata


def main():

    args = get_args()
    print(f"Running on patient {args.patient}")
    patient_name = args.patient

    spatial_dir = pl.Path(
        "/add/path/here/"
    )

    inf_aver = pd.read_csv(
        "/add/path/here/mean-expression-per-celltype.csv",
        index_col=0,
    )

    results_folder = f"/add/path/here/Cell2Location_results/{patient_name}"

    # create paths and names to results folders for reference regression and cell2location models
    run_name = f"{results_folder}/cell2location_map"

    sample_path = spatial_dir / patient_name

    adata = get_preprocessed_sample(
        sample_path=sample_path, min_counts=5000, pct_mt=30, min_cells=10
    )

    tissue_path = spatial_dir / patient_name / "spatial/tissue_positions_list.csv"
    tissue_position = pd.read_csv(tissue_path, index_col=0)
    tissue_position = tissue_position.loc[adata.obs_names]

    # Set coordinates
    x_array = tissue_position["array_row"].tolist()
    y_array = tissue_position["array_col"].tolist()
    x_pixel = tissue_position["pxl_row_in_fullres"].tolist()
    y_pixel = tissue_position["pxl_col_in_fullres"].tolist()

    x_min, x_max = np.min(x_pixel), np.max(x_pixel)
    y_min, y_max = np.min(y_pixel), np.max(y_pixel)

    adata.obs["sample"] = patient_name

    # find shared genes and subset both anndata and reference signatures
    intersect = np.intersect1d(adata.var_names, inf_aver.index)
    adata = adata[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    # prepare anndata for cell2location model
    cell2location.models.Cell2location.setup_anndata(adata=adata, batch_key="sample")

    # create and train the model
    mod = cell2location.models.Cell2location(
        adata,
        cell_state_df=inf_aver,
        # the expected average cell abundance: tissue-dependent
        # hyper-prior which can be estimated from paired histology:
        N_cells_per_location=5,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection:
        detection_alpha=20,
    )
    mod.view_anndata_setup()

    mod.train(
        max_epochs=20000,
        # train using full data (batch_size=None)
        batch_size=None,
        # use all data points in training because
        # we need to estimate cell abundance at all locations
        train_size=1,
        use_gpu=False,
    )

    # plot ELBO loss history during training, removing first 1000 epochs from the plot
    mod.plot_history(1000)
    plt.legend(labels=["full data training"])

    # Save model
    mod.save(f"{run_name}", overwrite=True)

    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata = mod.export_posterior(
        adata,
        sample_kwargs={
            "num_samples": 1000,
            "batch_size": mod.adata.n_obs,
            "use_gpu": False,
        },
    )

    # Save anndata object with results
    adata_file = f"{run_name}/sp.h5ad"
    adata.write(adata_file)
    adata_file

    mod.plot_QC()

    # Compute expected expression per cell type
    expected_dict = mod.module.model.compute_expected_per_cell_type(
        mod.samples["post_sample_q05"], mod.adata_manager
    )

    # Add to anndata layers
    for i, n in enumerate(mod.factor_names_):
        adata.layers[n] = expected_dict["mu"][i]

    # add 5% quantile, representing confident cell abundance, 'at least this amount is present',
    # to adata.obs with nice names for plotting
    adata.obs[adata.uns["mod"]["factor_names"]] = adata.obsm["q05_cell_abundance_w_sf"]

    # Save anndata object with results
    adata_file = f"{run_name}/sp.h5ad"
    adata.write(adata_file)
    adata_file

    adata.obs[
        ["Carcinoma", "Endothelial", "Fibroblast", "Lymphoid", "Myeloid", "Muscle"]
    ].to_csv(run_name + "/celltype_abundance.csv")


if __name__ == "__main__":

    main()
