#!/usr/bin/env python3
"""
msi_sensor.py
Created on March 30, 2025,
Author: Nick Wlodychak

This module provides functions for parsing single cell data into a format for msi-sensor from
https://github.com/xjtu-omics/msisensor-rna

msi-sensor requires gene expression data with a SampleID column, here we use the cell barcodes as the SampleID.

- First we read in the scanpy dataframe and filter for empty cells.
- Then we replace 0 values with extremely small values as SVM cannot tolerate 0s.
- We verify model compatibility with the expression data, and then subset the genes to only those required for the model.
- This is then written to a csv file that can be finally fed into msi-sensor.

input file format: 10x expression data
output file format: .csv containing MSI predictions

:copyright: (c) 2025 by Nick Wlodychak
"""

import argparse
import subprocess
import pandas as pd
import os
from pathlib import Path
import scanpy as sc
from sklearn.impute import KNNImputer
from sklearn.impute import SimpleImputer
import numpy as np
import pickle


def load_data_scanpy(expression_data: Path) -> sc.AnnData:
    """
    This loads the gene expression data into scanpy and doe the normalization and QC steps.
    :param expression_data: Where the 10x expression data lives
    :return: patient_filtered: sc.AnnData opject with the filtered genes
    """

    print("Loading expression data...")
    rna_adata = sc.read_10x_mtx(expression_data)

    # patient id extraction
    rna_adata.obs['patient'] = rna_adata.obs_names.str.split('_').str[0]
    print("Unique patients:", rna_adata.obs['patient'].unique())

    # masking by patients
    patient = ["CID3586", "CID3838", "CDI3921", "CID3941", "CID3946", "CID3948", "CID3963", "CID4398", "CID44971"]
    # patient = ["CID3586"]
    patient_mask = rna_adata.obs['patient'].isin(patient)
    patient_filtered = rna_adata[patient_mask].copy()

    # basic qc to filter out genes / cells
    sc.pp.filter_cells(patient_filtered, min_counts = 200)  # Filter out empty droplets
    # sc.pp.filter_genes(patient_filtered, min_cells = 1)  # Filter rarely expressed genes
    # sc.pp.normalize_total(patient_filtered, target_sum = 1e4)
    # sc.pp.log1p(patient_filtered)

    print("Expression data filtered!")
    return patient_filtered


def prepare_msisensor_input(rna_data: sc.AnnData, required_genes: set, outpath: Path) -> Path:
    """
    This function parses the 10x adata table and generates a msisensor input
    :param rna_data: the adata table to be passed to MSI-sensor
    :param required_genes: genes needed for the msi-sensor model
    :param outpath: The output file location - ready for MSI-sensor
    :param impute_method: Method to impute missing data (0: mean replacement, 1: fill with 0, 2: KNN)
    :return: output_path: where the msisensor input file is
    """
    outfile = "msisensor_input.csv"
    output_path = Path(outpath, outfile)

    # convert to dataframe
    df = rna_data.to_df().astype(np.float64)

    # santitaztion
    df = df.replace([np.inf, -np.inf], np.nan).fillna(1e-10)
    df = df.replace(0, 1e-10)  # replace 0s for svm
    df = df.filter(items = required_genes)

    # formatting
    df.index.name = 'SampleID'
    df.reset_index().to_csv(output_path, index = False, float_format = '%.16f')

    #
    # print("Converting expression data into MSI...")
    # expr_matrix = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X
    # expr_matrix = expr_matrix.astype(np.float64)
    #
    # print(f"Matrix dtype: {expr_matrix.X.dtype}")
    # print(f"Sparsity: {100 * (expr_matrix.X.nnz / np.prod(expr_matrix.X.shape)):.2f}%")
    #
    # # Sanitize matrix before DataFrame creation
    # expr_matrix = np.nan_to_num(expr_matrix, nan = 0, posinf = 0, neginf = 0)
    # print("NaNs in expr_matrix:", np.isnan(expr_matrix).sum())
    # gene_names = adata.var_names
    #
    # # parse adata based on the cell barcodes
    # df = pd.DataFrame(expr_matrix,
    #                   index = adata.obs_names,
    #                   columns = gene_names).astype(np.float64)
    # nan_count = df.isna().sum().sum()
    # print("NaNs in expr_matrix:", nan_count)
    # try:
    #     if nan_count > 0:
    #         print(f"Found {nan_count} NaN values. Imputing...")
    #         if impute_method == 0:
    #             # fill nan with column means (gene-wise)
    #             df = df.fillna(df.mean())
    #
    #         elif impute_method == 1:
    #             # fill na with 0
    #             df = df.fillna(0)
    #
    #         elif impute_method == 2:
    #             # impute with knn (intensive computationally)
    #             imputer = KNNImputer(n_neighbors = 5)
    #             expr_matrix_imputed = imputer.fit_transform(expr_matrix)
    #             df = pd.DataFrame(expr_matrix_imputed, index = adata.obs_names, columns = gene_names)
    #
    #         elif impute_method == 3:
    #             imputer = SimpleImputer(missing_values=np.nan, strategy='mean')
    #             expr_matrix_imputed = imputer.fit_transform(expr_matrix)
    #             df = pd.DataFrame(expr_matrix_imputed, index = adata.obs_names, columns = gene_names)
    #         else:
    #             raise ValueError(f"Imputation method Unknown!")
    #
    # except KeyError:
    #     print(f"Imputation method Error!")
    #
    # df.index.name = 'SampleID'
    # df = df.reset_index()
    #
    # # Create output path
    # output_path = Path(outpath, outfile)
    # df.to_csv(output_path, index = False, na_rep = '0')
    print(f"Saved MSIsensor-RNA input to {output_path}")

    return output_path


def run_msi_sensor(input: Path, outpath: Path, model: Path):
    """
    Runs the msi sensor pipeline on the input file
    :param input: a csv formatted file
    :param output_path: the msi sensor output directory
    :return: Nothing
    """

    outfile = "msisensor_output.csv"
    print("Initializing msi sensor pipeline...")
    try:
        output_path = Path(outpath, outfile)
        print("Running MSI-sensor pipeline")
        call = f"msisensor-rna detection -i {input} -o {output_path} -m {model} -d TRUE"
        subprocess.call(call, shell = True)
    except Exception as e:
        print(f"An error occurred while processing msi sensor! {e}")
        exit(1)
    print(f"msi sensor pipeline complete!\n Data at: {output_path}")


def verify_model_compatibility(model: Path, rna_data: sc.AnnData) -> set:
    """
    Check model compatibility with adata frame
    :param model: path to the pickle mode
    :param rna_data: adata read into memory
    :return: required_genes: list of genes needed for SVM
    """

    with open(model, 'rb') as f:
        out_model = pickle.load(f)
    required_genes = out_model.feature_names_in_
    missing = set(required_genes) - set(adata.var_names) # set genes from both frames
    if missing:
        raise ValueError(f"Missing {len(missing)} model-required genes")
    else:
        print(f"Model compatible with {len(required_genes)} Genes")
    return required_genes


if __name__ == "__main__":
    params_dict = {
        "rna_data_path": Path("data/SCP1039/expression"),
        "outpath"      : Path("msisensor_pipeline"),
        "model"        : Path("msisensor-rna/model/TCGA.MSIPopular.model.pkl")
    }

    # housekeeping
    os.chdir(Path(__file__).resolve().parent.parent)
    os.makedirs(params_dict["outpath"], exist_ok = True)

    # main pipeline
    adata = load_data_scanpy(expression_data = params_dict["rna_data_path"])
    required_genes = verify_model_compatibility(model = params_dict["model"],
                                                rna_data = adata)
    msi_input = prepare_msisensor_input(rna_data = adata,
                                        required_genes = required_genes,
                                        outpath = params_dict["outpath"])
    run_msi_sensor(input = msi_input,
                   outpath = params_dict["outpath"],
                   model = params_dict["model"])
