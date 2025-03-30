import subprocess
import scipy
import pandas as pd
import os
from pathlib import Path
import scanpy as sc


def load_data_scanpy(expression_data) -> sc.AnnData:
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
    patient = ["CID3586", "CID3838", "CDI3921", "CID3941", "CID3946", "CID3948", "CID3963"]
    patient_mask = rna_adata.obs['patient'].isin(patient)
    patient_filtered = rna_adata[patient_mask].copy()

    # basic qc to filter out genes / cells
    sc.pp.filter_cells(patient_filtered, min_counts = 200)  # Filter out empty droplets
    # sc.pp.filter_genes(patient_filtered, min_cells = 10)  # Filter rarely expressed genes
    print("Expression data filtered!")
    return patient_filtered


def prepare_msisensor_input(adata, outfile = "msisensor_input.csv") -> Path:
    """
    This function parses the 10x adata table and generates a msisensor input
    :param adata: the adata table to be passed to MSI-sensor
    :param outfile: The output csv file ready for MSI-sensor
    :return:
    """
    print("Converting expression data into MSI...")
    expr_matrix = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X
    gene_names = adata.var_names

    # parse adata based on the cell barcodes
    df = pd.DataFrame(expr_matrix,
                      index = adata.obs_names,
                      columns = gene_names)

    df.index.name = 'SampleID'
    df = df.reset_index()

    # Create output path
    output_path = Path(outpath, outfile)
    df.to_csv(output_path, index = False)
    print(f"Saved MSIsensor-RNA input to {output_path}")

    return output_path


def run_msi_sensor(input, outfile = "msisensor_output"):
    """
    Runs the msi sensor pipeline on the input file
    :param input: a csv formatted file
    :param output_path: the msi sensor output directory
    :return:
    """
    print("Initializing msi sensor pipeline...")
    try:
        output_path = Path(outpath, outfile)
        print("Running MSI-sensor pipeline")
        output = Path(output_path, "msi_predictions.csv")
        call = f"msisensor-rna detection -i {input} -m msisensor-rna/model/TCGA.MSIPopular.model.pkl -o {output}"
        subprocess.call(call, shell = True)
    except Exception as e:
        print("An error occurred while processing msi sensor!")
        exit(1)
    print(f"msi sensor pipeline complete!\n Data at: {output}")

if __name__ == "__main__":
    os.chdir(Path(__file__).resolve().parent.parent)
    outpath = "msisensor_pipeline"
    os.makedirs(outpath, exist_ok = True)
    rna_data_path = "data/SCP1039/expression"

    adata = load_data_scanpy(rna_data_path)
    msi_input = prepare_msisensor_input(adata)
    # msi_input = Path(outpath, "msisensor_input.csv")
    run_msi_sensor(msi_input)
