import pandas as pd
import numpy as np
import os
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from goatools import obo_parser
from goatools.associations import read_gaf
from goatools.go_enrichment import GOEnrichmentStudy
import plotly.express as px
import matplotlib.pyplot as plt

# path to metadata files
metadata_dir = os.path.join(os.path.dirname(os.getcwd()),'metadata')
# Define file paths
GO_OBO = os.path.join(metadata_dir, "go-basic.obo")
GAF_FILE = os.path.join(metadata_dir, "ecocyc.gaf")  # Replace with your GAF file
GTF_FILE = os.path.join(metadata_dir, "genomic.gtf")

def get_ID_conversion(path_to_gtf):
    df = pd.read_csv(path_to_gtf, sep='\t', comment='#', header=None)
    df.columns = ['seqid', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    df = df[df['feature'] == 'CDS']
    # create dictionary to store gene id and gene name
    gene_id_name = {}
    for index, row in df.iterrows():
        gene_id = row['attributes'].split(';')[2].split(' ')[2].replace('"', '')
        # remove 'GeneID:' from gene_id
        gene_id = gene_id.split(':')[-1]
        gene_name = row['attributes'].split(';')[4].split(' ')[2].replace('"', '')
        gene_id_name[gene_name.lower()] = gene_id

    return gene_id_name


def get_lfc_thresh(deseq_df, target_size, fold, p_val_thresh=0.01):
    """
    Function to get the log fold change threshold for a given p-value threshold and target number of DEGs
    :param deseq_df: DataFrame containing DESeq2 results
    :param target_size: target number of DEGs
    :param p_val_thresh: p-value threshold
    :return: log fold change threshold
    """
    df_sorted = deseq_df.copy()
    # filter out genes that do not pass the p-value threshold
    df_sorted = df_sorted[df_sorted['padj'] < p_val_thresh]
    # Sort the DESeq2 results by LFC
    df_sorted = df_sorted.sort_values('log2FoldChange', ascending=False)
    # Get the LFC threshold for the target number of DEGs
    if fold == 'up':
        lfc_thresh = np.maximum(df_sorted.iloc[target_size]['log2FoldChange'],1)
    elif fold == 'down':
        lfc_thresh = np.minimum(df_sorted.iloc[-target_size]['log2FoldChange'],-1)
    else:
        raise ValueError("fold must be 'up' or 'down'")
    return lfc_thresh


def remove_unidentified_genes(deseq_df, gene_id_name):
    """
    Function to remove genes that are not identified
    :param df: DataFrame containing the gene names
    :return: DataFrame with the unidentified genes removed
    """
    # Remove genes that are not identified
    ind = [val.lower() in gene_id_name.keys() for val in deseq_df.index.values]
    df_filtered = deseq_df.iloc[ind]
    return df_filtered


def run_deseq(file_path, metadata_path, output_dir, contrast):
    """
    Function to run DESeq2 analysis
    :param file_path: path to the count data
    :param metadata_path: path to the metadata file
    :param output_dir: path to save the DESeq2 results
    :param contrast: contrast for the DESeq2 analysis: e.g., ["condition", "casp-t2", "exp-t0"]
    :return:
    """

    count_matrix = pd.read_csv(file_path, index_col=0, header=0)
    count_matrix = count_matrix.T
    # Load your count matrix and metadata
    # Assume 'count_matrix.csv' is your RNA-seq count data
    # and 'metadata.csv' contains the sample info with a 'condition' column
    metadata = pd.read_csv(metadata_path, index_col=0, header=0)  # contains sample info (e.g., conditions)
    # Ensure that the samples in the count matrix are in the same order as in the metadata
    assert all(count_matrix.index == metadata.index), "Mismatch between count matrix columns and metadata index"

    # Initialize the DESeq2 dataset
    dds = DeseqDataSet(counts=count_matrix, metadata=metadata, design_factors="condition")

    # Preprocess data: normalization and dispersion estimation
    dds.deseq2()

    # Run the differential expression analysis
    dds_results = DeseqStats(dds, contrast=contrast)

    # Compute the results
    dds_results.summary()  # Summarize the results

    # Extract the DEG table
    deg_results = dds_results.results_df

    # Save the DEG results to a CSV file
    deg_results.to_csv(os.path.join(output_dir, f"deseq2_results_{contrast[1]}_vs_{contrast[2]}.csv"))
    print("Differential expression analysis completed and results saved!")


def plot_deseq_results(deseq_file, lfc_cutoff=1, p_cutoff=0.01, output_dir=None):
    # plot volcano plot of deseq2 results
    deg_results = pd.read_csv(deseq_file, index_col=0)
    # Plot the volcano plot
    # Add a new column for -log10(padj)
    deg_results['-log10(padj)'] = -np.log10(deg_results['padj'].replace(0, 1e-300))  # Avoid log(0)
    # select significant genes
    significance = (deg_results['padj'] < p_cutoff) & (np.abs(deg_results['log2FoldChange']) > lfc_cutoff)
    deg_results['significant'] = significance

    # Create the volcano plot
    fig = px.scatter(
        deg_results,
        x='log2FoldChange',
        y='-log10(padj)',
        title='Volcano Plot',
        color='significant',
        color_discrete_sequence=['red', 'grey'],
        labels={'log2FoldChange': 'Log2 Fold Change', '-log10(padj)': '-Log10 Adjusted P-value'}
    )
    # background color
    # Highlight significant genes
    fig.add_hline(y=-np.log10(p_cutoff), line_dash='dash', line_color='black', annotation_text='P-value cutoff',
                  annotation_position='bottom right')
    fig.add_vline(x=lfc_cutoff, line_dash='dash', line_color='black', annotation_text='Log2FC cutoff',
                  annotation_position='bottom right')
    fig.add_vline(x=-lfc_cutoff, line_dash='dash', line_color='black')
    # color points below the cutoff in red

    # set fig size
    fig.update_layout(width=800, height=800)
    # set xlimits:
    fig.update_xaxes(range=[-10, 10])
    # Show the gene names on hover

    fig.update_traces(text=deg_results.index,
                      hovertemplate='Gene: %{text}<br>Log2 Fold Change: %{x}<br>-Log10 Adjusted P-value: %{y}')
    # Save the plot, if an output directory is provided
    if output_dir:
        fig.write_image(os.path.join(output_dir, 'volcano_plot.png'))
    # Show the plot
    fig.show()


def run_go_enrichment(DEG_file, p_cutoff, target_size, fold, output_dir):
    # path to metadata files
    # Load GO ontology
    go_dag = obo_parser.GODag(GO_OBO)
    # Load Gene Association File
    geneid2gos = read_gaf(GAF_FILE)
    # load gene id conversion
    gene_id_name = get_ID_conversion(GTF_FILE)
    # load DEG results
    deg_results = pd.read_csv(DEG_file, index_col=0)
    # remove unidentified genes
    deg_results = remove_unidentified_genes(deg_results, gene_id_name)
    # Read DEG and background gene lists
    # get the log fold change threshold
    lfc_cutoff = get_lfc_thresh(deg_results, target_size, fold, p_val_thresh=p_cutoff)
    if fold == 'up':
        deg_genes = set(
            deg_results.index[(deg_results['padj'] < p_cutoff) & (deg_results['log2FoldChange'] > lfc_cutoff)])
    elif fold == 'down':
        deg_genes = set(
            deg_results.index[(deg_results['padj'] < p_cutoff) & (deg_results['log2FoldChange'] < lfc_cutoff)])

    background_genes = set(deg_results.index)
    # translate to gene ID:
    deg_genes_ID = []
    for gene in deg_genes:
        try:
            deg_genes_ID.append(gene_id_name[gene.lower()])
        except:
            pass
    background_genes_ID = []
    for gene in background_genes:
        try:
            background_genes_ID.append(gene_id_name[gene.lower()])
        except:
            pass
    # Initialize GOEnrichmentStudy object
    goeaobj = GOEnrichmentStudy(
        background_genes_ID,  # List of background genes
        geneid2gos,  # geneid/GO associations
        go_dag,  # GO DAG
        propagate_counts=False,
        alpha=0.05,  # Significance cut-off
        methods=['fdr_bh']  # Multiple testing correction
    )

    # Run GO enrichment analysis
    goea_results = goeaobj.run_study(deg_genes_ID)
    # Filter significant results
    significant_results = [r for r in goea_results if r.p_fdr_bh < 0.05]

    # Convert to DataFrame for easier handling
    results_df = pd.DataFrame({
        "GO_ID": [r.GO for r in significant_results],
        "Term": [r.name for r in significant_results],
        "Category": [r.NS for r in significant_results],
        "p-value": [r.p_uncorrected for r in significant_results],
        "FDR": [r.p_fdr_bh for r in significant_results],
        "Genes": [",".join(r.study_items) for r in significant_results],
        "Fold Enrichment": [r.enrichment for r in significant_results],
        "Depth": [r.depth for r in significant_results],
        "Ratio in Study": [r.ratio_in_study for r in significant_results],
        "Ratio in Population": [r.ratio_in_pop for r in significant_results],
        "Parents": [r.goterm._parents for r in significant_results],
    })
    results_df.to_csv(os.path.join(output_dir, f"GO_enrichment_{DEG_file}.csv"),
                      index=False)
    # Print top 10 enriched GO terms
    print(results_df.head(10))


def run_go_single_cell(DEG_file, cluster, p_cutoff, fold, output_dir):
    # path to metadata files
    # Load GO ontology
    go_dag = obo_parser.GODag(GO_OBO)
    # Load Gene Association File
    geneid2gos = read_gaf(GAF_FILE)
    # load gene id conversion
    gene_id_name = get_ID_conversion(GTF_FILE)
    # load DEG results
    deg_results = pd.read_excel(os.path.join(os.path.dirname(os.getcwd()),'scanpy','marker_genes_per_cluster_shx_scaled.xlsx'), sheet_name=cluster, index_col=0, header=0)
    # remove unidentified genes
    deg_results = remove_unidentified_genes(deg_results, gene_id_name)
    # Read DEG and background gene lists
    # get the log fold change threshold
    lfc_cutoff = 0
    if fold == 'up':
        deg_genes = set(
            deg_results.index[(deg_results['adjusted_pval'] < p_cutoff) & (deg_results['log2_fold_change'] > lfc_cutoff)])
    elif fold == 'down':
        deg_genes = set(
            deg_results.index[(deg_results['adjusted_pval'] < p_cutoff) & (deg_results['log2_fold_change'] < lfc_cutoff)])

    background_genes = set(deg_results.index)
    # translate to gene ID:
    deg_genes_ID = []
    for gene in deg_genes:
        try:
            deg_genes_ID.append(gene_id_name[gene.lower()])
        except:
            pass
    background_genes_ID = []
    for gene in background_genes:
        try:
            background_genes_ID.append(gene_id_name[gene.lower()])
        except:
            pass
    # Initialize GOEnrichmentStudy object
    goeaobj = GOEnrichmentStudy(
        background_genes_ID,  # List of background genes
        geneid2gos,  # geneid/GO associations
        go_dag,  # GO DAG
        propagate_counts=False,
        alpha=0.05,  # Significance cut-off
        methods=['fdr_bh']  # Multiple testing correction
    )

    # Run GO enrichment analysis
    goea_results = goeaobj.run_study(deg_genes_ID)
    # Filter significant results
    significant_results = [r for r in goea_results if r.p_fdr_bh < 0.05]

    # Convert to DataFrame for easier handling
    results_df = pd.DataFrame({
        "GO_ID": [r.GO for r in significant_results],
        "Term": [r.name for r in significant_results],
        "Category": [r.NS for r in significant_results],
        "p-value": [r.p_uncorrected for r in significant_results],
        "FDR": [r.p_fdr_bh for r in significant_results],
        "Genes": [",".join(r.study_items) for r in significant_results],
        "Fold Enrichment": [r.enrichment for r in significant_results],
        "Depth": [r.depth for r in significant_results],
        "Ratio in Study": [r.ratio_in_study for r in significant_results],
        "Ratio in Population": [r.ratio_in_pop for r in significant_results],
        "Parents": [r.goterm._parents for r in significant_results],
    })
    results_df.to_csv(os.path.join(output_dir, f"GO_enrichment_cluster_{cluster}.csv"),
                      index=False)
    # Print top 10 enriched GO terms
    print(results_df.head(10))