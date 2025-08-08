import pandas as pd
import os
import glob

def get_chromosomes():
    """
    Returns a list of chromosomes.
    
    Returns:
        list: List of chromosome names.
    """
    return [f'chr{i}' for i in range(1, 23)] # + ['chrX', 'chrY']

def read_gene_sets(file_path):
    """
    Reads gene sets from a file and returns a DataFrame.
    
    Args:
        file_path (str): Path to the gene sets file.
        
    Returns:
        pd.DataFrame: DataFrame containing gene sets.
    """
    # if ends with obj, read as pickle
    if file_path.endswith('.obj'):
        df = pd.read_pickle(file_path)
    # if ends with csv, read as csv
    elif file_path.endswith('.csv'):
        df = pd.read_csv(file_path)

    # reset index
    df.reset_index(inplace=True)
    df.columns = ['gene_set', 'gene']
    # genes is list of genes, explode it
    df = df.explode('gene')
    return df

def read_setlist(file_path):
    """
    Reads a gene set file in the regenie format and returns a DataFrame.
    
    Args:
        file_path (str): Path to the regenie gene set file.
        
    Returns:
        pd.DataFrame: DataFrame containing gene sets.
    """
    df = pd.read_csv(file_path, sep='\t', header=None)
    df.columns = ['transcript', 'chr', 'pos', 'snp']

    return df

def read_regenie_annotation(file_path):
    """
    Reads an annotation file and returns a DataFrame.
    
    Args:
        file_path (str): Path to the annotation file.
        
    Returns:
        pd.DataFrame: DataFrame containing annotations.
    """
    df = pd.read_csv(file_path, sep='\t', header=None)
    df.columns = ['snp', 'transcript', 'snp_set']
    
    return df

def convert_annotation(annot_filepath, geneset_df):
    """
    Converts an annotation file to a DataFrame with gene sets.
    
    Args:
        annot_filepath (str): Path to the annotation file.
        geneset_df (pd.DataFrame): DataFrame containing gene sets.
        
    Returns:
        pd.DataFrame: DataFrame with gene sets and annotations.
    """
    annot_df = read_regenie_annotation(annot_filepath)
    annot_df_merged = annot_df.merge(geneset_df, on='gene', how='left')

    # report if there are any genes in the annotation that are not in the gene set
    missing_genes = annot_df_merged[annot_df_merged['gene_set'].isna()]
    if not missing_genes.empty:
        print(f"Warning: {missing_genes.gene.nunique()} genes in the annotation file are not present in the gene set.")
        print("These genes will be dropped from the final DataFrame.")
    
    return annot_df_merged

def combine_annotations(directory):
    """
    Combines all annotation files in a directory into a single DataFrame.
    
    Args:
        directory (str): Path to the directory containing annotation files.
        
    Returns:
        pd.DataFrame: Combined DataFrame with all annotations.
    """
    all_files = glob.glob(os.path.join(directory, '*.REGENIE.annotationFile.txt'))
    all_dfs = []

    for file in all_files:
        print(f"Reading annotation file: {file}")
        df = read_regenie_annotation(file)
        all_dfs.append(df)

    combined_df = pd.concat(all_dfs, ignore_index=True)
    return combined_df

def combine_setlists(directory):
    """
    Combines all setlist files in a directory into a single DataFrame.
    
    Args:
        directory (str): Path to the directory containing setlist files.
        
    Returns:
        pd.DataFrame: Combined DataFrame with all setlists.
    """
    all_files = glob.glob(os.path.join(directory, '*.REGENIE.setListFile.txt'))
    all_dfs = []

    for file in all_files:
        print(f"Reading setlist file: {file}")
        df = read_setlist(file)
        all_dfs.append(df)

    combined_df = pd.concat(all_dfs, ignore_index=True)
    return combined_df

def double_occurances(combined_df, col="snp"):
    duplicates = combined_df[combined_df.duplicated(subset=col, keep=False)]
    print(f"Found {duplicates[col].nunique()} {col}s that appear more than once across all files.")
    return duplicates


def main():
    annot_filepath = 'data/PTV_test/PTV_test.chr1.REGENIE.annotationFile.txt'
    geneset_df = read_gene_sets('data/burden_test_modules.obj')
    regenie_setlist = read_setlist('data/PTV_test/PTV_test.chr1.REGENIE.setListFile.txt')
    annot_df = read_regenie_annotation(annot_filepath)
    combined_annot_df = combine_annotations('data/PTV_test')
    double_occurances(combined_annot_df)
    combined_setlist_df = combine_setlists('data/PTV_test')
    double_occurances(combined_setlist_df, col='transcript')

    # transcript to gene mapping
    transcript_to_gene = pd.read_csv('data/transcript_gene_map.csv')
    transcript_to_gene.columns = ['chr','transcript', 'gene', 'gene_symbol']

    # compare transcript to gene mapping with the geneset_df
    gs = geneset_df.groupby("gene")
    # to df
    gs = gs.agg(list).reset_index()
    merged_transcript_df = transcript_to_gene.merge(gs, on='gene', how='inner')
    print(f"Found {merged_transcript_df.gene.nunique()} matching genes between transcript and gene set.")
    # print how many not found
    not_found = gs[~gs['gene'].isin(transcript_to_gene['gene'])]
    print(f"Found {not_found.gene.nunique()} genes in the gene set that are not in the transcript to gene mapping.")
    # not found the other way around
    not_found_transcript = transcript_to_gene[~transcript_to_gene['gene'].isin(gs['gene'])]
    print(f"Found {not_found_transcript.gene.nunique()} genes in the transcript to gene mapping that are not in the gene set.")

    # are there duplicates in the merged_transcript_df?
    duplicates = merged_transcript_df[merged_transcript_df.duplicated(subset=['transcript'], keep=False)]
    if not duplicates.empty:
        print(f"Warning: Found {duplicates.transcript.nunique()} transcripts that appear more than once in the merged transcript to gene mapping.")
    else:
        print("No duplicates found in the merged transcript to gene mapping.")


    # map the setlist to the transcript to gene mapping
    combined_setlist_df = combined_setlist_df.merge(merged_transcript_df, on='transcript', how='left')
    # report if there are any transcripts in the setlist that are not in the transcript to gene mapping
    missing_transcripts = combined_setlist_df[combined_setlist_df['gene'].isna()]
    if not missing_transcripts.empty:
        print(f"Warning: {missing_transcripts.transcript.nunique()} transcripts in the setlist are not present in the transcript to gene mapping.")
    else:
        print("All transcripts in the setlist are present in the transcript to gene mapping.")

    print(combined_setlist_df.head())

    # same for annot_df
    combined_annot_df = combined_annot_df.merge(merged_transcript_df, on='transcript', how='left')
    # report if there are any transcripts in the annotation that are not in the transcript to gene mapping
    missing_transcripts_annot = combined_annot_df[combined_annot_df['gene'].isna()]
    if not missing_transcripts_annot.empty:
        print(f"Warning: {missing_transcripts_annot.transcript.nunique()} transcripts in the annotation are not present in the transcript to gene mapping.")
    else:
        print("All transcripts in the annotation are present in the transcript to gene mapping.")

    # save combined_annot_df to a file
    #combined_annot_df.to_csv('data/PTV_all_snps.csv', index=False)
    # save combined_setlist_df to a file
    #combined_setlist_df.to_csv('data/PTV_all_transcripts.csv', index=False)
    #converted_df = convert_annotation(annot_filepath, geneset_df)
    #print(converted_df.head())



if __name__ == "__main__":
    main()