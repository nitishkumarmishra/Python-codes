import pandas as pd

def get_largest_cds(gtf_file, output_file):
    # Specify the data types for each column
    dtype = {
        'seqname': str,
        'source': str,
        'feature': str,
        'start': int,
        'end': int,
        'score': str,
        'strand': str,
        'frame': str,
        'attribute': str
    }

    # Read the GTF file into a pandas DataFrame
    df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, dtype=dtype, low_memory=False,
                     names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

    # Filter for CDS features
    cds_df = df[df['feature'] == 'CDS']

    # Extract gene_id from the attribute column
    cds_df['gene_id'] = cds_df['attribute'].str.extract('gene_id "([^"]+)"')

    
    cds_df['gene_name'] = cds_df['attribute'].str.extract('gene_name "([^"]+)"')


    cds_df['transcript_id'] = cds_df['attribute'].str.extract('transcript_id "([^"]+)"')

    
    # Calculate the length of each CDS
    cds_df['length'] = cds_df['end'] - cds_df['start'] + 1


    # Rename seqname with Chr  
    cds_df['Chr'] = cds_df['seqname']

    # Find the largest CDS for each gene
    largest_cds_df = cds_df.loc[cds_df.groupby('gene_id')['length'].idxmax()]
    #largest_cds_df = cds_df.loc[cds_df.groupby('transcript_id')['length'].idxmax()]

    # Write the output to a text file
    #largest_cds_df.to_csv(output_file, sep='\t', index=False)

    # Select specific columns to write to the output file
    selected_columns = ['transcript_id', 'gene_name', 'gene_id', 'Chr', 'start', 'end', 'length']
    largest_cds_df[selected_columns].to_csv(output_file, sep='\t', index=False)


# Example usage
gtf_file = 'Mus_musculus.GRCm39.104.rdna_rn18s.gtf'
output_file = 'python_largest_cds.txt'
get_largest_cds(gtf_file, output_file)
print(f"Largest CDS information has been written to {output_file}")


