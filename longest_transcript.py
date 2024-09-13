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


    # Filter for transcript features
    transcript_df = df[df['feature'] == 'transcript']

    # Extract gene_id from the attribute column
    transcript_df['gene_id'] = transcript_df['attribute'].str.extract('gene_id "([^"]+)"')


    transcript_df['gene_name'] = transcript_df['attribute'].str.extract('gene_name "([^"]+)"')


    transcript_df['transcript_id'] = transcript_df['attribute'].str.extract('transcript_id "([^"]+)"')


    # Calculate the length of each transcript
    transcript_df['length'] = transcript_df['end'] - transcript_df['start'] + 1


    # Rename seqname with Chr
    transcript_df['Chr'] = transcript_df['seqname']

    # Find the largest transcript for each gene
    largest_transcript_df = transcript_df.loc[transcript_df.groupby('gene_id')['length'].idxmax()]
    #largest_cds_df = cds_df.loc[cds_df.groupby('transcript_id')['length'].idxmax()]

    # Write the output to a text file
    #largest_cds_df.to_csv(output_file, sep='\t', index=False)

    # Select specific columns to write to the output file
    selected_columns = ['transcript_id', 'gene_name', 'gene_id', 'Chr', 'start', 'end', 'length']
    largest_transcript_df[selected_columns].to_csv(output_file, sep='\t', index=False)



# Example usage
gtf_file = 'Mus_musculus.GRCm39.104.rdna_rn18s.gtf'
output_file = 'python_largest_transcript.txt'
get_largest_cds(gtf_file, output_file)
print(f"Largest CDS information has been written to {output_file}")


