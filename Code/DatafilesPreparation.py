import pandas as pd
import numpy as np
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
from Bio.Blast import NCBIXML



blast_dir = "\"C:\\Program Files\\NCBI\\blast-2.7.1+\\bin\""



def ENSG (name, sep) :
    a=name.split(sep)
    return a[0]

def ENST (name, sep) :
    a=name.split(sep)
    return a[1]


def generate_fasta (fasta_filename, mrna_df):
    # Generate fasta file with all mRNA
    sequences = []
    for index, row in mrna_df.iterrows():
        mRNA_seq = row.mRNA_seq_extended
        mRNA_name = row.mRNA_name
        record = SeqRecord(Seq(mRNA_seq), description=mRNA_name)
        sequences.append(record)
    SeqIO.write(sequences, fasta_filename, "fasta")


def make_blast_db (db_title, mrna_df):
    sequences = []
    for index, row in mrna_df.iterrows():
        mRNA_seq = row.sequence
        mRNA_name = row.description
        record = SeqRecord(Seq(mRNA_seq), id=str(index),   description=mRNA_name)
        sequences.append(record)
    fasta_filename=db_title + ".fasta"
    SeqIO.write(sequences, fasta_filename, "fasta")

    cmd = "makeblastdb -in {fasta} -parse_seqids -dbtype nucl -out {out}".format(dir=blast_dir, fasta=fasta_filename, out=db_title)
    os.system(cmd)
    os.remove(fasta_filename)


def run_blastn (mrna, db_title, blast_output_filname):
    mRNA_seq = mrna
    mRNA_name = "mrna_to_find"
    filename = "mrna_to_find.fasta"
    record = SeqRecord(Seq(mRNA_seq), description=mRNA_name)
    SeqIO.write(record, filename, "fasta")

    cline = NcbiblastnCommandline(query=filename, db=db_title, strand="plus", evalue=0.001,
                                  out=blast_output_filname, outfmt=5, max_hsps=1, max_target_seqs=1) # it should return the best result
    cmd = str(cline)
    print (cmd)
    os.system(cmd)

def parse_blast (blast_result, mrna_seq_candidate):
    #assuming only the best result stored in blast_result
    E_VALUE_THRESH = 10
    blast_xml_handle = open(blast_result, "r")
    records = NCBIXML.parse(blast_xml_handle)
    for b in records:
        if b.alignments:  # skip queries with no matches
            for align in b.alignments:
                for hsp in align.hsps:
                    if hsp.expect < E_VALUE_THRESH:
                        identities = hsp.identities
                        enst = ENST(align.hit_def,"|")
    try:
        full_mrna = mrna_seq_candidate[mrna_seq_candidate.enst == ENST(align.hit_def, "|")].sequence.item()
        return enst, identities, full_mrna
    except UnboundLocalError:
        return np.nan, np.nan, np.nan



def CLASH_3UTR_filter (clash_file, output_file):
    """
    This function read Human CLASH file, filters out the non 3'utr and generate a CVS file.
    Args:
        clash_file:
        output_file:

    Returns:

    """
    # Read clash data and print its summary information
    human_clash_data = pd.read_csv(clash_file, sep="\t", skiprows=30)
    print ("The CLASH file contains {} rows".format(human_clash_data.shape[0]))

    # filtter out all the rows without 3'UTR
    human_clash_data_utr3 = human_clash_data[human_clash_data['3\'UTR'].notnull()]
    human_clash_data_utr3['ensg'] = human_clash_data_utr3.apply(lambda row: ENSG(row['mRNA_name'], "_"), axis=1)
    human_clash_data_utr3['enst'] = human_clash_data_utr3.apply(lambda row: ENST(row['mRNA_name'], "_"), axis=1)
    print ("After filtering all the non 3utr, the CLASH file contains {} rows".format(human_clash_data_utr3.shape[0]))

    human_clash_data_utr3.to_csv(output_file)



def CLASH_BIOMART_MERGE_FULL_MATCH (clash_file, biomart_file, output_file):
    print ("Start merging CLASH<-BIOMART")
    clash_data = pd.read_csv(clash_file)
    biomart_df = pd.read_csv(biomart_file)

    clash_data['full_mrna_seq'] = np.nan
    clash_data['full_mrna_source'] = np.nan

    for clash_index, clash_row in clash_data.iterrows():
        #find candidates with the same ensg
        a=biomart_df['ensg']==clash_row.ensg
        if sum(a)>0:
            full_mrna_seq_candidate=biomart_df[a]
            # Choose the candidates that contains the mRNA_seq_extended
            f_result=[]
            for candidate_index, candidate_row in full_mrna_seq_candidate.iterrows():
                try:
                    f = candidate_row.sequence.find(clash_row.mRNA_seq_extended)
                except Exception:
                    f = -1
                f_result.append(f)
                #Insert the candidates to a dedicated dataframe
            f_result=np.array(f_result)
            b=f_result!=-1
            if sum(b)>0:
                final_candidates=full_mrna_seq_candidate[b]
                #find the candidate with the longest full mrna seq
                sequence_len = final_candidates.sequence.str.len()
                row_with_the_longest_sequence = final_candidates[sequence_len==sequence_len.max()].iloc[0]
                clash_data.loc[clash_index,'full_mrna_seq'] = row_with_the_longest_sequence.sequence
                clash_data.loc[clash_index,'full_mrna_source'] = "BIOMART"


    clash_data.to_csv(output_file)
    print ("Finish merging. There are {} rows without full mrna seq".format(clash_data['full_mrna_seq'].isna().sum()))



def CLASH_BIOMART_MERGE_BLASTN (clash_file, biomart_file, output_file):
    """
    I run this function on my windows machine. I has some problems with libnsl.so.1 on my virtual machine
    This function handle the "problematic" mrna seq:
    seq that we havn't found the full mrna seq during the first step.
    the function runs blast against each "problemtic" row.
    Args:
        clash_file:
        biomart_file:
        output_file:

    Returns:

    """
    row_count = 0
    clash_data = pd.read_csv(clash_file)
    biomart_df = pd.read_csv(biomart_file)
    clash_data['full_mrna_identities'] = np.nan
    clash_data['full_mrna_enst'] = np.nan

    rows = clash_data[clash_data['full_mrna_seq'].isna()] # rows with no match

    for clash_index, clash_row in rows.iterrows():
        #find candidates with the same ensg
        a=biomart_df['ensg']==clash_row.ensg
        if sum(a)==0:
            clash_data.loc[clash_index, 'full_mrna_seq'] = "No ENSG match"
            continue

        if sum(a)>0:
            full_mrna_seq_candidate=biomart_df[a]
            db_title = "biomart" + str(clash_row.ensg)
            make_blast_db (db_title, full_mrna_seq_candidate)
            blast_output_filname = "blastout.xml"
            mrna = clash_row.mRNA_seq_extended
            run_blastn(mrna, db_title, blast_output_filname)
            enst, identities, full_mrna = parse_blast(blast_output_filname, full_mrna_seq_candidate)
            clash_data.loc[clash_index,'full_mrna_seq'] = full_mrna
            clash_data.loc[clash_index,'full_mrna_source'] = "BLASTN"
            clash_data.loc[clash_index,'full_mrna_identities'] = identities
            clash_data.loc[clash_index,'full_mrna_enst'] = enst
            # os.remove("biomart.*")
            # os.remove(blast_output_filname)
            # os.remove("mrna_to_find.fasta")
            row_count+=1

    clash_data.to_csv(output_file)


    print ("Finish merging. Merging {} rows.".format(row_count))



def main():
    os.sys.path.append(blast_dir)
    # CLASH_3UTR_filter ("Data/Human/Raw/1-s2.0-S009286741300439X-mmc1.txt",
    #                    "Data/Human/Parsed/human_clash_data_utr3.csv")
    #
    CLASH_BIOMART_MERGE_FULL_MATCH ("Data/Human/Parsed/human_clash_data_utr3.csv",
                         "Data/Human/Raw/biomart_3utr.csv",
                         "Data/Human/Parsed/human_clash_data_utr3_with_biomart_seq.csv")

    CLASH_BIOMART_MERGE_BLASTN ("Data/Human/Parsed/human_clash_data_utr3_with_biomart_seq.csv",
                         "Data/Human/Raw/biomart_3utr.csv",
                         "Data/Human/Parsed/human_clash_data_utr3_with_biomart_seq_step2.csv")





if __name__ == '__main__':
    main()

#

