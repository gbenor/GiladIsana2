import BlastUtils
from pathlib import Path
import JsonLog
from Bio import SeqIO
import datetime
import pandas as pd
from datetime import datetime
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



def filename_suffix_append(f, s):
    f = Path(f)
    return f.parent / Path(f.stem + s + f.suffix)

def filename_date_append (f):
    return filename_suffix_append (f, "_{}".format(datetime.now().strftime("%Y%m%d-%H%M%S")))



class Pipeline(object):

    def __init__(self, paper_name, organism, in_df, data_dir):
        self.necessary_columns = set(['GI_ID', 'microRNA_name', 'miRNA sequence', 'target sequence', 'number of reads'])

        self.paper_name = paper_name
        self.organism = organism
        self.data_dir = data_dir

        self.Pipeline_Files()
        self.inter_df = in_df
        assert self.necessary_columns.issubset(self.inter_df.columns), \
            "The input file doesn't contain all the necessary columns for the pipeline. \n " \
            "necessary columns: {} \n" \
            "file columns: {}".format(self.necessary_columns, self.inter_df.columns)


        JsonLog.add_to_json('Pipeline input samples count', self.inter_df.shape[0])
        print ("##################################################")
        print ("Starting the data prepare pipelne:")
        print ("Papar: \t{}".format(self.paper_name))
        print ("Organism:\t{}".format(self.organism))
        print ("Samples:\t{}".format(self.inter_df.shape[0]))
        print ("##################################################")


    def Pipeline_Files (self):
        self.blast_db_fasta_biomart = str(self.data_dir / Path (self.organism +"_biomart.fasta"))
        self.blast_db_fasta_biomart_valid = str(self.data_dir / Path (self.organism +"_biomart_valid.fasta"))
        self.blast_db_biomart = str(self.data_dir / "blast_files"  / "blastdb_biomart")

        self.blast_db_csv_ucsc = str(self.data_dir / Path(self.organism + "_ucsc.csv"))
        self.blast_db_fasta_ucsc_valid = str(self.data_dir / Path(self.organism + "_ucsc_valid.fasta"))
        self.blast_db_ucsc = str(self.data_dir / "blast_files"  / "blastdb_ucsc")
        self.blast_tmp = str(self.data_dir / "blast_files" /"tmp" / "tmp_result.xml")
        self.blast_result = str(self.data_dir / Path (self.organism +"_" +self.paper_name +"_blast.csv"))
        self.blast_result = filename_date_append(self.blast_result)
        self.blast_with_mRNA = str(self.data_dir / Path (self.organism +"_mRNA.csv"))
        self.final_output = str("Datafiles_Prepare/CSV"/ Path (self.organism +"_" +self.paper_name +"_Data.csv"))
        self.final_output = filename_date_append(self.final_output)

    def create_blast_db_biomart(self):
        fasta_sequences = SeqIO.parse(open(self.blast_db_fasta_biomart), 'fasta')
        valid_seq = []
        c = 0
        for s in fasta_sequences:
            c+=1
            if s.seq != 'Sequenceunavailable':
                valid_seq.append(s)
        SeqIO.write(valid_seq, self.blast_db_fasta_biomart_valid, 'fasta')
        print ("Finish filtering Biomart fasta")
        JsonLog.add_to_json("Total Biomart seq", c)
        JsonLog.add_to_json("Valid Biomart seq", len (valid_seq))

        BlastUtils.create_blast_db(self.blast_db_fasta_biomart_valid, self.blast_db_biomart)

    def create_blast_db_ucsc (self, seq_min_len=10):
        in_df = pd.read_csv(self.blast_db_csv_ucsc)

        valid_seq = []
        c = 0
        for index, row in in_df.iterrows():
            c+=1
            if len(row.sequence)>=seq_min_len:
                record = SeqRecord(Seq(row.sequence),
                                   id="{}_{}_{}".format(row.group_name, row.UCSCKG_ID, index))
                valid_seq.append(record)
        SeqIO.write(valid_seq, self.blast_db_fasta_ucsc_valid, 'fasta')
        print ("Finish filtering UCSC fasta")
        JsonLog.add_to_json("Total UCSC seq", c)
        JsonLog.add_to_json("Valid UCSC seq", len(valid_seq))

        BlastUtils.create_blast_db(self.blast_db_fasta_ucsc_valid, self.blast_db_ucsc)


    def blast_mRNA(self, db_title="biomart"):
        blast_db_title = self.blast_db_biomart if db_title=="biomart" else self.blast_db_ucsc
        output_file = self.blast_result
        blast_tmp = self.blast_tmp
        inter_df = self.inter_df

        for index, row in self.inter_df.iterrows():
            seq_to_blast = row['target sequence']
            BlastUtils.run_blastn(seq_to_blast, blast_db_title, blast_tmp)
            full_match, title, sbjct_start, sbjct_end, identities = BlastUtils.parse_blast(blast_tmp, seq_to_blast)
            inter_df.loc[index, db_title + '_full_match'] = full_match
            inter_df.loc[index, db_title + '_title'] = title
            inter_df.loc[index, db_title + '_sbjct_start'] = sbjct_start
            inter_df.loc[index, db_title + '_sbjct_end'] = sbjct_end
            inter_df.loc[index, db_title + '_identities'] = identities
            if index > 100:
                break

        inter_df.to_csv(output_file)

    def add_full_mRNA (self, db_title="biomart"):
        blast_result = self.blast_result
        blast_db_fasta = self.blast_db_fasta_biomart
        output_file = self.blast_with_mRNA

        record_dict = SeqIO.index(blast_db_fasta, "fasta")
        inter_df = pd.read_csv(blast_result)
        for index, row in inter_df.iterrows():
            try:
                key = row[db_title + '_title'].split()[0]
                full_mrna = str(record_dict[key].seq)
                inter_df.loc[index, db_title + '_full_mrna'] = full_mrna
                inter_df.loc[index, db_title + '_blast_status'] = "OK" if full_mrna.find(row['target sequence']) != -1 else "Not OK"
            except AttributeError:
                inter_df.loc[index, db_title + '_blast_status'] = "Not OK"

        inter_df.to_csv(output_file)



    def file_formatting (self):
        in_df = pd.read_csv(self.blast_with_mRNA)
        in_df['Source'] = self.paper_name
        in_df['Organism'] = self.organism
        in_df['Creation_time'] = JsonLog.get_creation_time()

        # Take only the valid rows: Status=OK
        valid_rows = in_df['biomart_blast_status']=="OK"
        JsonLog.add_to_json("Pipeline valid blast results", sum(valid_rows))
        JsonLog.add_to_json("Pipeline invalid blast results",(in_df.shape[0] - sum(valid_rows)))
        in_df = in_df[valid_rows]

        #Remove miRNA with XXX/stars
        rows_without_XXX = in_df['miRNA sequence'].apply(lambda x: x.find('X')==-1)
        JsonLog.add_to_json("Pipeline valid miRNA_no_xxx" ,sum(rows_without_XXX))
        JsonLog.add_to_json("Pipeline invalid miRNA (xxx)",(in_df.shape[0] - sum(rows_without_XXX)))
        in_df = in_df[rows_without_XXX]

        # Remove miRNA with stars
        rows_without_stars = in_df['microRNA_name'].apply(lambda x: x.find('star') == -1)
        JsonLog.add_to_json("Pipeline valid miRNA_no_***", sum(rows_without_stars))
        JsonLog.add_to_json("Pipeline invalid miRNA (star)", (in_df.shape[0] - sum(rows_without_stars)))
        in_df = in_df[rows_without_stars]

        # Choose the necessary columns
        in_df_filter = in_df.filter(
            ['Source', 'Organism', 'GI_ID', 'microRNA_name', 'miRNA sequence', 'target sequence', 'number of reads',
             'biomart_title', 'biomart_sbjct_start', 'biomart_sbjct_end', 'full_mrna'], axis=1)

        in_df_filter.rename(columns={'biomart_title': 'mRNA_name', 'biomart_sbjct_start': 'mRNA_start',
                              'biomart_sbjct_end': 'mRNA_end'}, inplace=True)


        # reset the index
        in_df_filter.reset_index(drop=True, inplace=True)
        # self.log.append("miRNA statistics")
        # self.log.append(dict(in_df_filter['microRNA_name'].value_counts()))

        # save to file
        in_df_filter.to_csv(filename_date_append(self.final_output))


    def run (self):
        #####################################################
        # build the blast DB
        #####################################################
        self.create_blast_db_biomart()
        # if self.organism == "human":
        #     self.create_blast_db_ucsc()


        #####################################################
        # Run the blast against biomart
        #####################################################
        self.blast_mRNA("biomart")
        #self.blast_mRNA("ucsc")


        #####################################################
        # Add the full mRNA seq and report the status
        #####################################################
        self.add_full_mRNA()

        #####################################################
        # Remove invalid rows and format the data
        #####################################################
        self.file_formatting()




