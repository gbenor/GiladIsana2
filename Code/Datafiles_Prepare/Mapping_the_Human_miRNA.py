from Pipeline import *
from pathlib import Path
import JsonLog
import pandas as pd

def read_paper_data(f):
    return pd.read_csv(f, sep="\t", skiprows=30)


def df_prepare (in_df):
    in_df['GI_ID']= range(len(in_df))
    in_df.rename(columns={'miRNA_seq': 'miRNA sequence', 'mRNA_seq_extended': 'target sequence',
                          'chimeras_decompressed': 'number of reads'}, inplace=True)

    #  in_df.rename(columns={'miRNA ID': 'microRNA_name'}, inplace=True)
    return in_df



def main ():
    interaction_file = str(Path("Raw/1-s2.0-S009286741300439X-mmc1.txt"))
    log_dir = "Datafiles_Prepare/Logs/"

    organisms = ["Human"]
    for organism in organisms:
        p_dir = Path(organism) / "Raw"
        JsonLog.set_filename(
            filename_date_append(Path(log_dir) / Path("Mapping_the_Human_miRNA_" + organism + ".json")))
        JsonLog.add_to_json('file name', interaction_file)
        JsonLog.add_to_json('paper',
                            "Mapping the Human miRNA Interactome by CLASH Reveals Frequent Noncanonical Binding")
        JsonLog.add_to_json('Organism', organism)
        JsonLog.add_to_json('paper_url', "https://www.sciencedirect.com/science/article/pii/S009286741300439X")
        p = Pipeline(paper_name="Mapping_the_Human_miRNA",
                     organism=organism,
                     in_df=df_prepare(read_paper_data(interaction_file)),
                     data_dir=p_dir)
        p.run()
        p.file_formatting()

if __name__ == "__main__":
    main()
