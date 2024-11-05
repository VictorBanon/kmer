import copy
import gzip
import random
import time
from collections import Counter
from pathlib import Path

# import polars as pl
import numpy as np
import pandas as pd
from Bio import SeqIO

# from scripts.count_kmer import count_kmers
from kmers import count_kmers_py as count_kmers

SPATH = Path(__file__).parent.parent


def _main(
    ipath=SPATH / "data",
    opath=SPATH / "result",
    number_replica=10,
    use_cache=False,
):
    assert number_replica > 0

    folders = [f.name for f in ipath.iterdir() if f.is_dir()]

    tax_path = ipath / "list_with_taxonomy.csv"

    # Import CSV file into DataFrame
    tax_df = pd.read_csv(tax_path)[["#Organism Name", "Replicons"]]  # family # #Organism Name

    # Step 1: Split the column 'terms' by ';' to create a list in each row
    tax_df["Replicons"] = tax_df["Replicons"].str.split(";")

    # Step 2: Use explode() to create a new row for each term
    tax_df = tax_df.explode("Replicons")

    for id, folder in enumerate(folders[:]):
        percent = 100 * id / len(folders)
        print(f"Completion percent: {percent}%")
        genomic = ipath / folder / f"{folder}_genomic.fna.gz"

        with gzip.open(genomic, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                file_output = opath / f"{record.id}.csv"

                if use_cache and file_output.exists():
                    print(f"{file_output} already exists. Skipping {genomic}.")
                else:
                    print(f"Computing {genomic}")

                    six_mer_obs = Counter()
                    six_mer_obs += count_kmers(str(record.seq), k=6)
                    six_mer_sim = copy.deepcopy(six_mer_obs)
                    for _ in range(number_replica - 1):
                        # Unir la lista de nuevo en una cadena
                        seq_perm = list(str(record.seq))
                        random.shuffle(seq_perm)
                        seq_perm = "".join(seq_perm)
                        six_mer_sim += count_kmers(seq_perm, k=6)

                    # Divide counter1 by counter2
                    six_mer_obs = Counter(
                        {
                            key: np.log(six_mer_obs[key] / (1 + six_mer_sim[key] / 10))
                            for key in six_mer_obs
                        }
                    )

                    # Write file
                    df = pd.DataFrame(six_mer_obs.items(), columns=["Item", "Count obs"])
                    df.to_csv(file_output, index=False)
                    print(f"Finished {genomic}")
                    print(f"The output can be found at {file_output}")

                break  # end record
        break  # end folders


def main():
    start_time = time.time()
    _main()
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Execution time: {execution_time:.2f} seconds")


if __name__ == "__main__":
    main()
