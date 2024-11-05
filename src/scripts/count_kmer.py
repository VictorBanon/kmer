import time  
from collections import Counter
from Bio import SeqIO
import gzip  
import pandas as pd  
import numpy as np
import itertools
import polars as pl
import os
from itertools import product
import random
import copy
# import plotly.express as px 
# import plotly.io as pio

def count_kmers(seq, k=6):
    cnt = Counter()
    if len(seq) < k:
        return cnt

    kmer = seq[:k] 
    cnt[kmer] += 1
    for i in range(k, len(seq)):
        kmer = kmer[1:] + seq[i]
        cnt[kmer] += 1

    return cnt 

def main():
    path = "src/data" 
    path_output = "src/result" 

    folders = [f.name for f in os.scandir(path) if f.is_dir()] 

    # number_replica
    number_replica = 10

    # Import list with taxonomy
    tax_path = os.path.join(path,"list_with_taxonomy.csv")

    # Import CSV file into DataFrame
    tax_df = pd.read_csv(tax_path)[["#Organism Name","Replicons"]] # family # #Organism Name 

    # Step 1: Split the column 'terms' by ';' to create a list in each row
    tax_df['Replicons'] = tax_df['Replicons'].str.split(';')

    # Step 2: Use explode() to create a new row for each term
    tax_df = tax_df.explode('Replicons')

    S=len(folders)

    for id,folder in enumerate(folders[:]): 
            
        percent = 100 * id / S
        print(f"{percent}%")
        file = os.path.join(path,folder,folder+"_genomic.fna.gz")
        with gzip.open(file, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                # Check if the file exists
                file_output =os.path.join(path_output,record.id+".csv")
                if not os.path.exists(file_output):
                    print(record.id) 
                    for i in range(number_replica):
                        if i == 0:
                            six_mer_obs = count_kmers(str(record.seq))
                            six_mer_sim = copy.deepcopy(six_mer_obs)

                            # Convert the Counter to a DataFrame
                        else:  
                            # Unir la lista de nuevo en una cadena
                            seq_perm = list(str(record.seq)) 
                            random.shuffle(seq_perm) 
                            seq_perm = ''.join(seq_perm)
                            six_mer_sim += count_kmers(seq_perm) 
                    six_mer_sim = Counter({key: value / 10 for key, value in six_mer_sim.items()}) 
                    # Divide counter1 by counter2 
                    six_mer_obs = Counter({key: np.log(six_mer_obs[key] / (six_mer_sim[key]+1)) for key in six_mer_obs})  
                    # Save the DataFrame to the CSV file
                    df = pd.DataFrame(six_mer_obs.items(), columns=['Item', 'Count obs'])
                    df.to_csv(file_output, index=False)
                    print("done")        
                else:
                    print(f'{file_output} already exists.')
                break
        break
 
def test_no_kmers():
    expected = Counter()
    assert count_kmers("AAA") == expected

def test_small_seq():
    expected = Counter({"AAACGT": 1, "AACGTG": 1})
    assert count_kmers("AAACGTG") == expected

if __name__ == "__main__":
    # tests
    test_no_kmers()
    test_small_seq()
    # Extract 6mer 
    start_time = time.time()  # Record the start time
    main()              # Call the function
    end_time = time.time()    # Record the end time

    execution_time = end_time - start_time  # Calculate execution time
    print(f"Execution time: {execution_time:.2f} seconds") 