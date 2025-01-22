use pyo3::prelude::*;
use rand::seq::SliceRandom;
use rand::thread_rng;
use std::collections::HashMap;

pub type KmerCounter<'a> = HashMap<&'a str, usize>;

pub fn count_kmers(seq: &str, k: usize) -> KmerCounter {
    let mut cnt = KmerCounter::new();

    if seq.len() < k {
        return cnt;
    }

    for i in 0..(seq.len() - k + 1) {
        let slice = &seq[i..i + k];
        *cnt.entry(slice).or_insert(0) += 1;
    }

    cnt
}

fn shuffle_str(seq: &str) -> String {
    let mut seq_perm: Vec<char> = seq.chars().collect();
    let mut rng = thread_rng();
    seq_perm.shuffle(&mut rng);
    seq_perm.into_iter().collect()
}

fn shuffle_str_kmer_r_method(seq: &str, k: i32) -> String {
    // Gueto implementation
    let mut comb_k: Vec<&str> = vec![]; 
    if k==2 {
        comb_k = vec!["A", "T", "C", "G"]; 
    } else if k==3 {
        comb_k = vec!["AA","AT","AC","AG","TA","TT","TC","TG","CA","CT","CC","CG","GA","GT","GC","GG"]; 
    } else {
        !todo!()// Raise an error
    }
    // comb_k should be shuffle also
    let mut rng = thread_rng();
    comb_k.shuffle(&mut rng); 

    //permute for every (k-1)mer
    let mut seq_p: String = seq.to_string();
    for kmer in &comb_k {
        //split  
        let mut perm: Vec<String> = seq_p.split(kmer).map(|s| s.to_string()).collect();
        //shuffle 1 to n-1, ignore 
        if perm.len() > 2 {
            let mut rng = thread_rng();
            let size = perm.len(); 
            perm[1..size - 1].shuffle(&mut rng); // Shuffle only the middle elements
        }
        //joint
        seq_p = perm.join(kmer);
    }
    seq_p 

}

#[pyfunction]
pub fn count_kmers_py(seq: &str, k: usize) -> PyResult<KmerCounter> {
    Ok(count_kmers(seq, k))
}

#[pyfunction]
fn shuffle_str_py(seq: &str) -> PyResult<String> {
    Ok(shuffle_str(seq))
}
#[pyfunction]
fn shuffle_str_kmer_r_method_py(seq: &str, k: i32) -> PyResult<String> {
    Ok(shuffle_str_kmer_r_method(seq,k))
}

#[pymodule]
fn kmers(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(count_kmers_py, m)?)?;
    m.add_function(wrap_pyfunction!(shuffle_str_py, m)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_no_kmers() {
        let expected = KmerCounter::new();
        assert_eq!(count_kmers("AAA", 6), expected);
    }

    #[test]
    fn test_small_seq() {
        let mut expected = KmerCounter::new();
        expected.insert("AAACGT", 1);
        expected.insert("AACGTG", 1);
        assert_eq!(count_kmers("AAACGTG", 6), expected);
    }

    #[test]
    fn test_shuffle_str() {
        let seq = "AACGGTA";
        let seq_perm = shuffle_str(seq);
        assert_eq!(count_kmers(seq, 1), count_kmers(&seq_perm, 1))
    }
    #[test]
    fn test_shuffle_str_r_2() {
        let seq = "CAACTGGGCACATAATGCGTACGCCCATCTAGTACACCCA";
        let seq_perm = shuffle_str_kmer_r_method(seq,2);
        // Print the original and shuffled sequences
        println!("Original sequence: {}", seq);
        println!("Shuffled sequence: {}", seq_perm);
        // Check that `seq_perm` and `seq` are different
        assert_ne!(seq, &seq_perm);
        // Check that the k-mer counts are still the same
        assert_eq!(count_kmers(seq, 2), count_kmers(&seq_perm, 2))
    }
    #[test]
    fn test_shuffle_str_r_3() {
        let seq = "CAACTGGGCACATAATGCGTACGCCCATCTAGTACACCCA";
        let seq_perm = shuffle_str_kmer_r_method(seq,3);
        // Print the original and shuffled sequences
        println!("Original sequence: {}", seq);
        println!("Shuffled sequence: {}", seq_perm);
        // Check that `seq_perm` and `seq` are different
        assert_ne!(seq, &seq_perm);
        // Check that the k-mer counts are still the same
        assert_eq!(count_kmers(seq, 3), count_kmers(&seq_perm, 3))
    }
}
