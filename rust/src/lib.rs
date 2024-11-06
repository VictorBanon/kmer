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

#[pyfunction]
pub fn count_kmers_py(seq: &str, k: usize) -> PyResult<KmerCounter> {
    Ok(count_kmers(seq, k))
}

#[pyfunction]
fn shuffle_str_py(seq: &str) -> PyResult<String> {
    Ok(shuffle_str(seq))
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
}
