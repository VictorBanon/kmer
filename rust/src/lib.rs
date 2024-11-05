use pyo3::prelude::*;
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

#[pyfunction]
pub fn count_kmers_py(seq: &str, k: usize) -> PyResult<KmerCounter> {
    Ok(count_kmers(seq, k))
}

#[pymodule]
fn kmers(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(count_kmers_py, m)?)
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
}
