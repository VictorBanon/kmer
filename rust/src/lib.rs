use std::collections::HashMap;

pub type KmerCounter = Vec<(String, usize)>;

pub fn count_kmers(seq: &str, k: usize) -> KmerCounter {
    let mut cnt = HashMap::new();

    if seq.len() < k {
        return Vec::new();
    }

    for i in 0..(seq.len() - k + 1) {
        let slice = &seq[i..i + k];
        *cnt.entry(slice).or_insert(0) += 1;
    }

    let mut sorted_counts: KmerCounter = cnt.iter().map(|(k, v)| (k.to_string(), *v)).collect();
    sorted_counts.sort_by(|a, b| {
        b.1.cmp(&a.1) // Counts in descending order
            .then_with(|| a.0.cmp(&b.0)) // Then keys in ascending order
    });

    sorted_counts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_no_kmers() {
        let expected: KmerCounter = Vec::new();
        assert_eq!(count_kmers("AAA", 6), expected);
    }

    #[test]
    fn test_small_seq() {
        let mut expected = KmerCounter::new();
        expected.push(("AAACGT".to_string(), 1));
        expected.push(("AACGTG".to_string(), 1));
        assert_eq!(count_kmers("AAACGTG", 6), expected);
    }
}
