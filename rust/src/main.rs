use csv::Writer;
use rand::Rng;
use std::collections::HashMap;
use std::time::Instant;

type KmerCounter = Vec<(String, usize)>;

fn count_kmers(seq: &str, k: usize) -> KmerCounter {
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

fn make_seq(size: usize) -> String {
    let choices = ['A', 'C', 'G', 'T'];
    let mut rng = rand::thread_rng();
    (0..size).map(|_| choices[rng.gen_range(0..4)]).collect()
}

fn write_to_csv(file_path: &str, data: &KmerCounter) {
    let mut wtr = Writer::from_path(file_path).expect("Unable to create CSV file");

    wtr.write_record(&["kmer", "count"])
        .expect("Unable to write header");

    for (kmer, count) in data {
        wtr.write_record(&[kmer, &count.to_string()])
            .expect("Unable to write record");
    }

    wtr.flush().expect("Unable to flush CSV writer");
}

fn main() {
    for exp in 1..=8 {
        let seq_size = 10_usize.pow(exp);

        // Measure time to generate the sequence
        let seq_start = Instant::now();
        let seq = make_seq(seq_size);
        let seq_elapsed = seq_start.elapsed();

        // Measure time to count kmers
        let count_start = Instant::now();
        let cnt = count_kmers(&seq, 6);
        let count_elapsed = count_start.elapsed();

        println!(
            "exp = {}, sequence generation = {} sec, k-mer counting = {} sec",
            exp,
            seq_elapsed.as_secs_f32(),
            count_elapsed.as_secs_f32()
        );

        // Write to CSV file for some exponent to test the writter
        if exp == 6 {
            write_to_csv("kmer_count.csv", &cnt);
        }
    }
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
