use csv::Writer;
use kmers::{count_kmers, KmerCounter};
use rand::Rng;
use std::time::Instant;

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
fn make_seq(size: usize) -> String {
    let choices = ['A', 'C', 'G', 'T'];
    let mut rng = rand::thread_rng();
    (0..size).map(|_| choices[rng.gen_range(0..4)]).collect()
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

        // Write to CSV file for some exponent to test the writer
        if exp == 6 {
            write_to_csv("kmer_count.csv", &cnt);
        }
    }
}
