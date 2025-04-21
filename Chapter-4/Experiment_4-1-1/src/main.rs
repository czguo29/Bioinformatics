use bio::io::fastq;
use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::collections::HashMap;

/// Fungsi untuk membangun Position Weight Matrix (PWM) dari urutan FASTQ.
fn build_pwm(sequences: &Vec<String>) -> Vec<HashMap<char, f64>> {
    let seq_length = sequences[0].len();  // Mendapatkan panjang urutan pertama
    let mut counts_per_position = vec![HashMap::new(); seq_length];  // Inisialisasi penghitung untuk setiap posisi

    // Menghitung jumlah nukleotida di setiap posisi
    for seq in sequences {
        for (i, base) in seq.chars().enumerate() {
            // Pastikan tidak mengakses indeks yang tidak valid
            if i < counts_per_position.len() {
                let count = counts_per_position[i].entry(base).or_insert(0.0);
                *count += 1.0;
            }
        }
    }

    // Mengubah jumlah menjadi probabilitas (normalisasi berdasarkan posisi)
    for position in 0..seq_length {
        let total_count: f64 = counts_per_position[position].values().sum();
        for base in &['A', 'C', 'G', 'T'] {
            let count = counts_per_position[position].get(base).unwrap_or(&0.0);
            let probability = if total_count > 0.0 {
                *count / total_count
            } else {
                0.0
            };
            counts_per_position[position].insert(*base, probability);
        }
    }

    counts_per_position
}

/// Fungsi untuk membangun Markov Random Field (MRF) dari urutan FASTQ.
fn build_mrf(sequences: &Vec<String>) -> HashMap<(char, char), f64> {
    let mut transition_counts = HashMap::new();
    let mut total_pairs = 0.0;

    // Inisialisasi penghitung transisi untuk pasangan (X->Y)
    for &base1 in &['A', 'C', 'G', 'T'] {
        for &base2 in &['A', 'C', 'G', 'T'] {
            transition_counts.insert((base1, base2), 0.0);
        }
    }

    // Menghitung transisi antar basa yang berurutan
    for seq in sequences {
        let chars: Vec<char> = seq.chars().collect();
        for i in 0..chars.len() - 1 {
            if let Some(count) = transition_counts.get_mut(&(chars[i], chars[i + 1])) {
                *count += 1.0;
                total_pairs += 1.0;
            }
        }
    }

    // Mengubah jumlah transisi menjadi probabilitas
    let mut transition_probabilities = HashMap::new();
    for (&pair, &count) in &transition_counts {
        let prob = if total_pairs > 0.0 {
            count / total_pairs
        } else {
            0.0
        };
        transition_probabilities.insert(pair, prob);
    }

    transition_probabilities
}

fn main() {
    // Menerima argumen command line:
    // 1) input FASTQ
    // 2) file output PWM
    // 3) file output MRF
    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        eprintln!("Usage: {} <input_fastq> <pwm_output> <mrf_output>", args[0]);
        std::process::exit(1);
    }

    let input_fastq = &args[1];
    let pwm_output_path = &args[2];
    let mrf_output_path = &args[3];

    // Membaca urutan dari file FASTQ
    let mut seqs = Vec::new();
    let reader = fastq::Reader::from_file(input_fastq).expect("Could not open FASTQ file");
    for record in reader.records() {
        let rec = record.expect("Error reading record");
        seqs.push(rec.seq().to_vec());
    }

    // Mengonversi urutan byte menjadi String
    let string_seqs: Vec<String> = seqs.iter()
                                       .map(|s| String::from_utf8_lossy(s).into_owned())
                                       .collect();

    // Membangun PWM
    let pwm = build_pwm(&string_seqs);

    // Membangun MRF (transisi orde pertama)
    let mrf = build_mrf(&string_seqs);

    // Menulis hasil PWM ke file output
    let pwm_file = File::create(pwm_output_path).expect("Cannot create PWM output file");
    let mut pwm_writer = BufWriter::new(pwm_file);

    writeln!(pwm_writer, "Position Weight Matrix (probabilities)").unwrap();
    for (pos, position_map) in pwm.iter().enumerate() {
        println!("Position {}: A={}, C={}, G={}, T={}",
            pos,
            position_map.get(&'A').unwrap_or(&0.0),
            position_map.get(&'C').unwrap_or(&0.0),
            position_map.get(&'G').unwrap_or(&0.0),
            position_map.get(&'T').unwrap_or(&0.0)
        ); // Debugging: Memeriksa nilai PWM
        writeln!(pwm_writer,
            "Position {}: A={:.3}, C={:.3}, G={:.3}, T={:.3}",
            pos,
            position_map[&'A'],
            position_map[&'C'],
            position_map[&'G'],
            position_map[&'T']
        ).unwrap();
    }

    // Menulis hasil MRF ke file output
    let mrf_file = File::create(mrf_output_path).expect("Cannot create MRF output file");
    let mut mrf_writer = BufWriter::new(mrf_file);

    writeln!(mrf_writer, "1st-order Markov Random Field (transition probabilities)").unwrap();
    for base1 in &['A', 'C', 'G', 'T'] {
        for base2 in &['A', 'C', 'G', 'T'] {
            let probability = mrf[&(*base1, *base2)];
            println!("Transition {}->{}: {:.4}", base1, base2, probability); // Debugging: Memeriksa transisi MRF
            writeln!(mrf_writer, "{}->{}: {:.4}", base1, base2, probability).unwrap();
        }
    }

    println!("PWM dan MRF telah berhasil dihitung dan ditulis ke file output!");
}
