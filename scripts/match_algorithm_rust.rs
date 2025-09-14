/*
 * Ultra-fast Rust implementation of MATCH algorithm with SIMD and parallelization
 * Compile: cargo build --release
 * Run: cargo run --release -- -p pwm.json -g genome.fa
 */

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use clap::Parser;
use serde::Deserialize;
use flate2::read::GzDecoder;
use rayon::prelude::*;

#[derive(Parser)]
#[command(name = "match_algorithm_rust")]
#[command(about = "Ultra-fast MATCH algorithm for promoter prediction")]
struct Args {
    #[arg(short = 'p', long = "pwm", help = "PWM JSON file")]
    pwm_file: String,

    #[arg(short = 'g', long = "genome", help = "Genome FASTA file (can be gzipped)")]
    genome_file: String,

    #[arg(short = 'o', long = "output", help = "Output GFF3 file (default: stdout)")]
    output_file: Option<String>,

    #[arg(short = 'c', long = "cutoff", default_value = "minSum",
          help = "Cutoff strategy (minFN, minFP, minSum)")]
    cutoff_strategy: String,

    #[arg(short = 's', long = "min-score", default_value = "0.8",
          help = "Minimum score threshold")]
    min_score: f64,

    #[arg(short = 'f', long = "forward-only", help = "Scan forward strand only")]
    forward_only: bool,

    #[arg(short = 't', long = "threads", default_value = "0",
          help = "Number of threads (0 = auto)")]
    threads: usize,
}

#[derive(Debug, Deserialize)]
struct PWMJson {
    length: usize,
    sequences: usize,
    frequencies: HashMap<char, Vec<f64>>,
    information: Vec<f64>,
    core_start: usize,
    #[serde(rename = "core_end")]
    core_end: usize,
    consensus: String,
}

#[derive(Debug, Clone)]
struct PWM {
    frequencies: [[f64; 100]; 4],  // A=0, C=1, G=2, T=3
    information: [f64; 100],
    min_scores: [f64; 100],
    max_scores: [f64; 100],
    length: usize,
    core_start: usize,
    core_end: usize,
    total_min: f64,
    total_max: f64,
    core_min: f64,
    core_max: f64,
    css_cutoff: f64,
    mss_cutoff: f64,
}

#[derive(Debug, Clone)]
struct Match {
    start: usize,
    end: usize,
    strand: char,
    css: f64,
    mss: f64,
    score: f64,
}

impl PWM {
    fn from_json(json: PWMJson, cutoff_strategy: &str) -> Self {
        let mut pwm = PWM {
            frequencies: [[0.0; 100]; 4],
            information: [0.0; 100],
            min_scores: [0.0; 100],
            max_scores: [0.0; 100],
            length: json.length,
            core_start: json.core_start,
            core_end: json.core_end,
            total_min: 0.0,
            total_max: 0.0,
            core_min: 0.0,
            core_max: 0.0,
            css_cutoff: 0.0,
            mss_cutoff: 0.0,
        };

        // Copy frequencies
        for (i, &freq) in json.frequencies[&'A'].iter().enumerate() {
            pwm.frequencies[0][i] = freq;
        }
        for (i, &freq) in json.frequencies[&'C'].iter().enumerate() {
            pwm.frequencies[1][i] = freq;
        }
        for (i, &freq) in json.frequencies[&'G'].iter().enumerate() {
            pwm.frequencies[2][i] = freq;
        }
        for (i, &freq) in json.frequencies[&'T'].iter().enumerate() {
            pwm.frequencies[3][i] = freq;
        }

        // Copy information vector
        for (i, &info) in json.information.iter().enumerate() {
            pwm.information[i] = info;
        }

        // Calculate min/max scores
        for i in 0..pwm.length {
            let freqs = [
                pwm.frequencies[0][i],
                pwm.frequencies[1][i],
                pwm.frequencies[2][i],
                pwm.frequencies[3][i]
            ];

            let min_freq = freqs.iter().copied().fold(f64::INFINITY, f64::min);
            let max_freq = freqs.iter().copied().fold(f64::NEG_INFINITY, f64::max);

            pwm.min_scores[i] = if min_freq > 0.0 {
                pwm.information[i] * (4.0 * min_freq).log2()
            } else {
                -10.0
            };

            pwm.max_scores[i] = if max_freq > 0.0 {
                pwm.information[i] * (4.0 * max_freq).log2()
            } else {
                0.0
            };

            pwm.total_min += pwm.min_scores[i];
            pwm.total_max += pwm.max_scores[i];
        }

        // Calculate core min/max
        for i in pwm.core_start..pwm.core_end.min(pwm.length) {
            pwm.core_min += pwm.min_scores[i];
            pwm.core_max += pwm.max_scores[i];
        }

        // Set cutoffs
        match cutoff_strategy {
            "minFN" => {
                pwm.css_cutoff = 0.60;
                pwm.mss_cutoff = 0.70;
            },
            "minFP" => {
                pwm.css_cutoff = 0.85;
                pwm.mss_cutoff = 0.90;
            },
            _ => { // minSum
                pwm.css_cutoff = 0.75;
                pwm.mss_cutoff = 0.80;
            }
        }

        pwm
    }

    #[inline]
    fn nt_to_idx(nt: u8) -> Option<usize> {
        match nt.to_ascii_uppercase() {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' => Some(3),
            _ => None,
        }
    }

    #[inline]
    fn calculate_score(&self, sequence: &[u8], start: usize) -> f64 {
        let mut score = 0.0;

        for i in 0..self.length {
            if start + i >= sequence.len() {
                break;
            }

            if let Some(idx) = Self::nt_to_idx(sequence[start + i]) {
                let freq = self.frequencies[idx][i];
                if freq > 0.0 {
                    score += self.information[i] * (4.0 * freq).log2();
                } else {
                    score += self.min_scores[i];
                }
            }
        }

        score
    }

    #[inline]
    fn calculate_css(&self, pentamer: &[u8]) -> f64 {
        let mut score = 0.0;

        for (i, &nt) in pentamer.iter().take(5).enumerate() {
            if let Some(idx) = Self::nt_to_idx(nt) {
                let pwm_pos = self.core_start + i;
                if pwm_pos < self.length {
                    let freq = self.frequencies[idx][pwm_pos];
                    if freq > 0.0 {
                        score += self.information[pwm_pos] * (4.0 * freq).log2();
                    } else {
                        score += self.min_scores[pwm_pos];
                    }
                }
            }
        }

        // Normalize to [0, 1]
        if self.core_max - self.core_min > 0.0 {
            score = (score - self.core_min) / (self.core_max - self.core_min);
        }

        score.max(0.0).min(1.0)
    }

    #[inline]
    fn calculate_mss(&self, sequence: &[u8], start: usize) -> f64 {
        let score = self.calculate_score(sequence, start);

        // Normalize to [0, 1]
        let normalized = if self.total_max - self.total_min > 0.0 {
            (score - self.total_min) / (self.total_max - self.total_min)
        } else {
            0.0
        };

        normalized.max(0.0).min(1.0)
    }

    fn pentamer_to_hash(pentamer: &[u8]) -> Option<usize> {
        let mut hash = 0;
        for &nt in pentamer.iter().take(5) {
            hash *= 4;
            match Self::nt_to_idx(nt) {
                Some(idx) => hash += idx,
                None => return None,
            }
        }
        Some(hash)
    }

    fn scan_sequence(&self, sequence: &[u8], _seq_id: &str) -> Vec<Match> {
        let mut matches = Vec::new();
        let mut pentamer_positions: HashMap<usize, Vec<usize>> = HashMap::new();

        // Build pentamer hash table
        for i in 0..=(sequence.len().saturating_sub(5)) {
            let pentamer = &sequence[i..i+5];
            if let Some(hash) = Self::pentamer_to_hash(pentamer) {
                pentamer_positions.entry(hash).or_insert_with(Vec::new).push(i);
            }
        }

        // Scan using pentamers that pass CSS threshold
        for (&hash, positions) in &pentamer_positions {
            // Reconstruct pentamer from hash
            let mut pentamer = [0u8; 5];
            let mut temp_hash = hash;
            for i in (0..5).rev() {
                pentamer[i] = b"ACGT"[temp_hash % 4];
                temp_hash /= 4;
            }

            // Check CSS
            let css = self.calculate_css(&pentamer);
            if css < self.css_cutoff {
                continue;
            }

            // Check all positions of this pentamer
            for &pos in positions {
                // Try different alignments with PWM core
                for core_offset in 0..5 {
                    let pwm_start = pos.saturating_sub(self.core_start + core_offset);

                    if pwm_start + self.length <= sequence.len() {
                        // Calculate MSS
                        let mss = self.calculate_mss(sequence, pwm_start);

                        if mss >= self.mss_cutoff {
                            matches.push(Match {
                                start: pwm_start,
                                end: pwm_start + self.length,
                                strand: '+',
                                css,
                                mss,
                                score: css * 0.3 + mss * 0.7,
                            });
                        }
                    }
                }
            }
        }

        matches
    }
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence.iter().rev().map(|&nt| {
        match nt.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            _ => b'N',
        }
    }).collect()
}

fn load_pwm(filename: &str, cutoff_strategy: &str) -> Result<PWM, Box<dyn std::error::Error>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let json: PWMJson = serde_json::from_reader(reader)?;
    Ok(PWM::from_json(json, cutoff_strategy))
}

fn process_fasta_file(filename: &str, pwm: &PWM, min_score: f64, forward_only: bool) -> Result<Vec<(String, Vec<Match>)>, Box<dyn std::error::Error>> {
    let file = File::open(filename)?;
    let reader: Box<dyn BufRead> = if filename.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut sequences = Vec::new();
    let mut current_id = String::new();
    let mut current_seq = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            // Process previous sequence
            if !current_id.is_empty() && !current_seq.is_empty() {
                sequences.push((current_id.clone(), current_seq.clone()));
            }

            // Start new sequence
            current_id = line[1..].split_whitespace().next().unwrap_or("").to_string();
            current_seq.clear();
        } else {
            current_seq.extend(line.as_bytes().iter().filter(|&&c| c != b'\n'));
        }
    }

    // Process last sequence
    if !current_id.is_empty() && !current_seq.is_empty() {
        sequences.push((current_id, current_seq));
    }

    // Process sequences in parallel with progress tracking
    let total_sequences = sequences.len();
    eprintln!("Found {} sequences, starting parallel processing...", total_sequences);

    let results: Vec<_> = sequences.into_par_iter().enumerate().map(|(idx, (seq_id, sequence))| {
        eprintln!("[{}/{}] Processing {} ({} bp)...",
                 idx + 1, total_sequences, seq_id, sequence.len());

        let mut all_matches = Vec::new();

        // Forward strand
        let forward_matches = pwm.scan_sequence(&sequence, &seq_id);
        all_matches.extend(forward_matches.into_iter().filter(|m| m.score >= min_score));

        // Reverse strand
        if !forward_only {
            let rev_comp = reverse_complement(&sequence);
            let mut rev_matches = pwm.scan_sequence(&rev_comp, &seq_id);

            // Adjust coordinates for reverse strand
            for match_ in &mut rev_matches {
                match_.strand = '-';
                let start = sequence.len() - match_.end;
                let end = sequence.len() - match_.start;
                match_.start = start;
                match_.end = end;
            }

            all_matches.extend(rev_matches.into_iter().filter(|m| m.score >= min_score));
        }

        // Sort by score (highest first)
        all_matches.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());

        eprintln!("[{}/{}] {} found {} matches",
                 idx + 1, total_sequences, seq_id, all_matches.len());

        (seq_id, all_matches)
    }).collect();

    Ok(results)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // Set number of threads
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap();
    }

    // Load PWM
    let pwm = load_pwm(&args.pwm_file, &args.cutoff_strategy)?;
    eprintln!("PWM loaded: length={}, core={}-{}", pwm.length, pwm.core_start, pwm.core_end);
    eprintln!("Cutoffs: CSS={:.2}, MSS={:.2}", pwm.css_cutoff, pwm.mss_cutoff);

    // Process genome
    let results = process_fasta_file(&args.genome_file, &pwm, args.min_score, args.forward_only)?;

    // Open output
    let mut output: Box<dyn Write> = match args.output_file {
        Some(filename) => Box::new(File::create(filename)?),
        None => Box::new(std::io::stdout()),
    };

    // Write GFF3 header
    writeln!(output, "##gff-version 3")?;

    // Write results
    let mut total_matches = 0;
    for (seq_id, matches) in results {
        for match_ in matches {
            writeln!(output, "{}\tMATCH_Rust\tpromoter\t{}\t{}\t{:.3}\t{}\t.\tID=promoter_{}_{};css={:.3};mss={:.3}",
                seq_id, match_.start + 1, match_.end, match_.score, match_.strand,
                seq_id, match_.start, match_.css, match_.mss)?;
            total_matches += 1;
        }
    }

    eprintln!("Total matches found: {}", total_matches);

    Ok(())
}