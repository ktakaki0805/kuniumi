//! Evaluates `kuniumi_math::det` functions on inputs read from standard input,
//! for `tools/det_error.py`, which measures their error against correctly
//! rounded values (ticket M1-01).
//!
//! Each input line is a function name followed by one or two argument bit
//! patterns in hexadecimal (`atan2` takes `y x`, `pow` takes `x y`). Each
//! output line is the bit pattern of the result in hexadecimal.

use std::io::{self, BufRead, BufWriter, Write};

use kuniumi_math::det;

fn parse_bits(word: Option<&str>, line: &str) -> f64 {
    let word = word.unwrap_or_else(|| panic!("missing argument in {line:?}"));
    let bits = u64::from_str_radix(word, 16).unwrap_or_else(|e| panic!("{e} in {line:?}"));
    f64::from_bits(bits)
}

fn main() -> io::Result<()> {
    let stdin = io::stdin();
    let mut out = BufWriter::new(io::stdout().lock());
    for line in stdin.lock().lines() {
        let line = line?;
        let mut words = line.split_whitespace();
        let Some(name) = words.next() else { continue };
        let a = parse_bits(words.next(), &line);
        let result = match name {
            "sin" => det::sin(a),
            "cos" => det::cos(a),
            "tan" => det::tan(a),
            "atan" => det::atan(a),
            "asin" => det::asin(a),
            "acos" => det::acos(a),
            "exp" => det::exp(a),
            "ln" => det::ln(a),
            "atan2" => det::atan2(a, parse_bits(words.next(), &line)),
            "pow" => det::pow(a, parse_bits(words.next(), &line)),
            _ => panic!("unknown function in {line:?}"),
        };
        writeln!(out, "{:016x}", result.to_bits())?;
    }
    out.flush()
}
