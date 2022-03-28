fn main() {
    for _ in 0..1_000 {
        if let Some(largest) = find_largest(1_000_000_000) {
            println!("{largest}");
        }
    }
}

fn find_largest(end: usize) -> Option<usize> {
    let mut n = end - if end % 2 == 0 { 1 } else { 0 };

    while n >= 2 {
        if is_prime(n) {
            return Some(n);
        }
        n -= 2
    }

    None
}

//#[inline]
fn is_prime(n: usize) -> bool {
    for i in 2..(n as f32).sqrt() as usize {
        if n % i == 0 {
            return false;
        }
    }

    true
}
