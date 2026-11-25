//! Kleiner, fester Zufallsgenerator (xoshiro256**), damit ein Lauf mit derselben Saat auf
//! jedem System dieselben Reads erzeugt. Keine fremde Crate: die Folge darf sich nie durch ein
//! Update ändern.

pub struct Rng {
    s: [u64; 4],
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        // SplitMix64 füllt den Zustand aus der Saat.
        let mut x = seed;
        let mut next = || {
            x = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
            let mut z = x;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            z ^ (z >> 31)
        };
        Self {
            s: [next(), next(), next(), next()],
        }
    }

    pub fn next_u64(&mut self) -> u64 {
        let result = self.s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }

    /// Gleichverteilt in [0, 1).
    pub fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Gleichverteilt in [0, n).
    pub fn below(&mut self, n: u64) -> u64 {
        ((self.next_u64() as u128 * n as u128) >> 64) as u64
    }

    /// Normalverteilt (Box-Muller).
    pub fn gauss(&mut self, mean: f64, sd: f64) -> f64 {
        let u1 = 1.0 - self.unit();
        let u2 = self.unit();
        mean + sd * (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }

    /// Abstand bis zum nächsten Ereignis bei Rate `p` je Position (geometrisch, mindestens 1).
    pub fn gap(&mut self, p: f64) -> u64 {
        let u = 1.0 - self.unit();
        1 + (u.ln() / (1.0 - p).ln()) as u64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn same_seed_same_numbers() {
        let mut a = Rng::new(7);
        let mut b = Rng::new(7);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
        // Feste Werte: ändert sich die Folge, sind alte Messungen nicht mehr nachrechenbar.
        let mut c = Rng::new(1);
        assert_eq!(c.next_u64(), 0xb3f2_af6d_0fc7_10c5);
    }

    #[test]
    fn below_stays_in_range() {
        let mut r = Rng::new(3);
        for _ in 0..1000 {
            assert!(r.below(10) < 10);
        }
    }
}
