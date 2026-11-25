//! Vektoren mit 16 Bytes ohne Vorzeichen (`U8x16`) oder 8 Werten mit 16 Bit (`I16x8`), für das
//! gestreifte Smith-Waterman. Auf x86_64 mit SSE2, auf aarch64 mit NEON (beide sind auf diesen
//! Prozessoren immer da), sonst mit gewöhnlichen Feldern. Jede Operation liefert auf allen drei
//! Wegen dasselbe; die Tests prüfen das gegen die Felder.

/// Ein Vektor fester Breite. Werte von `I16x8` sind in allen Rechnungen nie negativ; `subs`
/// begrenzt dann unten bei 0.
pub trait Lanes: Copy {
    type Elem: Copy + Ord + Default;
    const N: usize;
    fn splat(x: Self::Elem) -> Self;
    fn zero() -> Self;
    /// Gesättigte Addition (Bytes: oben bei 255, 16 Bit: oben bei 32767).
    fn adds(self, o: Self) -> Self;
    /// Subtraktion, unten bei 0 begrenzt.
    fn subs(self, o: Self) -> Self;
    fn max(self, o: Self) -> Self;
    /// Ob in irgendeiner Spur `self > o` gilt.
    fn any_gt(self, o: Self) -> bool;
    /// Alle Spuren um eine nach oben verschoben, Spur 0 wird 0, die letzte fällt heraus.
    fn shift_in_zero(self) -> Self;
    /// Größter Wert aller Spuren.
    fn hmax(self) -> Self::Elem;
    /// Aus den ersten `N` Werten von `v`.
    fn load(v: &[Self::Elem]) -> Self;
    /// In die ersten `N` Werte von `out`.
    fn store(self, out: &mut [Self::Elem]);
}

#[cfg(target_arch = "x86_64")]
mod imp {
    use super::Lanes;
    use std::arch::x86_64::*;

    #[derive(Clone, Copy)]
    pub struct U8x16(__m128i);
    #[derive(Clone, Copy)]
    pub struct I16x8(__m128i);

    // SAFETY (für alle Blöcke hier): SSE2 gehört zu jedem x86_64-Prozessor; die Lade- und
    // Speicherbefehle sind die ungebundenen (`loadu`/`storeu`), die Längen werden vorher geprüft.
    impl Lanes for U8x16 {
        type Elem = u8;
        const N: usize = 16;
        #[inline(always)]
        fn splat(x: u8) -> Self {
            unsafe { Self(_mm_set1_epi8(x as i8)) }
        }
        #[inline(always)]
        fn zero() -> Self {
            unsafe { Self(_mm_setzero_si128()) }
        }
        #[inline(always)]
        fn adds(self, o: Self) -> Self {
            unsafe { Self(_mm_adds_epu8(self.0, o.0)) }
        }
        #[inline(always)]
        fn subs(self, o: Self) -> Self {
            unsafe { Self(_mm_subs_epu8(self.0, o.0)) }
        }
        #[inline(always)]
        fn max(self, o: Self) -> Self {
            unsafe { Self(_mm_max_epu8(self.0, o.0)) }
        }
        #[inline(always)]
        fn any_gt(self, o: Self) -> bool {
            unsafe {
                let d = _mm_subs_epu8(self.0, o.0);
                _mm_movemask_epi8(_mm_cmpeq_epi8(d, _mm_setzero_si128())) != 0xffff
            }
        }
        #[inline(always)]
        fn shift_in_zero(self) -> Self {
            unsafe { Self(_mm_slli_si128::<1>(self.0)) }
        }
        #[inline(always)]
        fn hmax(self) -> u8 {
            unsafe {
                let mut m = _mm_max_epu8(self.0, _mm_srli_si128::<8>(self.0));
                m = _mm_max_epu8(m, _mm_srli_si128::<4>(m));
                m = _mm_max_epu8(m, _mm_srli_si128::<2>(m));
                m = _mm_max_epu8(m, _mm_srli_si128::<1>(m));
                (_mm_cvtsi128_si32(m) & 0xff) as u8
            }
        }
        #[inline(always)]
        fn load(v: &[u8]) -> Self {
            assert!(v.len() >= 16);
            unsafe { Self(_mm_loadu_si128(v.as_ptr().cast())) }
        }
        #[inline(always)]
        fn store(self, out: &mut [u8]) {
            assert!(out.len() >= 16);
            unsafe { _mm_storeu_si128(out.as_mut_ptr().cast(), self.0) }
        }
    }

    impl Lanes for I16x8 {
        type Elem = i16;
        const N: usize = 8;
        #[inline(always)]
        fn splat(x: i16) -> Self {
            unsafe { Self(_mm_set1_epi16(x)) }
        }
        #[inline(always)]
        fn zero() -> Self {
            unsafe { Self(_mm_setzero_si128()) }
        }
        #[inline(always)]
        fn adds(self, o: Self) -> Self {
            unsafe { Self(_mm_adds_epi16(self.0, o.0)) }
        }
        #[inline(always)]
        fn subs(self, o: Self) -> Self {
            unsafe { Self(_mm_subs_epu16(self.0, o.0)) }
        }
        #[inline(always)]
        fn max(self, o: Self) -> Self {
            unsafe { Self(_mm_max_epi16(self.0, o.0)) }
        }
        #[inline(always)]
        fn any_gt(self, o: Self) -> bool {
            unsafe { _mm_movemask_epi8(_mm_cmpgt_epi16(self.0, o.0)) != 0 }
        }
        #[inline(always)]
        fn shift_in_zero(self) -> Self {
            unsafe { Self(_mm_slli_si128::<2>(self.0)) }
        }
        #[inline(always)]
        fn hmax(self) -> i16 {
            unsafe {
                let mut m = _mm_max_epi16(self.0, _mm_srli_si128::<8>(self.0));
                m = _mm_max_epi16(m, _mm_srli_si128::<4>(m));
                m = _mm_max_epi16(m, _mm_srli_si128::<2>(m));
                (_mm_cvtsi128_si32(m) & 0xffff) as u16 as i16
            }
        }
        #[inline(always)]
        fn load(v: &[i16]) -> Self {
            assert!(v.len() >= 8);
            unsafe { Self(_mm_loadu_si128(v.as_ptr().cast())) }
        }
        #[inline(always)]
        fn store(self, out: &mut [i16]) {
            assert!(out.len() >= 8);
            unsafe { _mm_storeu_si128(out.as_mut_ptr().cast(), self.0) }
        }
    }
}

#[cfg(target_arch = "aarch64")]
mod imp {
    use super::Lanes;
    use std::arch::aarch64::*;

    #[derive(Clone, Copy)]
    pub struct U8x16(uint8x16_t);
    #[derive(Clone, Copy)]
    pub struct I16x8(int16x8_t);

    // SAFETY (für alle Blöcke hier): NEON gehört zu jedem aarch64-Prozessor; die Längen beim
    // Laden und Speichern werden vorher geprüft.
    impl Lanes for U8x16 {
        type Elem = u8;
        const N: usize = 16;
        #[inline(always)]
        fn splat(x: u8) -> Self {
            unsafe { Self(vdupq_n_u8(x)) }
        }
        #[inline(always)]
        fn zero() -> Self {
            unsafe { Self(vdupq_n_u8(0)) }
        }
        #[inline(always)]
        fn adds(self, o: Self) -> Self {
            unsafe { Self(vqaddq_u8(self.0, o.0)) }
        }
        #[inline(always)]
        fn subs(self, o: Self) -> Self {
            unsafe { Self(vqsubq_u8(self.0, o.0)) }
        }
        #[inline(always)]
        fn max(self, o: Self) -> Self {
            unsafe { Self(vmaxq_u8(self.0, o.0)) }
        }
        #[inline(always)]
        fn any_gt(self, o: Self) -> bool {
            unsafe { vmaxvq_u8(vcgtq_u8(self.0, o.0)) != 0 }
        }
        #[inline(always)]
        fn shift_in_zero(self) -> Self {
            unsafe { Self(vextq_u8::<15>(vdupq_n_u8(0), self.0)) }
        }
        #[inline(always)]
        fn hmax(self) -> u8 {
            unsafe { vmaxvq_u8(self.0) }
        }
        #[inline(always)]
        fn load(v: &[u8]) -> Self {
            assert!(v.len() >= 16);
            unsafe { Self(vld1q_u8(v.as_ptr())) }
        }
        #[inline(always)]
        fn store(self, out: &mut [u8]) {
            assert!(out.len() >= 16);
            unsafe { vst1q_u8(out.as_mut_ptr(), self.0) }
        }
    }

    impl Lanes for I16x8 {
        type Elem = i16;
        const N: usize = 8;
        #[inline(always)]
        fn splat(x: i16) -> Self {
            unsafe { Self(vdupq_n_s16(x)) }
        }
        #[inline(always)]
        fn zero() -> Self {
            unsafe { Self(vdupq_n_s16(0)) }
        }
        #[inline(always)]
        fn adds(self, o: Self) -> Self {
            unsafe { Self(vqaddq_s16(self.0, o.0)) }
        }
        #[inline(always)]
        fn subs(self, o: Self) -> Self {
            unsafe {
                Self(vreinterpretq_s16_u16(vqsubq_u16(
                    vreinterpretq_u16_s16(self.0),
                    vreinterpretq_u16_s16(o.0),
                )))
            }
        }
        #[inline(always)]
        fn max(self, o: Self) -> Self {
            unsafe { Self(vmaxq_s16(self.0, o.0)) }
        }
        #[inline(always)]
        fn any_gt(self, o: Self) -> bool {
            unsafe { vmaxvq_u16(vcgtq_s16(self.0, o.0)) != 0 }
        }
        #[inline(always)]
        fn shift_in_zero(self) -> Self {
            unsafe { Self(vextq_s16::<7>(vdupq_n_s16(0), self.0)) }
        }
        #[inline(always)]
        fn hmax(self) -> i16 {
            unsafe { vmaxvq_s16(self.0) }
        }
        #[inline(always)]
        fn load(v: &[i16]) -> Self {
            assert!(v.len() >= 8);
            unsafe { Self(vld1q_s16(v.as_ptr())) }
        }
        #[inline(always)]
        fn store(self, out: &mut [i16]) {
            assert!(out.len() >= 8);
            unsafe { vst1q_s16(out.as_mut_ptr(), self.0) }
        }
    }
}

/// Dieselben Operationen mit gewöhnlichen Feldern: auf anderen Prozessoren im Einsatz, überall
/// als Maßstab für die Tests.
pub mod plain {
    use super::Lanes;

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub struct U8x16(pub [u8; 16]);
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub struct I16x8(pub [i16; 8]);

    macro_rules! plain_lanes {
        ($t:ident, $e:ty, $n:expr, $subs:expr) => {
            impl Lanes for $t {
                type Elem = $e;
                const N: usize = $n;
                #[inline(always)]
                fn splat(x: $e) -> Self {
                    Self([x; $n])
                }
                #[inline(always)]
                fn zero() -> Self {
                    Self([0; $n])
                }
                #[inline(always)]
                fn adds(self, o: Self) -> Self {
                    Self(std::array::from_fn(|k| self.0[k].saturating_add(o.0[k])))
                }
                #[inline(always)]
                fn subs(self, o: Self) -> Self {
                    let f: fn($e, $e) -> $e = $subs;
                    Self(std::array::from_fn(|k| f(self.0[k], o.0[k])))
                }
                #[inline(always)]
                fn max(self, o: Self) -> Self {
                    Self(std::array::from_fn(|k| self.0[k].max(o.0[k])))
                }
                #[inline(always)]
                fn any_gt(self, o: Self) -> bool {
                    (0..$n).any(|k| self.0[k] > o.0[k])
                }
                #[inline(always)]
                fn shift_in_zero(self) -> Self {
                    Self(std::array::from_fn(|k| if k == 0 { 0 } else { self.0[k - 1] }))
                }
                #[inline(always)]
                fn hmax(self) -> $e {
                    self.0.iter().copied().max().unwrap_or(0)
                }
                #[inline(always)]
                fn load(v: &[$e]) -> Self {
                    Self(std::array::from_fn(|k| v[k]))
                }
                #[inline(always)]
                fn store(self, out: &mut [$e]) {
                    out[..$n].copy_from_slice(&self.0);
                }
            }
        };
    }
    plain_lanes!(U8x16, u8, 16, |a, b| a.saturating_sub(b));
    plain_lanes!(I16x8, i16, 8, |a, b| a.saturating_sub(b).max(0));
}

#[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
pub use imp::{I16x8, U8x16};
#[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
pub use plain::{I16x8, U8x16};

#[cfg(test)]
mod tests {
    use super::*;

    fn next(s: &mut u64) -> u64 {
        *s ^= *s << 13;
        *s ^= *s >> 7;
        *s ^= *s << 17;
        *s
    }

    fn check<V: Lanes, P: Lanes<Elem = V::Elem>>(s: &mut u64, make: impl Fn(u64) -> V::Elem)
    where
        V::Elem: std::fmt::Debug,
    {
        let n = V::N;
        let get = |v: V| {
            let mut o = vec![V::Elem::default(); n];
            v.store(&mut o);
            o
        };
        let getp = |v: P| {
            let mut o = vec![V::Elem::default(); n];
            v.store(&mut o);
            o
        };
        for _ in 0..20_000 {
            let a: Vec<V::Elem> = (0..n).map(|_| make(next(s))).collect();
            let b: Vec<V::Elem> = (0..n).map(|_| make(next(s))).collect();
            let (va, vb) = (V::load(&a), V::load(&b));
            let (pa, pb) = (P::load(&a), P::load(&b));
            assert_eq!(get(va.adds(vb)), getp(pa.adds(pb)));
            assert_eq!(get(va.subs(vb)), getp(pa.subs(pb)));
            assert_eq!(get(va.max(vb)), getp(pa.max(pb)));
            assert_eq!(get(va.shift_in_zero()), getp(pa.shift_in_zero()));
            assert_eq!(va.any_gt(vb), pa.any_gt(pb));
            assert_eq!(va.hmax(), pa.hmax());
            assert_eq!(get(V::splat(a[0])), getp(P::splat(a[0])));
            assert_eq!(get(V::zero()), getp(P::zero()));
        }
    }

    #[test]
    fn vectors_equal_plain_arrays() {
        let mut s = 0x9e37_79b9_7f4a_7c15;
        check::<U8x16, plain::U8x16>(&mut s, |x| x as u8);
        // Kleine Werte, damit Gleichheit und Grenzen oft vorkommen.
        check::<U8x16, plain::U8x16>(&mut s, |x| (x % 6) as u8);
        // 16 Bit: nie negativ, auch nahe der oberen Grenze.
        check::<I16x8, plain::I16x8>(&mut s, |x| (x % 32768) as i16);
        check::<I16x8, plain::I16x8>(&mut s, |x| (x % 6) as i16);
        check::<I16x8, plain::I16x8>(&mut s, |x| (32760 + x % 8) as i16);
    }
}
