//! Sortieren in genau der Reihenfolge von bwa (`ks_introsort` aus klib, MIT-Lizenz).
//!
//! Das Verfahren ist nicht stabil: Gleichrangige Einträge landen in einer Reihenfolge, die vom
//! Ablauf abhängt. bwa verlässt sich an mehreren Stellen stillschweigend darauf, etwa wenn von
//! mehreren gleich schweren Ketten "die erste" behalten wird. Wer dieselben Ergebnisse will,
//! muss deshalb genau so sortieren. `lt(a, b)` heißt: a gehört vor b.

fn insertion_sort<T>(a: &mut [T], lt: &impl Fn(&T, &T) -> bool) {
    for i in 1..a.len() {
        let mut j = i;
        while j > 0 && lt(&a[j], &a[j - 1]) {
            a.swap(j, j - 1);
            j -= 1;
        }
    }
}

fn comb_sort<T>(a: &mut [T], lt: &impl Fn(&T, &T) -> bool) {
    const SHRINK: f64 = 1.247_330_950_103_978_6;
    let n = a.len();
    let mut gap = n;
    loop {
        if gap > 2 {
            gap = (gap as f64 / SHRINK) as usize;
            if gap == 9 || gap == 10 {
                gap = 11;
            }
        }
        let mut swapped = false;
        for i in 0..n.saturating_sub(gap) {
            if lt(&a[i + gap], &a[i]) {
                a.swap(i, i + gap);
                swapped = true;
            }
        }
        if !(swapped || gap > 2) {
            break;
        }
    }
    if gap != 1 {
        insertion_sort(a, lt);
    }
}

/// Sortiert `a` wie `ks_introsort`: Quicksort mit Median aus drei, Teile bis 17 Einträge bleiben
/// liegen und werden am Ende gemeinsam durch Einfügen sortiert.
pub fn introsort<T: Clone>(a: &mut [T], lt: impl Fn(&T, &T) -> bool) {
    let n = a.len();
    if n < 1 {
        return;
    }
    if n == 2 {
        if lt(&a[1], &a[0]) {
            a.swap(0, 1);
        }
        return;
    }
    let mut d = 2u32;
    while (1usize << d) < n {
        d += 1;
    }
    d <<= 1;
    let mut stack: Vec<(isize, isize, u32)> = Vec::new();
    let (mut s, mut t) = (0isize, n as isize - 1);
    loop {
        if s < t {
            d -= 1;
            if d == 0 {
                comb_sort(&mut a[s as usize..=t as usize], &lt);
                t = s;
                continue;
            }
            let (mut i, mut j) = (s, t);
            let mut k = i + ((j - i) >> 1) + 1;
            if lt(&a[k as usize], &a[i as usize]) {
                if lt(&a[k as usize], &a[j as usize]) {
                    k = j;
                }
            } else {
                k = if lt(&a[j as usize], &a[i as usize]) {
                    i
                } else {
                    j
                };
            }
            let rp = a[k as usize].clone();
            if k != t {
                a.swap(k as usize, t as usize);
            }
            loop {
                loop {
                    i += 1;
                    if !lt(&a[i as usize], &rp) {
                        break;
                    }
                }
                loop {
                    j -= 1;
                    if !(i <= j && lt(&rp, &a[j as usize])) {
                        break;
                    }
                }
                if j <= i {
                    break;
                }
                a.swap(i as usize, j as usize);
            }
            a.swap(i as usize, t as usize);
            if i - s > t - i {
                if i - s > 16 {
                    stack.push((s, i - 1, d));
                }
                s = if t - i > 16 { i + 1 } else { t };
            } else {
                if t - i > 16 {
                    stack.push((i + 1, t, d));
                }
                t = if i - s > 16 { i - 1 } else { s };
            }
        } else if let Some((left, right, depth)) = stack.pop() {
            s = left;
            t = right;
            d = depth;
        } else {
            insertion_sort(a, &lt);
            return;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sorts_like_a_stable_sort_when_keys_are_distinct_and_keeps_all_entries_otherwise() {
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for n in [0usize, 1, 2, 3, 16, 17, 18, 40, 500, 5000] {
            let mut v: Vec<u64> = (0..n).map(|_| next()).collect();
            let mut want = v.clone();
            want.sort_unstable();
            introsort(&mut v, |a, b| a < b);
            assert_eq!(v, want, "n = {n}");
            // viele gleiche Schlüssel: sortiert nach Schlüssel, kein Eintrag geht verloren
            let mut v: Vec<(u8, usize)> = (0..n).map(|i| ((next() % 4) as u8, i)).collect();
            introsort(&mut v, |a, b| a.0 < b.0);
            assert!(v.windows(2).all(|w| w[0].0 <= w[1].0));
            let mut ids: Vec<usize> = v.iter().map(|x| x.1).collect();
            ids.sort_unstable();
            assert_eq!(ids, (0..n).collect::<Vec<_>>());
        }
    }

    /// Festgehaltene Reihenfolge bei Gleichstand, von Hand am Verfahren nachvollzogen:
    /// drei gleiche Schlüssel, der Median-Schritt wählt den letzten als Trennwert, die Teilung
    /// bleibt beim zweiten stehen und tauscht ihn mit dem letzten.
    #[test]
    fn tie_order_is_the_one_of_the_method_not_the_input_order() {
        let mut v = vec![(1, 'a'), (1, 'b'), (1, 'c')];
        introsort(&mut v, |a, b| a.0 < b.0);
        assert_eq!(v.iter().map(|x| x.1).collect::<String>(), "acb");
    }
}
