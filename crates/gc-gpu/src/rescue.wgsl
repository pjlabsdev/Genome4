// Gestreifter Byte-Lauf der lokalen Ausrichtung (wie gc_align::local::kernel_reference),
// eine Ausrichtung je Aufruf. 16 Spuren je Vektor, jede Spur ein Wert 0..255 in einem u32.

struct Params {
    n_jobs: u32,
    matching: u32,
    mismatch: u32,
    ambiguous: u32,
    e_del: u32,
    oe_del: u32,
    e_ins: u32,
    oe_ins: u32,
}

struct Job {
    q_off: u32,
    q_len: u32,
    t_off: u32,
    t_len: u32,
    rows_off: u32,
    pad0: u32,
    pad1: u32,
    pad2: u32,
}

@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read> jobs: array<Job>;
@group(0) @binding(2) var<storage, read> bases: array<u32>;
@group(0) @binding(3) var<storage, read_write> result: array<u32>;
@group(0) @binding(4) var<storage, read_write> rows: array<atomic<u32>>;

const L: u32 = 16u;
const MAXV: u32 = 256u;

var<private> h0: array<u32, 256>;
var<private> h1: array<u32, 256>;
var<private> e_row: array<u32, 256>;
var<private> h_max: array<u32, 256>;
var<private> q: array<u32, 256>;

fn base_at(i: u32) -> u32 {
    return (bases[i >> 2u] >> ((i & 3u) * 8u)) & 0xffu;
}

fn sat_add(a: u32, b: u32) -> u32 {
    return min(a + b, 255u);
}

fn sat_sub(a: u32, b: u32) -> u32 {
    return select(0u, a - b, a > b);
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let id = gid.x;
    if (id >= P.n_jobs) {
        return;
    }
    let job = jobs[id];
    let qlen = job.q_len;
    let slen = (qlen + L - 1u) / L;
    let shift = max(P.mismatch, P.ambiguous);
    for (var p = 0u; p < slen * L; p++) {
        if (p < qlen) {
            q[p] = base_at(job.q_off + p);
        } else {
            q[p] = 9u;
        }
        h0[p] = 0u;
        h1[p] = 0u;
        e_row[p] = 0u;
        h_max[p] = 0u;
    }
    var gmax = 0u;
    var te = -1;
    var cur = 0u; // 0: h0 ist die vorige Zeile, 1: h1
    for (var i = 0u; i < job.t_len; i++) {
        let b = min(base_at(job.t_off + i), 4u);
        var h: array<u32, 16>;
        var f: array<u32, 16>;
        var row_max: array<u32, 16>;
        h[0] = 0u;
        for (var k = 1u; k < L; k++) {
            let x = (slen - 1u) * L + k - 1u;
            h[k] = select(h1[x], h0[x], cur == 0u);
        }
        for (var k = 0u; k < L; k++) {
            f[k] = 0u;
            row_max[k] = 0u;
        }
        for (var j = 0u; j < slen; j++) {
            for (var k = 0u; k < L; k++) {
                let x = j * L + k;
                let qb = q[j + k * slen];
                var s: u32;
                if (qb == 9u) {
                    s = shift;
                } else if (b > 3u || qb > 3u) {
                    s = shift - P.ambiguous;
                } else if (b == qb) {
                    s = shift + P.matching;
                } else {
                    s = shift - P.mismatch;
                }
                let t = max(sat_sub(sat_add(h[k], s), shift), e_row[x]);
                let v = max(t, f[k]);
                row_max[k] = max(row_max[k], v);
                if (cur == 0u) {
                    h[k] = h0[x];
                    h1[x] = v;
                } else {
                    h[k] = h1[x];
                    h0[x] = v;
                }
                e_row[x] = max(sat_sub(e_row[x], P.e_del), sat_sub(v, P.oe_del));
                f[k] = max(sat_sub(f[k], P.e_ins), sat_sub(t, P.oe_ins));
            }
        }
        // Nachkorrektur: die offene Insertion von Spur zu Spur weiterreichen.
        var done = false;
        for (var round = 0u; round < 16u && !done; round++) {
            for (var k = L - 1u; k > 0u; k--) {
                f[k] = f[k - 1u];
            }
            f[0] = 0u;
            for (var j = 0u; j < slen; j++) {
                var any = false;
                for (var k = 0u; k < L; k++) {
                    let x = j * L + k;
                    var v: u32;
                    if (cur == 0u) {
                        v = max(h1[x], f[k]);
                        h1[x] = v;
                    } else {
                        v = max(h0[x], f[k]);
                        h0[x] = v;
                    }
                    f[k] = sat_sub(f[k], P.e_ins);
                    any = any || (f[k] > sat_sub(v, P.oe_ins));
                }
                if (!any) {
                    done = true;
                    break;
                }
            }
        }
        var imax = 0u;
        for (var k = 0u; k < L; k++) {
            imax = max(imax, row_max[k]);
        }
        let ri = job.rows_off + i;
        atomicOr(&rows[ri >> 2u], imax << ((ri & 3u) * 8u));
        if (imax > gmax) {
            gmax = imax;
            te = i32(i);
            for (var x = 0u; x < slen * L; x++) {
                h_max[x] = select(h0[x], h1[x], cur == 0u);
            }
        }
        cur = 1u - cur;
    }
    var mx = -1;
    var qe = -1;
    for (var j = 0u; j < slen; j++) {
        for (var k = 0u; k < L; k++) {
            let v = i32(h_max[j * L + k]);
            let p = i32(j + k * slen);
            if (v > mx || (v == mx && p < qe)) {
                mx = v;
                qe = p;
            }
        }
    }
    result[id * 4u] = gmax;
    result[id * 4u + 1u] = bitcast<u32>(te);
    result[id * 4u + 2u] = bitcast<u32>(qe);
}
