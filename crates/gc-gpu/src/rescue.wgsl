// Gestreifter Byte-Lauf der lokalen Ausrichtung (wie gc_align::local::kernel_reference).
// Eine Ausrichtung je Arbeitsgruppe (= eine Untergruppe von 32 Fäden): Faden k rechnet Spur
// k mod 16 des 16er-Vektors; die Fäden 16..31 rechnen dieselben Spuren noch einmal, so laufen
// alle Fäden einer Untergruppe denselben Weg. Weitergeben zwischen Spuren mit
// subgroupShuffleUp, "irgendeine Spur" mit subgroupAny, Zeilenmaximum mit subgroupMax.
diagnostic(off, subgroup_uniformity);

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

fn base_at(i: u32) -> u32 {
    return (bases[i >> 2u] >> ((i & 3u) * 8u)) & 0xffu;
}

fn sat_add(a: u32, b: u32) -> u32 {
    return min(a + b, 255u);
}

fn sat_sub(a: u32, b: u32) -> u32 {
    return select(0u, a - b, a > b);
}

@compute @workgroup_size(32)
fn main(
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(num_workgroups) nwg: vec3<u32>,
    @builtin(subgroup_invocation_id) sl: u32,
    @builtin(subgroup_size) sg: u32,
) {
    let id = wid.x + wid.y * nwg.x;
    if (id >= P.n_jobs) {
        return;
    }
    // Zu schmale Untergruppen: melden, der Prozessor rechnet.
    if (sg < 16u || sg % 16u != 0u) {
        if (sl == 0u) {
            result[id * 4u + 3u] = sg;
        }
        return;
    }
    let k = sl % L;
    let job = jobs[id];
    let qlen = job.q_len;
    let slen = (qlen + L - 1u) / L;
    let shift = max(P.mismatch, P.ambiguous);
    var q: array<u32, 16>;
    var h0: array<u32, 16>;
    var h1: array<u32, 16>;
    var e_row: array<u32, 16>;
    var h_max: array<u32, 16>;
    for (var j = 0u; j < slen; j++) {
        let p = j + k * slen;
        q[j] = select(9u, base_at(job.q_off + p), p < qlen);
        h0[j] = 0u;
        h1[j] = 0u;
        e_row[j] = 0u;
        h_max[j] = 0u;
    }
    var gmax = 0u;
    var te = -1;
    for (var i = 0u; i < job.t_len; i++) {
        let b = min(base_at(job.t_off + i), 4u);
        // Vorige Zeile, letzte Read-Position der Spur davor.
        let up = subgroupShuffleUp(h0[slen - 1u], 1u);
        var h = select(up, 0u, k == 0u);
        var f = 0u;
        var row_max = 0u;
        for (var j = 0u; j < slen; j++) {
            let qb = q[j];
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
            let t = max(sat_sub(sat_add(h, s), shift), e_row[j]);
            let v = max(t, f);
            row_max = max(row_max, v);
            h = h0[j];
            h1[j] = v;
            e_row[j] = max(sat_sub(e_row[j], P.e_del), sat_sub(v, P.oe_del));
            f = max(sat_sub(f, P.e_ins), sat_sub(t, P.oe_ins));
        }
        // Nachkorrektur: die offene Insertion von Spur zu Spur weiterreichen.
        var done = false;
        for (var round = 0u; round < 16u; round++) {
            let fu = subgroupShuffleUp(f, 1u);
            f = select(fu, 0u, k == 0u);
            for (var j = 0u; j < slen; j++) {
                let v = max(h1[j], f);
                h1[j] = v;
                f = sat_sub(f, P.e_ins);
                if (!subgroupAny(f > sat_sub(v, P.oe_ins))) {
                    done = true;
                    break;
                }
            }
            if (done) {
                break;
            }
        }
        let imax = subgroupMax(row_max);
        if (sl == 0u) {
            let ri = job.rows_off + i;
            atomicOr(&rows[ri >> 2u], imax << ((ri & 3u) * 8u));
        }
        if (imax > gmax) {
            gmax = imax;
            te = i32(i);
            for (var j = 0u; j < slen; j++) {
                h_max[j] = h1[j];
            }
        }
        for (var j = 0u; j < slen; j++) {
            h0[j] = h1[j];
        }
    }
    // Read-Position des besten Werts in der besten Zeile: größter Wert, bei Gleichstand die
    // kleinste Position. Schlüssel Wert * 65536 + (65535 - Position), dann Maximum.
    var key = 0u;
    for (var j = 0u; j < slen; j++) {
        let p = j + k * slen;
        key = max(key, (h_max[j] << 16u) | (65535u - p));
    }
    key = subgroupMax(key);
    if (sl == 0u) {
        result[id * 4u] = gmax;
        result[id * 4u + 1u] = bitcast<u32>(te);
        result[id * 4u + 2u] = bitcast<u32>(i32(65535u - (key & 0xffffu)));
    }
}
