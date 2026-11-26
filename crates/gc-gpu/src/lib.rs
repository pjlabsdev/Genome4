//! Rechnen auf der Grafikeinheit über wgpu (Metal auf dem Mac, Vulkan oder DX12 sonst).
//!
//! Bisher ein Rechenkern: die lokalen Ausrichtungen der Partnersuche (gestreifter Byte-Lauf wie
//! `gc_align::local`, eine Ausrichtung je Aufruf auf der Grafikeinheit). Er liefert je Fenster
//! den besten Score, seine Zeile, die Read-Position und das Maximum jeder Zeile; daraus rechnet
//! `gc_align::local::from_rows` dasselbe Ergebnis wie der Prozessor. Alles bleibt auf dem Gerät:
//! Die Grafikeinheit sitzt im selben Rechner und teilt sich beim Mac den Arbeitsspeicher.

use bytemuck::{Pod, Zeroable};
use wgpu::util::DeviceExt;

/// Bewertung, wie sie der Rechenkern braucht (Abzüge positiv).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KernelScoring {
    pub matching: u32,
    pub mismatch: u32,
    pub ambiguous: u32,
    pub o_del: u32,
    pub e_del: u32,
    pub o_ins: u32,
    pub e_ins: u32,
}

/// Ergebnis eines Fensters: bester Score, Zeile (-1 = keine), Read-Position, Zeilenmaxima.
pub type Forward = (i32, i32, i32, Vec<u8>);

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
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

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct Job {
    q_off: u32,
    q_len: u32,
    t_off: u32,
    t_len: u32,
    rows_off: u32,
    pad: [u32; 3],
}

/// Längster Read, den der Rechenkern kann (Felder fester Größe auf der Grafikeinheit).
pub const MAX_QUERY: usize = 256;

/// Die Grafikeinheit mit dem übersetzten Rechenkern.
pub struct Gpu {
    device: wgpu::Device,
    queue: wgpu::Queue,
    pipeline: wgpu::ComputePipeline,
    /// Name und Schnittstelle der Grafikeinheit, für das Protokoll.
    pub name: String,
}

impl std::fmt::Debug for Gpu {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Gpu({})", self.name)
    }
}

/// Der Quelltext des Rechenkerns.
pub const RESCUE_WGSL: &str = include_str!("rescue.wgsl");

impl Gpu {
    /// Öffnet die schnellste Grafikeinheit. Fehler als Text, wenn es keine gibt.
    pub fn open() -> Result<Self, String> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            force_fallback_adapter: false,
            compatible_surface: None,
            ..Default::default()
        }))
        .map_err(|e| format!("keine Grafikeinheit: {e}"))?;
        let info = adapter.get_info();
        let name = format!("{} ({:?})", info.name, info.backend);
        let limits = adapter.limits();
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            label: Some("genome"),
            required_limits: limits,
            ..Default::default()
        }))
        .map_err(|e| format!("Grafikeinheit {name} nicht nutzbar: {e}"))?;
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("partnersuche"),
            source: wgpu::ShaderSource::Wgsl(RESCUE_WGSL.into()),
        });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("partnersuche"),
            layout: None,
            module: &module,
            entry_point: Some("main"),
            compilation_options: Default::default(),
            cache: None,
        });
        Ok(Self {
            device,
            queue,
            pipeline,
            name,
        })
    }

    /// Richtet jedes Paar (Read, Referenzstück) aus. Codes 0..3, 4 = N. Jeder Read höchstens
    /// [`MAX_QUERY`] Basen und nicht leer.
    pub fn forward(
        &self,
        jobs: &[(&[u8], &[u8])],
        sc: &KernelScoring,
    ) -> Result<Vec<Forward>, String> {
        if jobs.is_empty() {
            return Ok(Vec::new());
        }
        let mut bytes: Vec<u8> = Vec::new();
        let mut descs = Vec::with_capacity(jobs.len());
        let mut rows_len = 0usize;
        for (q, t) in jobs {
            assert!(
                !q.is_empty() && q.len() <= MAX_QUERY,
                "Read zu lang für die Grafikeinheit"
            );
            let q_off = bytes.len() as u32;
            bytes.extend_from_slice(q);
            let t_off = bytes.len() as u32;
            bytes.extend_from_slice(t);
            descs.push(Job {
                q_off,
                q_len: q.len() as u32,
                t_off,
                t_len: t.len() as u32,
                rows_off: rows_len as u32,
                pad: [0; 3],
            });
            rows_len += t.len();
        }
        bytes.resize(bytes.len().next_multiple_of(4).max(4), 0);
        let rows_bytes = rows_len.next_multiple_of(4).max(4);
        let params = Params {
            n_jobs: jobs.len() as u32,
            matching: sc.matching,
            mismatch: sc.mismatch,
            ambiguous: sc.ambiguous,
            e_del: sc.e_del,
            oe_del: sc.o_del + sc.e_del,
            e_ins: sc.e_ins,
            oe_ins: sc.o_ins + sc.e_ins,
        };
        let d = &self.device;
        let storage = wgpu::BufferUsages::STORAGE;
        let p_buf = d.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: None,
            contents: bytemuck::bytes_of(&params),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let j_buf = d.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: None,
            contents: bytemuck::cast_slice(&descs),
            usage: storage,
        });
        let b_buf = d.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: None,
            contents: &bytes,
            usage: storage,
        });
        let res_size = (jobs.len() * 16) as u64;
        let r_buf = d.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size: res_size,
            usage: storage | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        // Neue Puffer sind mit Nullen gefüllt; der Kern setzt die Zeilenmaxima mit Oder.
        let rows_buf = d.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size: rows_bytes as u64,
            usage: storage | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let read_buf = d.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size: res_size + rows_bytes as u64,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let layout = self.pipeline.get_bind_group_layout(0);
        let bind = d.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: p_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: j_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: b_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: r_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 4,
                    resource: rows_buf.as_entire_binding(),
                },
            ],
        });
        let mut enc = d.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        {
            let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: None,
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.pipeline);
            pass.set_bind_group(0, &bind, &[]);
            let groups = (jobs.len() as u32).div_ceil(64);
            // Höchstens 65.535 Gruppen je Richtung: darüber in die zweite Richtung ausweichen
            // ist nicht nötig, solange ein Block unter 4 Mio. Fenstern bleibt.
            assert!(groups <= 65_535, "zu viele Fenster für einen Aufruf");
            pass.dispatch_workgroups(groups, 1, 1);
        }
        enc.copy_buffer_to_buffer(&r_buf, 0, &read_buf, 0, res_size);
        enc.copy_buffer_to_buffer(&rows_buf, 0, &read_buf, res_size, rows_bytes as u64);
        self.queue.submit([enc.finish()]);
        let slice = read_buf.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r);
        });
        d.poll(wgpu::PollType::wait_indefinitely())
            .map_err(|e| format!("Grafikeinheit: {e}"))?;
        rx.recv()
            .map_err(|e| format!("Grafikeinheit: {e}"))?
            .map_err(|e| format!("Grafikeinheit: {e}"))?;
        let data = slice
            .get_mapped_range()
            .map_err(|e| format!("Grafikeinheit: {e}"))?;
        let res: &[u32] = bytemuck::cast_slice(&data[..res_size as usize]);
        let rows = &data[res_size as usize..];
        let out = descs
            .iter()
            .enumerate()
            .map(|(i, j)| {
                let r = j.rows_off as usize;
                (
                    res[4 * i] as i32,
                    res[4 * i + 1] as i32,
                    res[4 * i + 2] as i32,
                    rows[r..r + j.t_len as usize].to_vec(),
                )
            })
            .collect();
        drop(data);
        read_buf.unmap();
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_kernel_source_is_valid_wgsl() {
        let module = naga::front::wgsl::parse_str(RESCUE_WGSL).expect("WGSL lesbar");
        naga::valid::Validator::new(
            naga::valid::ValidationFlags::all(),
            naga::valid::Capabilities::all(),
        )
        .validate(&module)
        .expect("WGSL gültig");
    }
}
