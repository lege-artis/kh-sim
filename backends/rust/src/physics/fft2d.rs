/// fft2d.rs — 2-D DFT/IDFT wrappers built on rustfft.
///
/// Convention matches NumPy rfft2 / irfft2 with indexing='ij':
///   - axis-0 (length nx): full complex FFT  (fftfreq)
///   - axis-1 (length ny): full complex FFT  (fftfreq)
///
/// We use full complex FFT in both axes (vs NumPy's rfft2 which exploits
/// conjugate symmetry in axis-1). The extra memory cost is negligible for the
/// grid sizes used here (up to 256x128).
///
/// Output layout: row-major (nx, ny), index = i*ny + j
use num_complex::Complex64;
use rustfft::{Fft, FftPlanner};
use std::sync::Arc;

/// Cached FFT plans for a fixed (nx, ny) grid.
pub struct Fft2DPlans {
    nx: usize,
    ny: usize,
    fwd_x: Arc<dyn Fft<f64>>,
    fwd_y: Arc<dyn Fft<f64>>,
    inv_x: Arc<dyn Fft<f64>>,
    inv_y: Arc<dyn Fft<f64>>,
}

impl Fft2DPlans {
    pub fn new(nx: usize, ny: usize) -> Self {
        let mut planner = FftPlanner::<f64>::new();
        Self {
            nx,
            ny,
            fwd_x: planner.plan_fft_forward(nx),
            fwd_y: planner.plan_fft_forward(ny),
            inv_x: planner.plan_fft_inverse(nx),
            inv_y: planner.plan_fft_inverse(ny),
        }
    }

    /// Forward 2-D DFT.
    /// Input:  real slice, length nx*ny, row-major.
    /// Output: complex Vec, length nx*ny, row-major.
    pub fn fft2(&self, real: &[f64]) -> Vec<Complex64> {
        let nx = self.nx;
        let ny = self.ny;
        debug_assert_eq!(real.len(), nx * ny);

        // Promote to complex
        let mut buf: Vec<Complex64> = real.iter().map(|&r| Complex64::new(r, 0.0)).collect();

        // FFT along axis-1 (rows of length ny)
        for i in 0..nx {
            let row = &mut buf[i * ny..(i + 1) * ny];
            self.fwd_y.process(row);
        }

        // FFT along axis-0 (columns of length nx) — gather into temp, scatter back
        let mut col = vec![Complex64::new(0.0, 0.0); nx];
        for j in 0..ny {
            for i in 0..nx {
                col[i] = buf[i * ny + j];
            }
            self.fwd_x.process(&mut col);
            for i in 0..nx {
                buf[i * ny + j] = col[i];
            }
        }

        buf
    }

    /// Inverse 2-D DFT with 1/(nx*ny) normalisation.
    /// Input:  complex slice, length nx*ny, row-major.
    /// Output: real Vec (imaginary parts are discarded; they should be ~0).
    pub fn ifft2_real(&self, spec: &[Complex64]) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        debug_assert_eq!(spec.len(), nx * ny);

        let mut buf = spec.to_vec();

        // IFFT along axis-0 (columns)
        let mut col = vec![Complex64::new(0.0, 0.0); nx];
        for j in 0..ny {
            for i in 0..nx {
                col[i] = buf[i * ny + j];
            }
            self.inv_x.process(&mut col);
            for i in 0..nx {
                buf[i * ny + j] = col[i];
            }
        }

        // IFFT along axis-1 (rows)
        for i in 0..nx {
            let row = &mut buf[i * ny..(i + 1) * ny];
            self.inv_y.process(row);
        }

        // Normalise and return real part
        let norm = (nx * ny) as f64;
        buf.iter().map(|c| c.re / norm).collect()
    }
}

// ── Wavenumber arrays ──────────────────────────────────────────────────────────

/// NumPy fftfreq(n, d) scaled by 2*PI.
/// Returns angular frequencies in rad/unit.
#[allow(clippy::needless_range_loop)]
pub fn angular_fftfreq(n: usize, d: f64) -> Vec<f64> {
    use std::f64::consts::PI;
    let mut out = vec![0.0f64; n];
    // positive half: index `i` used both to write out[i] and compute `i as f64`
    for i in 0..=(n / 2) {
        out[i] = 2.0 * PI * i as f64 / (n as f64 * d);
    }
    // negative half: index `i` used both to write out[i] and compute `i as f64 - n as f64`
    // (mirrors Python: freq[i] = (i - n) / (n*d) for i > n//2)
    for i in (n / 2 + 1)..n {
        out[i] = 2.0 * PI * (i as f64 - n as f64) / (n as f64 * d);
    }
    out
}

/// Flat (nx*ny) wavenumber arrays KX, KY (row-major, meshgrid indexing='ij').
pub fn make_wavenumbers(
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
) -> (Vec<f64>, Vec<f64>) {
    let kx_1d = angular_fftfreq(nx, dx);
    let ky_1d = angular_fftfreq(ny, dy);

    let mut kx = vec![0.0f64; nx * ny];
    let mut ky = vec![0.0f64; nx * ny];

    for i in 0..nx {
        for j in 0..ny {
            kx[i * ny + j] = kx_1d[i];
            ky[i * ny + j] = ky_1d[j];
        }
    }
    (kx, ky)
}
