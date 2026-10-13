package khsim.physics

import org.jtransforms.fft.DoubleFFT_2D

/** 2-D DFT / IDFT backed by JTransforms 3.1 DoubleFFT_2D.
  *
  * Layout: row-major (nx, ny), flat index = i*ny + j.
  * JTransforms interleaved complex: buf(2k) = real, buf(2k+1) = imag.
  *
  * Convention matches NumPy fft2 / ifft2 with indexing='ij':
  *   axis-0 (length nx): full complex FFT  (fftfreq)
  *   axis-1 (length ny): full complex FFT  (fftfreq)
  */
final class Fft2D(val nx: Int, val ny: Int):
  // JTransforms DoubleFFT_2D(rows, cols) — nx rows, ny columns
  private val engine = new DoubleFFT_2D(nx.toLong, ny.toLong)

  /** Forward 2-D DFT.
    * @param real row-major real array of length nx*ny
    * @return interleaved complex array of length 2*nx*ny
    */
  def fft2(real: Array[Double]): Array[Double] =
    require(real.length == nx * ny)
    val buf = new Array[Double](2 * nx * ny)
    var k = 0
    while k < nx * ny do
      buf(2 * k) = real(k)
      // buf(2*k+1) = 0.0 already
      k += 1
    engine.complexForward(buf)
    buf

  /** Inverse 2-D DFT with 1/(nx*ny) normalisation; returns real part.
    * @param complex interleaved array of length 2*nx*ny
    * @return real array of length nx*ny
    */
  def ifft2Real(complex: Array[Double]): Array[Double] =
    require(complex.length == 2 * nx * ny)
    val buf = complex.clone()
    engine.complexInverse(buf, /*scale=*/true)
    val out = new Array[Double](nx * ny)
    var k = 0
    while k < nx * ny do
      out(k) = buf(2 * k)
      k += 1
    out

// ── Wavenumber arrays ──────────────────────────────────────────────────────────

object Wavenumbers:
  /** Angular fftfreq(n, d) = 2π * k / (n*d) for k in [0..n/2, -(n/2-1)..-1] */
  def angularFftfreq(n: Int, d: Double): Array[Double] =
    val out = new Array[Double](n)
    val half = n / 2
    var i = 0
    while i <= half do
      out(i) = 2.0 * math.Pi * i / (n * d)
      i += 1
    i = half + 1
    while i < n do
      out(i) = 2.0 * math.Pi * (i - n) / (n * d)
      i += 1
    out

  /** Flat (nx*ny) KX, KY arrays in row-major meshgrid indexing='ij'. */
  def make(nx: Int, ny: Int, dx: Double, dy: Double): (Array[Double], Array[Double]) =
    val kx1d = angularFftfreq(nx, dx)
    val ky1d = angularFftfreq(ny, dy)
    val kx   = new Array[Double](nx * ny)
    val ky   = new Array[Double](nx * ny)
    var i = 0
    while i < nx do
      var j = 0
      while j < ny do
        val idx = i * ny + j
        kx(idx) = kx1d(i)
        ky(idx) = ky1d(j)
        j += 1
      i += 1
    (kx, ky)
