package khsim.physics

import khsim.{Diagnostics, SimulationRequest, SimulationResult}

/** KH instability physics kernel — Scala port of kh_physics.py.
  *
  * Physics: 2D incompressible Navier-Stokes, vorticity-streamfunction form.
  * Numerics: pseudo-spectral (JTransforms 2D FFT) + RK4 time integration.
  *
  * Reference: kh-sim/shared/physics/KH-PHYSICS.md
  * Canonical: kh-sim/shared/physics/kh_physics.py
  */
object Solver:

  // ── Initial conditions ─────────────────────────────────────────────────────

  private def initialConditions(
    nx: Int, ny: Int, lx: Double, ly: Double,
    u0: Double, delta: Double, amp: Double, mode: Int,
  ): Array[Double] =
    val dx    = lx / nx
    val dy    = ly / ny
    val omega = new Array[Double](nx * ny)
    var i = 0
    while i < nx do
      val x = i * dx
      val pertX = amp * (2.0 * math.Pi * mode / lx) * math.cos(2.0 * math.Pi * mode * x / lx)
      var j = 0
      while j < ny do
        val y     = j * dy
        val z     = (y - ly / 2.0) / delta
        val shear = u0 / delta / (math.cosh(z) * math.cosh(z))
        omega(i * ny + j) = shear + pertX
        j += 1
      i += 1
    omega

  // ── Spectral Poisson solve ─────────────────────────────────────────────────

  /** psi_hat = omega_hat / k^2 ; psi_hat[0] = 0 (zero mean). */
  private def solvePoisson(
    omegaHat: Array[Double],  // interleaved complex, length 2*nx*ny
    kx: Array[Double],
    ky: Array[Double],
    nx: Int, ny: Int,
  ): Array[Double] =
    val psiHat = omegaHat.clone()
    var i = 0
    while i < nx do
      var j = 0
      while j < ny do
        val idx = i * ny + j
        val k2  = kx(idx) * kx(idx) + ky(idx) * ky(idx)
        if i == 0 && j == 0 then
          psiHat(0) = 0.0
          psiHat(1) = 0.0
        else
          val re = omegaHat(2 * idx)
          val im = omegaHat(2 * idx + 1)
          psiHat(2 * idx)     = re / k2
          psiHat(2 * idx + 1) = im / k2
        j += 1
      i += 1
    psiHat

  // ── Velocity from streamfunction ───────────────────────────────────────────

  /** u_hat = i*ky*psi_hat,  v_hat = -i*kx*psi_hat */
  private def velocityFromPsiHat(
    psiHat: Array[Double],
    kx: Array[Double],
    ky: Array[Double],
    fft: Fft2D,
  ): (Array[Double], Array[Double]) =
    val n     = fft.nx * fft.ny
    val uHat  = new Array[Double](2 * n)
    val vHat  = new Array[Double](2 * n)
    var k = 0
    while k < n do
      val re = psiHat(2 * k)
      val im = psiHat(2 * k + 1)
      // i * ky * (re + i*im) = ky*(-im) + i*ky*re
      uHat(2 * k)     = -ky(k) * im
      uHat(2 * k + 1) =  ky(k) * re
      // -i * kx * (re + i*im) = kx*im + i*(-kx*re)
      vHat(2 * k)     =  kx(k) * im
      vHat(2 * k + 1) = -kx(k) * re
      k += 1
    (fft.ifft2Real(uHat), fft.ifft2Real(vHat))

  // ── Vorticity RHS ──────────────────────────────────────────────────────────

  private def vorticityRhs(
    omega: Array[Double],
    u: Array[Double],
    v: Array[Double],
    nu: Double,
    kx: Array[Double],
    ky: Array[Double],
    fft: Fft2D,
    forcing: Array[Double],  // F(x,y,t_stage); null => no forcing (pre-7f path)
    drag: Double,            // Ekman linear drag alpha; 0.0 => no drag (pre-7f path)
  ): Array[Double] =
    val n        = fft.nx * fft.ny
    val omegaHat = fft.fft2(omega)

    // spectral gradients: domega/dx, domega/dy
    val dxSpec = new Array[Double](2 * n)
    val dySpec = new Array[Double](2 * n)
    val lapSpec = new Array[Double](2 * n)
    var k = 0
    while k < n do
      val re = omegaHat(2 * k)
      val im = omegaHat(2 * k + 1)
      // i*kx*(re+i*im) = -kx*im + i*kx*re
      dxSpec(2 * k)     = -kx(k) * im
      dxSpec(2 * k + 1) =  kx(k) * re
      // i*ky*(re+i*im)
      dySpec(2 * k)     = -ky(k) * im
      dySpec(2 * k + 1) =  ky(k) * re
      // -(kx^2+ky^2) * (re+i*im)
      val k2 = kx(k) * kx(k) + ky(k) * ky(k)
      lapSpec(2 * k)     = -k2 * re
      lapSpec(2 * k + 1) = -k2 * im
      k += 1

    val dOmegaDx = fft.ifft2Real(dxSpec)
    val dOmegaDy = fft.ifft2Real(dySpec)
    val lapOmega = fft.ifft2Real(lapSpec)

    // Nonlinear term with 2/3-rule dealiasing (KH-NM-4, Orszag 1971)
    val nlPhys = new Array[Double](n)
    var i = 0
    while i < n do
      nlPhys(i) = -u(i) * dOmegaDx(i) - v(i) * dOmegaDy(i)
      i += 1
    val nonlinear = fft.ifft2Real(dealiasTwoThirds(fft.fft2(nlPhys), kx, ky, n))

    val rhs = new Array[Double](n)
    i = 0
    while i < n do
      // domega/dt = -(u.grad)omega + nu*lap(omega) - alpha*omega + F(x,y,t)
      //            (PHASE-7F-FORCING-API section 2 + section 7)
      var r = nonlinear(i) + nu * lapOmega(i)
      if drag != 0.0 then r -= drag * omega(i)
      if forcing != null then r += forcing(i)
      rhs(i) = r
      i += 1
    rhs

  // ── 2/3-rule dealiasing + CFL bound ───────────────────────────────────────

  /** KH-NM-4: 2/3-rule dealiasing (Orszag 1971). Returns dealiased copy of interleaved spectrum. */
  private def dealiasTwoThirds(
    fieldHat: Array[Double],
    kx: Array[Double],
    ky: Array[Double],
    n: Int,
  ): Array[Double] =
    val kxMax = kx.map(_.abs).max
    val kyMax = ky.map(_.abs).max
    val kxCut = (2.0 / 3.0) * kxMax
    val kyCut = (2.0 / 3.0) * kyMax
    val out = fieldHat.clone()
    var k = 0
    while k < n do
      if kx(k).abs > kxCut || ky(k).abs > kyCut then
        out(2 * k)     = 0.0
        out(2 * k + 1) = 0.0
      k += 1
    out

  /** KH-NM-5: CFL-limited time step. Returns min of advective and diffusive bounds. */
  private def cflDtBound(uMax: Double, vMax: Double, dx: Double, dy: Double, nu: Double): Double =
    val dtAdv  = math.min(dx / (uMax + 1e-30), dy / (vMax + 1e-30))
    // nu = 0 (inviscid) is a valid limit: no diffusive bound (modelling rules v0.1, D).
    val dtDiff =
      if nu == 0.0 then Double.PositiveInfinity
      else math.min(dx * dx / (2.0 * nu), dy * dy / (2.0 * nu))
    math.min(dtAdv, dtDiff)

  // ── Helpers ────────────────────────────────────────────────────────────────

  private def addScaled(a: Array[Double], b: Array[Double], s: Double): Array[Double] =
    val out = new Array[Double](a.length)
    var i = 0
    while i < a.length do
      out(i) = a(i) + s * b(i)
      i += 1
    out

  // ── RK4 step ───────────────────────────────────────────────────────────────

  // -- Forcing field (phase-7f v0.2; PHASE-7F-FORCING-API-DESIGN-v0.1 s3+s7) --
  //
  // F is added to the vorticity RHS. Two constructions, RNG-free + bit-portable
  // across all six backends (deterministic phases from frac(ix*sqrt(p)+iy*sqrt(q))).
  //   "mode": F(x,y) = A*cos(2*pi*kx*x/lx)*cos(2*pi*ky*y/ly)  (strictly steady)
  //   "band": F(x,y,t) = A/sqrt(M) * sum_m cos(2*pi*(ix*x/lx + iy*y/ly) + phi0_m + Omega_m*t)
  //           phi0_m  = 2*pi*frac(ix*sqrt2 + iy*sqrt3)        (v0.1 static phase)
  //           Omega_m = 2*pi*frac(ix*sqrt11 + iy*sqrt13)*omega0   (v0.2 rotation)
  // Half-plane mode set: ix>=0, skip (ix==0 && iy<=0), kmin <= sqrt(ix^2+iy^2) <= kmax.
  // "band" requires lx == ly (integer-mode-radius isotropy).
  private object Forcing:
    private def frac(x: Double): Double = x - math.floor(x)

    final case class Band(
      cosField: Array[Array[Double]],
      sinField: Array[Array[Double]],
      omega:    Array[Double],
      ampScale: Double,
    )

    /** Build the static "mode" forcing field (time-independent). */
    def buildMode(nx: Int, ny: Int, lx: Double, ly: Double,
                  a: Double, kx: Int, ky: Int): Array[Double] =
      val f  = new Array[Double](nx * ny)
      val dx = lx / nx
      val dy = ly / ny
      var i = 0
      while i < nx do
        val cx = math.cos(2.0 * math.Pi * kx * (i * dx) / lx)
        var j = 0
        while j < ny do
          f(i * ny + j) = a * cx * math.cos(2.0 * math.Pi * ky * (j * dy) / ly)
          j += 1
        i += 1
      f

    /** Build the band-forcing mode set (per-mode precomputed fields). */
    def buildBand(nx: Int, ny: Int, lx: Double, ly: Double,
                  a: Double, kmin: Double, kmax: Double, omega0: Double): Band =
      require(math.abs(lx - ly) < 1e-12,
        "band forcing requires lx == ly (integer-mode-radius isotropy)")
      val dx = lx / nx
      val dy = ly / ny
      val s2  = math.sqrt(2.0);  val s3  = math.sqrt(3.0)
      val s11 = math.sqrt(11.0); val s13 = math.sqrt(13.0)
      val kmaxI = math.ceil(kmax).toInt
      val cosBuf = collection.mutable.ArrayBuffer[Array[Double]]()
      val sinBuf = collection.mutable.ArrayBuffer[Array[Double]]()
      val omBuf  = collection.mutable.ArrayBuffer[Double]()
      var ix = 0
      while ix <= kmaxI do
        var iy = -kmaxI
        while iy <= kmaxI do
          val krad = math.sqrt((ix * ix + iy * iy).toDouble)
          val keep = !(ix == 0 && iy <= 0) && krad >= kmin && krad <= kmax
          if keep then
            val phi0 = 2.0 * math.Pi * frac(ix * s2 + iy * s3)
            val om   = 2.0 * math.Pi * frac(ix * s11 + iy * s13) * omega0
            val cf = new Array[Double](nx * ny)
            val sf = new Array[Double](nx * ny)
            var i = 0
            while i < nx do
              val xterm = ix * (i * dx) / lx
              var j = 0
              while j < ny do
                val arg = 2.0 * math.Pi * (xterm + iy * (j * dy) / ly) + phi0
                cf(i * ny + j) = math.cos(arg)
                sf(i * ny + j) = math.sin(arg)
                j += 1
              i += 1
            cosBuf += cf; sinBuf += sf; omBuf += om
          iy += 1
        ix += 1
      val m = cosBuf.length
      val ampScale = if m > 0 then a / math.sqrt(m.toDouble) else 0.0
      Band(cosBuf.toArray, sinBuf.toArray, omBuf.toArray, ampScale)

    /** Evaluate the band forcing field at time t:
      * F(t) = (A/sqrt(M)) * sum_m [ C_m*cos(Omega_m t) - S_m*sin(Omega_m t) ]. */
    def evalBand(b: Band, t: Double, n: Int): Array[Double] =
      val out = new Array[Double](n)
      var m = 0
      while m < b.omega.length do
        val ct = math.cos(b.omega(m) * t)
        val st = math.sin(b.omega(m) * t)
        val cf = b.cosField(m); val sf = b.sinField(m)
        var i = 0
        while i < n do
          out(i) += b.ampScale * (cf(i) * ct - sf(i) * st)
          i += 1
        m += 1
      out
  end Forcing

  // forceAt: maps a stage time t -> F(x,y,t) field (or null for no forcing).
  // For "mode" / no-forcing the field is time-independent (forceAt ignores t).
  private def rk4Step(
    omega: Array[Double],
    nu: Double,
    kx: Array[Double],
    ky: Array[Double],
    fft: Fft2D,
    dt: Double,
    t: Double,
    forceAt: Double => Array[Double],  // stage-time -> forcing field (null allowed)
    drag: Double,
  ): Array[Double] =
    def f(w: Array[Double], ts: Double): Array[Double] =
      val wHat   = fft.fft2(w)
      val psiHat = solvePoisson(wHat, kx, ky, fft.nx, fft.ny)
      val (u, v) = velocityFromPsiHat(psiHat, kx, ky, fft)
      vorticityRhs(w, u, v, nu, kx, ky, fft, forceAt(ts), drag)

    // RK4 stage times: k1 @ t, k2/k3 @ t+dt/2, k4 @ t+dt (PHASE-7F-FORCING s7).
    val k1 = f(omega, t)
    val k2 = f(addScaled(omega, k1, 0.5 * dt), t + 0.5 * dt)
    val k3 = f(addScaled(omega, k2, 0.5 * dt), t + 0.5 * dt)
    val k4 = f(addScaled(omega, k3, dt), t + dt)

    val out = new Array[Double](omega.length)
    var i = 0
    while i < omega.length do
      out(i) = omega(i) + (dt / 6.0) * (k1(i) + 2*k2(i) + 2*k3(i) + k4(i))
      i += 1
    out

  // ── Diagnostics ────────────────────────────────────────────────────────────

  private def mean(a: Array[Double]): Double =
    var s = 0.0; var i = 0
    while i < a.length do { s += a(i); i += 1 }
    s / a.length

  private def meanSq(a: Array[Double]): Double =
    var s = 0.0; var i = 0
    while i < a.length do { s += a(i) * a(i); i += 1 }
    s / a.length

  private def maxAbs(a: Array[Double]): Double =
    var m = 0.0; var i = 0
    while i < a.length do { val v = math.abs(a(i)); if v > m then m = v; i += 1 }
    m

  private def computeDiagnostics(
    omega: Array[Double],
    u: Array[Double],
    v: Array[Double],
    kx: Array[Double],
    ky: Array[Double],
    fft: Fft2D,
    cflPeak: Double,
  ): Diagnostics =
    val n   = fft.nx * fft.ny
    val ke  = 0.5 * mean(Array.tabulate(n)(i => u(i)*u(i) + v(i)*v(i)))
    val ens = 0.5 * meanSq(omega)
    val mv  = maxAbs(omega)

    // Divergence via spectral differentiation
    val uHat = fft.fft2(u)
    val vHat = fft.fft2(v)
    val divSpec = new Array[Double](2 * n)
    var k = 0
    while k < n do
      val dxRe = -kx(k) * uHat(2*k+1)
      val dxIm =  kx(k) * uHat(2*k)
      val dyRe = -ky(k) * vHat(2*k+1)
      val dyIm =  ky(k) * vHat(2*k)
      divSpec(2*k)   = dxRe + dyRe
      divSpec(2*k+1) = dxIm + dyIm
      k += 1
    val div    = fft.ifft2Real(divSpec)
    val divRms = math.sqrt(meanSq(div))

    Diagnostics(ke, ens, mv, divRms, cflPeak)

  // ── Main entry point ───────────────────────────────────────────────────────

  def simulate(req: SimulationRequest): SimulationResult =
    val t0     = System.nanoTime()
    val nx     = req.nx
    val ny     = req.ny
    val lx     = req.lx
    val ly     = req.ly
    val dt     = req.dtVal
    val nSteps = req.nSteps
    // nu = |U| / Re: a reversed shear is the mirrored flow, not negative viscosity
    // (modelling rules v0.1, K3). Identical for U >= 0.
    val nu     = math.abs(req.u0) / req.re
    val delta  = 0.05 * ly

    val dx = lx / nx
    val dy = ly / ny

    val fft         = new Fft2D(nx, ny)
    val (kx, ky)    = Wavenumbers.make(nx, ny, dx, dy)

    var omega: Array[Double] = req.initialOmega match
      case Some(v) =>
        require(v.length == nx * ny, "initialOmega length must equal nx*ny")
        v.toArray
      case None =>
        initialConditions(nx, ny, lx, ly, req.u0, delta, req.amp, req.mode)

    // CFL guard (KH-NM-5): warn once if user-supplied dt exceeds stability bound
    val omegaHat0      = fft.fft2(omega)
    val psiHat0        = solvePoisson(omegaHat0, kx, ky, nx, ny)
    val (u0, v0)       = velocityFromPsiHat(psiHat0, kx, ky, fft)
    val uMax0          = u0.map(_.abs).max
    val vMax0          = v0.map(_.abs).max
    val dtBound        = cflDtBound(uMax0, vMax0, dx, dy, nu)
    if dt > dtBound then
      System.err.println(f"KH-NM-5: dt=$dt%.3e exceeds CFL bound $dtBound%.3e; simulation may be unstable")

    // Time loop — track peak CFL across all steps (KH-NM-5 diagnostic; mirrors kh_physics.py)
    // Forcing setup (phase-7f v0.2). "none" => forceAt returns null => zero-cost,
    // bit-identical to pre-7f RHS. "mode" => steady field; "band" => rotating (s7).
    val nGrid = nx * ny
    val forceAt: Double => Array[Double] = req.forceType match
      case "mode" =>
        val fField = Forcing.buildMode(nx, ny, lx, ly, req.forceAmp, req.forceKx, req.forceKy)
        (_: Double) => fField
      case "band" =>
        val band = Forcing.buildBand(nx, ny, lx, ly, req.forceAmp,
          req.forceKmin, req.forceKmax, req.forceOmega0)
        (ts: Double) => Forcing.evalBand(band, ts, nGrid)
      case _ =>
        (_: Double) => null
    val drag = req.drag

    var t       = req.tStartVal   // a continued run keeps the forcing phase (2026-10-10)
    var step    = 0
    var cflPeak = 0.0
    while step < nSteps do
      val wHatStep        = fft.fft2(omega)
      val psiHatStep      = solvePoisson(wHatStep, kx, ky, nx, ny)
      val (uStep, vStep)  = velocityFromPsiHat(psiHatStep, kx, ky, fft)
      val cflStep         = dt * (maxAbs(uStep) / dx + maxAbs(vStep) / dy)
      if cflStep > cflPeak then cflPeak = cflStep
      omega  = rk4Step(omega, nu, kx, ky, fft, dt, t, forceAt, drag)
      t     += dt
      step  += 1

    // Final field recovery
    val omegaHat    = fft.fft2(omega)
    val psiHat      = solvePoisson(omegaHat, kx, ky, nx, ny)
    val psi         = fft.ifft2Real(psiHat)
    val (u, v)      = velocityFromPsiHat(psiHat, kx, ky, fft)
    val diagnostics = computeDiagnostics(omega, u, v, kx, ky, fft, cflPeak)

    val elapsedMs = (System.nanoTime() - t0) / 1e6
    val tRounded  = math.round(t * 1e10) / 1e10

    SimulationResult(
      backend        = "scala-http4s",
      language       = "Scala",
      stepsCompleted = step,   // what the loop ran: 0 for a negative steps (was nSteps)
      tFinal         = tRounded,
      gridNx         = nx,
      gridNy         = ny,
      uVelocity      = u.toVector,
      vVelocity      = v.toVector,
      vorticity      = omega.toVector,
      pressure       = psi.toVector,
      diagnostics    = diagnostics,
      computeTimeMs  = math.round(elapsedMs * 100) / 100.0,
    )
