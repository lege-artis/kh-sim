package khsim

/** CanonicalPhysicsTest.scala -- phase 7f canonical-physics testbeds, Scala backend.
  *
  * Every CP-7F testbed now runs once, for the oracle and all five backends, from a
  * shared harness in tests/integration (ruling 2, 2026-10-07):
  *   CP-7F-A  test_canonical_physics_cp7fa.py   (MASTER-BRIEF s3.1, E9; moved 2026-10-12)
  *   CP-7F-B  test_canonical_physics_cp7fbd.py  (s3.2, E7; moved 2026-10-10)
  *   CP-7F-C  test_canonical_physics_cp7fc.py   (s3.3, E6; moved 2026-10-07)
  *   CP-7F-D  test_canonical_physics_cp7fbd.py  (s3.4, E7; moved 2026-10-10)
  *   CP-7F-E  the Pascal-binding anchor (s3.5, E8): the oracle, Pascal and Fortran.
  * The Scala CP-7F-A that stood here (12 seeded modes at 256^2, gate [0.180, 0.200],
  * cancelled outside native mode) ran on the double layer at delta = Ly/20, where the
  * two layers couple (assessment 2026-09-15, s7); its text is in the git history.
  * The two tests below only state where the testbeds run; they need no server.
  */

import org.scalatest.funsuite.AnyFunSuite

class CanonicalPhysicsTest extends AnyFunSuite:

  test("CP-7F-A moved to the shared harness (MASTER s3.1 E9)") {
    info("CP-7F-A: SKIP -- moved to tests/integration/test_canonical_physics_cp7fa.py (MASTER s3.1 E9)")
  }

  test("CP-7F-E is the Pascal-binding anchor (MASTER s3.5 E8)") {
    info("CP-7F-E: SKIP -- Pascal-binding anchor (MASTER s3.5 E8): the oracle, Pascal and Fortran assert it")
  }
