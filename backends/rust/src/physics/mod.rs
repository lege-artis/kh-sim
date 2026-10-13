/// physics — KH instability physics module
///
/// Public re-exports:
///   - simulate(req) -> SimulationResult   (main entry point used by handlers)
pub mod fft2d;
pub mod solver;

pub use solver::simulate;
