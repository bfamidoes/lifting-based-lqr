# Lifting-Based LQR Design for Discrete-Time Periodic Systems: Attitude Control and Formation-Flying

This repository includes a minimalistic implementation of a solution to the periodic LQR problem via a time-lifting formulation, with subsequent computation of the periodic gains for stepwise (i.e., memoryless, or static) feedback. The performance of lifted feedback is compared with stepwise feedback in two examples: satellite attitude control and satellite formation flying.

## Requirements

MATLAB with the Control System Toolbox (for `idare`).

## Examples

Run either script. Each calls `setup.m`, which adds the folders to the path.

| Script | Example |
|---|---|
| `runAttitudeExample` | Correction of initial satellite attitude error with magnetic actuation. |
| `runFormationExample` | Formation reconfiguration (or station-keeping) in a Tundra orbit. |

## Structure

| Folder | Content |
|---|---|
| `lqrDesign/` | Main controller synthesis algorithms. `computePeriodicLQR` is the entry point. |
| `systems/` | Periodic system matrices of the two examples. |
| `simulation/` | Closed loops with either feedback law: linear periodic model, and nonlinear two-body dynamics for formation flying. |
| `simulation/orbitalMechanics/` | Auxiliary orbital mechanics functions (orbital elements, frame/time conversions). |
| `plotting/` | Postprocessing functions. |

## License

[MIT](LICENSE)
