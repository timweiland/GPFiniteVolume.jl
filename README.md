# GPFiniteVolume.jl

[![Build Status](https://github.com/timweiland/GPFiniteVolume.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/timweiland/GPFiniteVolume.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![arXiv](https://img.shields.io/badge/arXiv-2605.31127-b31b1b.svg)](https://arxiv.org/abs/2605.31127)

Scalable Bayesian inference for nonlinear conservation laws via Gaussian process finite volume methods (GP-FVM) with sparse Cholesky approximation.

This is the reference implementation for:

> **Scalable Bayesian Inference for Nonlinear Conservation Laws**  
> Tim Weiland and Philipp Hennig.  
> *International Conference on Machine Learning (ICML)*, 2026.  
> [arXiv:2605.31127](https://arxiv.org/abs/2605.31127)

## What this package does

GPFiniteVolume.jl casts the finite volume method as Bayesian inference under a structured Gaussian process prior. The posterior provides not just a numerical solution to a PDE but a calibrated uncertainty estimate over solution, fluxes, and (in inverse problems) unknown source fields. The key ideas are:

- **Mixed linear functionals.** Cell integrals, point evaluations, and derivatives are all linear functionals of the GP and can be conditioned jointly.
- **The Finite Volume Method is inference on a joint state over linear functionals.** See the paper for details.
- **Sparse Cholesky via Vecchia approximation.** A KL-minimizing sparse approximation of the prior precision yields `O(N_s^{3/2})` factorization in 2D.
- **Marginal moment-matching in time.** A linear-in-time-steps filter that preserves spatial sparsity through moment-matching, avoiding the dense fill-in of standard Kalman prediction.

## Installation

GPFiniteVolume.jl depends on [GaussianMarkovRandomFields.jl](https://github.com/timweiland/GaussianMarkovRandomFields.jl) (registered) and [FunctionalGPs.jl](https://github.com/timweiland/FunctionalGPs.jl) (unregistered).  Install via:

```julia
using Pkg
Pkg.add(url="https://github.com/timweiland/FunctionalGPs.jl")
Pkg.add(url="https://github.com/timweiland/GPFiniteVolume.jl")
```

Or, to develop locally and reproduce the paper figures, clone this repository and instantiate:

```bash
git clone https://github.com/timweiland/GPFiniteVolume.jl
cd GPFiniteVolume.jl
julia --project=. -e 'using Pkg; Pkg.add(url="https://github.com/timweiland/FunctionalGPs.jl"); Pkg.instantiate()'
```

## Minimal API example

```julia
using GPFiniteVolume
using FunctionalGPs, GaussianMarkovRandomFields

# Matérn 5/2 kernel
k = HalfIntegerMaternKernel(2, [0.15])
endpoints = collect(range(0, 1, length=51))
intervals = intervals_from_endpoints(endpoints)

# Build a sparse GMRF over function values, derivatives, and cell integrals
x0, layout, approx = sparse_gmrf([
    :f       => EvaluationFunctional(endpoints),
    :f_dx    => EvaluationFunctional(endpoints) ∘ PartialDerivative((1,)),
    :f_int   => VectorizedLebesgueIntegral(intervals),
], k; ρ=2.0, ordering=:integrals_coarsest)

# Condition on a subset of function-value observations
obs_indices = indices(layout, :f)[1:10:end]
x_cond = prescribe_indices(x0, obs_indices, observations; noise_std=1e-3)
```

See `experiments/` for full forward and inverse problem examples.

## Reproducing the paper figures

All paper figures can be regenerated from the scripts under `experiments/`. Default parameters reproduce the camera-ready figures. Most scripts use ArgParse and accept `--help`.

| Figure | Script | Notes |
|--------|--------|-------|
| Fig. 1 — Source identification overview (§1) | `experiments/source_identification/figure_upstream_downstream.jl` | Requires running `run_gpfvm.jl` on the relevant `problems/*.toml` configs first |
| Fig. 2 — Posterior adaptation under added observations (§4.1) | `experiments/source_identification/figure_posterior_adaptation.jl` | |
| Fig. 3 — Burgers benchmark accuracy + runtime (§4.2) | `experiments/accuracy_vs_compute/run.jl` then `paper_plots.jl` | |
| Fig. 4 — Poisson kernel-smoothness convergence (§4.3) | `experiments/poisson_convergence.jl` | |
| Fig. 5 — Nonlinear shallow water simulation (§4.4) | `experiments/nonlinear_shallow_water/run_showcase.jl` | Followed by `plot_results.jl` |
| Fig. 6 — Ordering comparison Pareto frontier (App. A.3) | `experiments/sparsity_accuracy_experiment.jl` | |
| Fig. 7 — Exact precision Cholesky factor (App. A.3) | `experiments/sparsity_pattern_2d.jl` | |
| Figs. 8, 9 — Calibration plot + spatial coverage map (App. C) | `experiments/source_identification/calibration_experiment.jl` | 200 random instances; long-running (hours) |
| Fig. 10 — Scalability timing (App. D) | `experiments/source_identification/scalability_study.jl` | Sweeps N = 11..101 |
| Figs. 11, 12 — Nonlinear Burgers source identification (App. E) | `experiments/burgers_source_identification/run.jl --log-source` | |
| Fig. 13 — Gauss--Newton ablation (App. F) | `experiments/accuracy_vs_compute/gn_ablation.jl` | |

Typical invocation:

```bash
julia --project=. experiments/source_identification/run_gpfvm.jl --help
julia --project=. experiments/burgers_source_identification/run.jl -N 16 --n-timesteps 11 --log-source
```

The calibration and scalability studies are long-running (hours on a single CPU). Other figures complete in seconds to minutes.

## Repository layout

```
src/                            Package source
test/                           Unit tests (julia --project=. -e 'using Pkg; Pkg.test()')
experiments/                    Paper-reproducing experiments
  source_identification/        §4.1 + Apps C, D
  burgers_source_identification/  App. E (nonlinear inverse problem)
  accuracy_vs_compute/          §4.2 + App. F
  nonlinear_shallow_water/      §4.3 shallow water demo
  poisson_convergence.jl        §4.3 + App. B.3
  sparsity_accuracy_experiment.jl  App. A.3 ordering comparison
  sparse_burgers_*.jl           §1, §4.2 forward Burgers demos
```

## Citation

If you use this work, please cite the paper. The arXiv preprint is available at [arXiv:2605.31127](https://arxiv.org/abs/2605.31127); this entry will be updated with the PMLR reference once the ICML 2026 proceedings are published.

```bibtex
@misc{weiland2026scalable,
  title         = {Scalable Bayesian Inference for Nonlinear Conservation Laws},
  author        = {Weiland, Tim and Hennig, Philipp},
  year          = {2026},
  eprint        = {2605.31127},
  archivePrefix = {arXiv},
  primaryClass  = {cs.LG},
}
```

## License

MIT. See [LICENSE](LICENSE).
