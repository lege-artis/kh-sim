# lege-artis/kh-sim — v0.1.4

**Tag:** `v0.1.4`
**Type:** contract change in all five HTTP backends: an unknown key in a request is now refused.

## Unknown keys are refused (KH-VAL-2 K7, step 2)

Previously, a key the request contract does not define was logged and ignored. A request that sent `nsteps` instead of `steps` therefore ran with the default number of steps and returned a plausible but wrong answer.

The C++, Fortran, Rust, Scala and Pascal backends now refuse such a request:

| request | result |
|---|---|
| a key not in `SimulationRequest` | **400**, `field` = that key |
| several unknown keys | 400 naming the first in byte order — the same answer from every backend |
| an unknown key together with an invalid value | 400 naming the unknown key: the parser checks keys before any value |
| a foreign name with a contract equivalent | the message adds a hint, e.g. `unknown key 'nsteps' (did you mean 'steps'?)` |

Hinted names:

| foreign name | hint |
|---|---|
| `nsteps`, `n_steps`, `num_steps` | `steps` |
| `nx`, `ny` | `grid_nx`, `grid_ny` |
| `lx`, `ly` | `domain_lx`, `domain_ly` |
| `re`, `Re`, `reynolds` | `reynolds_number` |
| `amplitude` | `perturbation_amplitude` |
| `mode` | `perturbation_mode` |

The legacy `init_omega` remains part of the contract.

`shared/api/openapi.yaml` (version 2.2.0) declares the request schema closed (`additionalProperties: false`). The `test_conformance.py` cases that pin K7 are:
- the refusal, with and without a hint;
- the precedence: an unknown key is named before an invalid value.

**Before the switch.** The backends' logs, covering each container's life, recorded every unknown key received. They held only the `nsteps` that the conformance suite sent on purpose. The harnesses of the test suite already translate the golden-vector `nsteps` to the contract's `steps`.

## Other

`backends/cpp/CMakeLists.txt`: the comment on the canonical-physics executable now says where the CP-7F testbeds run (the shared harnesses).
