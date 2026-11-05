# Stage 6 implementation checkpoint

Branch feature/time_step_resolution, worktree /private/tmp/mc3d-resolution, base d5421ed.

Implemented so far (not ready for merge): collision_courant config/default 0.1, resolution_diagnostics config/default true, effective config fields, RMS-relative-speed bound Cell::CollisionFrequencyBound, independent collision/transport/end-time minimum, reservoir transport cap for open boundaries. Removed numeric transport sentinel. No-progress check is in TimeStep, because initialization permits an infinite end time. Updated existing parallel-preparation oracle for new limit.

Build passed before latest small fixes. Initial tests exposed initialization with infinite end time; validation moved to execution. Latest changes need rebuild/retest.

Remaining required work:
- Dedicated analytic/direct-pair bound tests, low populations, omega=1, invariance, gas fractions, independent global step oracle, mutation tests.
- Cache per-cell diagnostics during parallel preparation; stable peculiar-speed sum, avoid redundant passes.
- Diagnostics files (summary each step, cells first/last/snapshot steps), equilibrium mean-free-path estimates, applicability statuses, all three axes and histogram bins per approved plan.
- Wire diagnostics through runner, SI generator, CLI help, manifest method/schema/output hashes and replay. Currently resolution_diagnostics is parsed but unused.
- New configs and ten-seed two-beam comparison for collision_courant 0.1/0.05, make run LOAD_MODE=balanced; comparison tolerance max(5 standard errors, 2% initial second moment).
- Full CTest/Python, repeated new tests, sanitizer CI; documentation and PHYSICS_PLAN.md status, commit/push and GitHub PR. Do not merge or start stage 7.

Approved definitions: N>=2, S2=sum|vi-u|², g_rms=sqrt(2S2/(N-1)), nu_bound=(N-1)Rate(g_rms)/(sqrt(2)Lref Kn N0 Vgas/Vcell). This bounds mean per-particle frequency, not every particle's rate. N<2 and none => zero. nu_eq=N E[Rate(g)]/(sqrt(2)Lref Kn N0 f), lambda_eq=sqrt(8T/pi)/nu_eq; unavailable for N<2,T<=0,solid,none.

Diagnostics summary bins: N 0,1,2-4,5-9,10-19,20-49,>=50; dx/lambda <.5,.5-1,1-2,2-5,>=5; dt*nu <=.05,.05-.1,.1-.2,>.2. Pre-step state with actual applied dt and explicit timestamps. No RNG consumption.
