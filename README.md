# TBRe25 Telemetry Validation and Vehicle Balance

Sensor validation, tyre-model fitting, and steady-state balance analysis for the **TBRe25** Formula Student car (Team Bath Racing Electric), built on the FSG25 endurance run.

The work is one chain, not four projects: **sensor validation → tyre model → vehicle model → on-track confirmation.** Each stage checks the one before it, and the final stage is what validates the whole chain.

**Dataset:** Formula Student Germany 2025 endurance, 71,444 samples @ 20 Hz, AiM logger. Confirmed as the competition run from the team run record.

---

## 1. A 14% yaw-rate scale error, found in my own published work

The primary yaw-rate channel under-reads by **~14%**. Required gain **1.141**, established four independent ways:

| Method | Implied gain |
|---|---|
| Closed-loop heading integral (gyro turns vs GPS course turns over the run) | **1.141** |
| Total least squares, lateral accel vs `V·r`, steady-state only | 1.168 |
| TLS, yaw rate vs GPS course rate | 1.167 |
| Gain that minimises spurious sideslip accumulation | 1.15–1.17 |

The accelerometer and GPS are the references, and both are independently anchored: the accelerometer triad reads **1.0079 g at rest** against gravity (correct to 0.8%), and GPS speed validates to **0.13%** against position differentiation. In steady cornering `a_y = V·r` is exact. Two validated inputs disagree with the gyro, so the gyro is the outlier. Over the run it accumulates **14.88 turns against 16.97 from GPS**.

It is a pure scale error: stationary bias is 0.005 °/s, there is no clipping, and the ratio is flat across magnitude bins.

![Channel validation](figures/channel_validation.png)

*Left: measured lateral acceleration against the kinematic reference `V·r`, which involves no accelerometer. The steady-state subset (orange) lies on the gain-corrected line, not the 1:1. Right: yaw rate against GPS course rate — the fitted gain sits below the ideal 1:1 throughout, in both directions of turn, which is the signature of a scale error rather than an offset or a clipped channel.*

**This invalidated a result I had already published.** A prior analysis attributed the same artifact to a 17 ms sensor timing lag. The two mechanisms are *degenerate* — both inject false slip proportional to yaw rate, and at K ≈ 0.44 a 14% gain error is equivalent to ~15 ms of lag. Re-running the estimator shows the gain correction alone removes the artifact (slope 0.62 → 0.036 °/g, R² 0.65 → 0.02), and applying **both** over-corrects to −0.43 °/g with R² back at 0.65. They cannot coexist. The scale error is the better-evidenced mechanism, and the timing hypothesis is retained in limitations as the alternative the data cannot fully exclude.

Also established here, from kinematics rather than channel names: `InlineAcc` is the true lateral channel (corr +0.974 vs `V·r`), `LateralAcc` is longitudinal and sign-inverted. Timing skew is **per-file** — 100 ms on this log — and must be re-derived for every file, never carried across.

`validate_channels.m`

---

## 2. Tyre model from FSAE TTC data

Magic Formula fitted per load bin to Calspan TIRF cornering data for the car's Hoosier, then reduced to load-dependent power laws:

- **Load sensitivity exponent k = 0.070** (μ falls as `Fz^-0.070`)
- **Camber: cornering stiffness −4.6%/deg, peak force only −0.8%/deg** — camber costs initial response far more than ultimate grip
- Cross-validated against Calspan's own tabulated cornering stiffness **within 6%** at every load

Fitted coefficients are **not published** (see Data, below).

---

## 3. Recovering the belt-to-track scaling from vehicle data

Rear axle lateral force predicted from the flat-trac tyre model through a full vehicle model — load transfer decomposed elastic/geometric/unsprung, speed-dependent aero download, per-axle camber — and compared against the same quantity measured on track.

**Recovered scaling 0.585 (telemetry aero) / 0.587 (handover aero)** against the documented FSAE belt-to-track factor of **0.667** — a 12% gap on the friction level.

The single number is not the evidence. **The ratio is flat across the whole slip-angle range** (standard deviation 0.031 on a mean of 0.598, over α_r = 1.3–6.7°). That can only happen if the curve *shape* is right, which requires the gyro correction, the channel identity, the sideslip estimate, the load-transfer split and the tyre fit to all be right simultaneously. A sloped ratio would indicate something upstream is broken. The remaining 12% gap on the friction *level* is stated, not explained away.

Including aerodynamic download matters: omitting it returns 0.627 and appears to agree with the documented factor to 6%. That agreement is an artifact of evaluating the tyre at too low a normal load. Both consistent aero source pairs are computed and neither is mixed with the other.

![Belt-to-track validation](figures/belt_to_track.png)

*Left: the unscaled flat-trac model over-predicts measured rear axle force by roughly 40%, and a single constant scaling brings it onto the measurement across the whole range. Right: the per-bin ratio. Its **flatness** is the result — a friction-level error scales every point equally, whereas an error in the tyre fit, the load transfer split or the sideslip estimate would tilt or curve this line. The 2.19° point sits 6% high; that bin has the fewest samples and the widest speed spread.*

**Roll stiffness is derived, not assumed.** From springs, motion ratios and track widths: front (no ARB) 534.8 and rear 309.1 N·m/deg sum to **843.9** against the team's stated 844; adding the 56 N·m/deg rear anti-roll bar residual gives **899.9** against the stated 899. Reproducing both from first principles is what makes the 0.594 elastic split citable.

`belt_to_track_validation.m`

---

## 4. Vehicle balance without a steering-angle channel

TBRe25 logs no steering angle — confirmed absent from the AiM channel list and from both CAN buses — so front axle slip angle, and the classical understeer gradient, are not directly measurable.

**Method:** in steady state both axle forces follow from geometry alone. Invert the validated tyre model to find the slip angle each axle needs at its own per-tyre loads, including lateral transfer, aero, camber and front toe-out.

The two halves are **not equally trustworthy**, and the README says so:
- **Rear: validated.** Model α_r matches telemetry-derived α_r to ~0.2° over 0.5–1.25 g.
- **Front: prediction.** No measurable α_f exists on this car.

![Sub-limit balance](figures/balance_sublimit.png)

*Sub-limit balance across toe convention, aero source and front weight fraction. Every as-run configuration sits below zero; at zero toe (blue) the same model returns mild understeer above 1.2 g.*

**Reading the curvature above 1.3 g.** The handover-aero curve (dashed) turns sharply upward and crosses into understeer, while the telemetry-aero curve keeps falling. This is the front axle approaching saturation under the lower-downforce assumption: as a tyre nears its peak, the slip angle required for each additional newton climbs steeply, so α_f rises away from α_r. The telemetry pair puts more download on the front, delaying that turn-up past the plotted range. Both curves are describing the same thing the utilisation table shows from the other direction — the two axles run out of grip at nearly the same lateral acceleration, so which one turns up first is decided by parameters that are not resolved. The steepening itself is real behaviour, not a numerical artifact, and it is the reason a balance verdict taken at 1.4 g would be far less robust than one taken at 1.0 g.

### Result: the limiting axle is indeterminate, and that is the finding

Expressed as **axle utilisation** (required force ÷ axle capacity), across the full uncertainty space:

| Front fraction | Aero source | ay_max | Front util | Rear util | Diff | First |
|---|---|---|---|---|---|---|
| 0.43 | handover | 1.438 | 98.1% | 100.0% | 1.9% | REAR |
| 0.43 | telemetry | 1.443 | 96.2% | 100.0% | 3.8% | REAR |
| **0.48** | handover | 1.448 | 100.0% | 99.6% | **0.4%** | FRONT |
| **0.48** | telemetry | 1.459 | 98.7% | 100.0% | **1.3%** | REAR |
| 0.53 | handover | 1.434 | 100.0% | 97.4% | 2.6% | FRONT |
| 0.53 | telemetry | 1.460 | 100.0% | 98.8% | 1.2% | FRONT |

**The axles sit within 0.4–3.8% of simultaneous saturation.** At the nominal parameter set the gap is 0.4% — a few newtons on two thousand. The limiting axle flips with both the front weight fraction and the choice of aero source, because the sign of a quantity smaller than its own uncertainty is not a result. The car is close to axle-balanced, and no axle verdict is supportable on this parameter set.

Sub-limit balance is weakly **oversteer-biased** (α_f − α_r = −0.20 to −0.26° at 1.0 g), and the **1° front toe-out is the mechanism, not a contributor**: at zero toe the same model returns mild understeer above 1.2 g. The balance is a setup choice, not a property of the car's mass and geometry.

![Predicted vs measured limit](figures/limit_vs_measured.png)

*Predicted limit against measured sustained maxima. Points above the model lines are binding; points below are not evidence of agreement.*

**Limit comparison is asymmetric evidence.** A measured sustained maximum is a *lower bound* on capability: bands where measurement falls below the model prove nothing, because the limit may not have been demanded. Only bands where measurement **exceeds** the model are binding. Three bands exceed the model, and the car demonstrably reached **1.589 g where the model predicts 1.464** — the model under-predicts by at least 9%. Matching the measured peak needs scaling ~0.62–0.64 against the 0.585 fitted at mid-range slip. Candidate cause: combined-slip contamination in endurance mid-range data suppressing measured force. **Open, unproven.**

`balance_analysis.m`

---

## 5. Corrected yaw acceleration envelope

An earlier version of this plot labelled envelope edges "front tyre limit" and "rear tyre limit" and concluded neutral-to-understeering behaviour. Both claims were wrong and have been removed.

- The edges were the outer hull of a scatter with **no mechanism behind them**. A Milliken Moment Method diagram separates axles because steer and body slip are swept independently; measured data gives only what the driver did.
- Yaw acceleration ≈ 0 at peak lateral g means **steady state, not neutrality**. Understeering, neutral and oversteering cars all settle to constant yaw rate mid-corner.

What the plot does support: the envelope **narrows 26%** between the 0.25–0.75 g and 1.25–1.75 g bands, consistent with tyres approaching lateral saturation leaving less capacity for yaw moment. Left–right asymmetry is **28% at Silverstone against 1.2% at Hockenheim** — same car, same sensors, which demonstrates the asymmetry is circuit layout rather than a car property.

![Yaw acceleration envelope](figures/yaw_envelope.png)

*Corrected envelope. The right-hand axis converts yaw acceleration to yaw moment through `N = I_z·ṙ`, peaking near ±880 N·m. The boundary is flat-topped through the mid-range and tapers at both ends — the earlier four-point convex hull imposed straight edges the data does not have, which is part of why the axle labels looked plausible.*

The unmeasured IMU mounting offset was tested rather than assumed away: sweeping it ±0.4 m longitudinally changes peak lateral acceleration by under 1.2% and envelope width by under 0.7%, so the result does not depend on it.

The envelope is left **open at the ends**: percentile bounds are only defined where samples exist, and closing it would mean drawing a boundary through a region the car never visited.

`plot_yaw_vs_lateral.m`

---

## Stated uncertainties

Nothing here is quoted without a provenance. The parameters that actually move conclusions:

| Parameter | Value used | Status |
|---|---|---|
| Front weight fraction | 0.48 | **Quoted verbally**; four-way document conflict (47.8/48/51/54). Flips the limiting axle, and moves the recovered scaling 0.60–0.65. **Corner weights would close it.** |
| Aero | 3.07 / 49% F and 3.58 / 53.7% F | Two consistent pairs, both computed, **never mixed**. Also flips the limiting axle. |
| Front toe | 1° toe-out | Total vs per-side unstated. Affects balance magnitude; drops out of the limit entirely. |
| Camber | −1.0° F / −0.645° R | Three-way conflict on record. |
| Roll centres | 48.6 / 53.3 mm | Two other pairs exist in team documents (23/34, 17/30). |
| Yaw inertia | ~100 kg·m² | Component summation from the team mass register, ±20%. Suppressed on the steady-state set. |
| IMU position | unmeasured | **Bounded and shown not to matter.** Sweeping the longitudinal offset ±0.4 m moves peak lateral acceleration only 1.988→2.012 g and envelope width 3.153→3.175 g. Measured lateral acceleration contains `x_s·ṙ`, which is correlated with the vertical axis and would distort the boundary rather than add noise — so this had to be checked, not assumed small. |

---

## Repository

```
validate_channels.m          sensor identity, scale, timing  — RUN FIRST
belt_to_track_validation.m   tyre model vs measured rear axle force
balance_analysis.m           balance and limiting axle by inversion
plot_yaw_vs_lateral.m        corrected yaw acceleration envelope
figures/
  channel_validation.png     gyro scale error, two independent references
  belt_to_track.png          flat-trac prediction vs measurement, ratio flat
  balance_sublimit.png       alpha_f - alpha_r across the uncertainty space
  limit_vs_measured.png      predicted limit vs measured sustained maxima
  yaw_envelope.png           corrected envelope, 26% taper
```

Run order: `validate_channels` → `belt_to_track_validation` → `balance_analysis`. The last two require `hoosier_16x75_10_R20.m` (not redistributed — see below). MATLAB R2016b+; no toolbox dependencies.

---

## Data

Tyre data from the **Formula SAE Tire Test Consortium**, tested at the **Calspan Tire Research Facility**. TTC data is licensed and **not redistributed here** — results appear in de-identified form only, per the participant agreement. The fitted tyre model file containing the coefficients is withheld for the same reason. TTC members can obtain the source data through their institution's membership.

Vehicle telemetry and the team parameter set are Team Bath Racing Electric property and are not included.

## Licence

MIT — covers repository code only, not the TTC data or team telemetry.
