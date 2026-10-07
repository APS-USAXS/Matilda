# Gaussian-process desmearing of Bonse–Hart USAXS data

*Draft chapter for the instrument paper. Implementation theory, formulas, and
validation of the desmearing methods in Matilda. No general introduction;
begins at the forward model. Bibliographic details of the newest references
should be double-checked before submission (flagged where relevant).*

## 1. Slit-smearing forward model

The USAXS instrument uses Bonse–Hart crystal optics (Bonse & Hart, 1965), which
integrate scattering over a long slit in the direction perpendicular to the
scan. In the infinite-slit approximation the measured (smeared) intensity
$\tilde I(q)$ relates to the true, pinhole-equivalent intensity $I(q)$ by a
one-dimensional collimation integral (Lake, 1967; Ilavsky & Jemian, 2009):

$$
\tilde I(q) \;=\; \frac{1}{L}\int_{0}^{L} I\!\left(\sqrt{q^{2}+u^{2}}\right)\,\mathrm{d}u ,
\tag{1}
$$

where $L$ is the slit length in reciprocal-space units ($\text{Å}^{-1}$),
determined by the instrument geometry and stored per dataset (`slitLength`,
typically $L \approx 0.02$–$0.03\ \text{Å}^{-1}$, with $q_{\max}\approx 0.3\
\text{Å}^{-1}$). Equation (1) is a linear operator on $I$: discretising the
true intensity on a $q$-grid gives

$$
\tilde{\mathbf y} \;=\; \mathbf M\,\mathbf x, \qquad
\mathbf x \equiv I(q_j),\;\; \tilde{\mathbf y}\equiv \tilde I(q_i),
\tag{2}
$$

with $\mathbf M$ a fixed matrix obtained by trapezoidal quadrature of (1) and
linear interpolation in $\ln q$. The working grid is the measured grid extended
above $q_{\max}$ to $\sqrt{q_{\max}^2+L^2}$ so the integral is complete without
ad-hoc high-$q$ extrapolation.

## 2. The desmearing problem

Recovering $I$ from $\tilde I$ is a Fredholm problem of the first kind and is
ill-posed: direct inversion of $\mathbf M$ amplifies noise without bound. The
established solution used here for 25 years is the **Lake iterative method**
(Lake, 1967): starting from $\hat I^{(0)}=\tilde I$, iterate the multiplicative
correction $\hat I^{(k+1)} = \hat I^{(k)}\,\tilde I / \mathcal S[\hat I^{(k)}]$,
where $\mathcal S$ is the smearing operator (1). Lake is model-free and robust
on strongly scattering samples, but it has no regularisation: it reproduces
measurement noise and, where the smeared intensity approaches or crosses zero
(weak scattering, or over-subtraction of the instrumental background at mid-$q$),
the multiplicative update divides by near-zero values and produces non-physical
spikes and negative intensities.

An alternative widely used in small-angle scattering is to fit a regularised
model to the smeared data and evaluate it unsmeared — the indirect Fourier
transform (Glatter, 1977) and its Bayesian variants (Hansen, 2000) being the
classic examples for $P(r)$. Recently, Tung, Huang, Do and co-workers introduced
a Bayesian Gaussian-process (GP) formulation specifically for Bonse–Hart
slit geometry (Tung *et al.*, 2026a, 2026b), building on their model-free central
moment-expansion desmearing (Huang *et al.*, 2023, 2025). The implementation
below adapts that GP approach to the wide dynamic range of synchrotron USAXS.

## 3. Bayesian Gaussian-process desmearing

### 3.1 Multiplicative model

USAXS curves span many decades in both $q$ and $I$, with power-law regions and
occasional near-zero or negative excursions. To handle this we parametrise the
desmeared intensity **multiplicatively**,

$$
I(q) \;=\; \mu_0(q)\,\exp\!\big(s(q)\big),
\tag{3}
$$

where $\mu_0(q)>0$ is a smooth base trend (below) and $s(q)$ is a dimensionless
log-correction. The multiplicative form is essential: at low $q$ a steep power
law is strongly suppressed by the slit, so the true intensity can exceed the
smeared value by two to three orders of magnitude, yet $s=\ln(I/\mu_0)$ remains a
smooth, $O(1)$–$O(10)$ function of $\ln q$. The parametrisation also guarantees
$I(q)>0$, so the reconstruction is always physical.

### 3.2 Prior

We place a zero-mean GP prior on the log-correction as a function of
$x=\ln q$ (Rasmussen & Williams, 2006):

$$
s(x)\sim \mathcal{GP}\!\big(0,\,k(x,x')\big).
\tag{4}
$$

Working in $\ln q$ makes the prior scale-free across decades. Two covariance
kernels are provided, with correlation length $\ell = \ell_{\mathrm{dec}}\ln 10$
(so $\ell_{\mathrm{dec}}$ is expressed in decades of $q$) and amplitude
$\sigma_f$:

$$
k_{\text{RBF}}(x,x') = \sigma_f^2\exp\!\Big(-\tfrac{1}{2}\big(\tfrac{|x-x'|}{\ell}\big)^2\Big),
\qquad
k_{3/2}(x,x') = \sigma_f^2\Big(1+\tfrac{\sqrt3\,|x-x'|}{\ell}\Big)\exp\!\Big(-\tfrac{\sqrt3\,|x-x'|}{\ell}\Big).
\tag{5}
$$

The squared-exponential (RBF) kernel yields infinitely differentiable sample
paths — a very strong smoothness assumption that suppresses noise aggressively
but erases genuine sharp features (correlation peaks, form-factor oscillations).
The Matérn-3/2 kernel yields once-differentiable paths, admitting local structure
while still rejecting noise, and is the default. $\ell_{\mathrm{dec}}$ (default
$0.5$) is the principal user control.

### 3.3 Likelihood and error model

The measured smeared intensities are taken as Gaussian about the forward model,

$$
\tilde y_i \sim \mathcal N\!\big((\mathbf M \mathbf x)_i,\ \sigma_i^2\big),
\tag{6}
$$

with the linear operator of (2) evaluated in intensity (not log) space, so
near-zero and negative $\tilde y_i$ are handled naturally. The variances use the
measured uncertainties $\mathrm{Idev}_i$ with two floors,

$$
\sigma_i \;=\; \max\!\Big(\mathrm{Idev}_i,\ \ f_{\mathrm{abs}}\,\widetilde{\mathrm{Idev}},\ \ f_{\mathrm{rel}}\,|\tilde y_i|\Big),
\tag{7}
$$

where $\widetilde{\cdot}$ is the median, $f_{\mathrm{abs}}=0.01$ prevents
zero-weight points, and $f_{\mathrm{rel}}=0.01$ imposes a **1% relative floor**.
The 1% floor reflects the instrument's true relative-error floor and is
essential: it correctly down-weights near-zero points (where $f_{\mathrm{rel}}|\tilde y_i|\to0$
and the absolute floor dominates, so noise is treated as noise, not signal), and
it regularises datasets whose stored uncertainties are unphysically small.

### 3.4 MAP estimate and uncertainty

With prior (4)–(5), covariance matrix $\mathbf K_{jk}=k(x_j,x_k)$, and noise
$\boldsymbol\Sigma=\mathrm{diag}(\sigma_i^2)$, the maximum-a-posteriori
log-correction minimises

$$
\Phi(\mathbf s) \;=\; \big(\tilde{\mathbf y}-\mathbf M\mathbf x\big)^{\!\top}\boldsymbol\Sigma^{-1}\big(\tilde{\mathbf y}-\mathbf M\mathbf x\big) \;+\; \mathbf s^{\!\top}\mathbf K^{-1}\mathbf s,
\qquad \mathbf x = \boldsymbol\mu_0\odot e^{\mathbf s}.
\tag{8}
$$

Because $\mathbf x$ depends nonlinearly on $\mathbf s$ through the exponential,
(8) is minimised by damped Gauss–Newton iteration. With Jacobian
$\mathbf J = \mathbf M\,\mathrm{diag}(\mathbf x)$ and residual
$\mathbf r = \tilde{\mathbf y}-\mathbf M\mathbf x$, each step solves

$$
\big(\mathbf J^{\!\top}\boldsymbol\Sigma^{-1}\mathbf J + \mathbf K^{-1}\big)\,\boldsymbol\delta
\;=\; \mathbf J^{\!\top}\boldsymbol\Sigma^{-1}\mathbf r - \mathbf K^{-1}\mathbf s,
\tag{9}
$$

followed by a backtracking line search on $\Phi$ (guaranteeing descent) and a
hard bound $|s|\le\ln(10^8)$. Iteration starts at $\mathbf s=\mathbf 0$ (i.e.
$I=\mu_0$) and typically converges in fewer than ten steps.

At the optimum the Laplace (Gaussian) posterior for $\mathbf s$ has covariance

$$
\mathrm{Cov}(\mathbf s) \;=\; \big(\mathbf J^{\!\top}\boldsymbol\Sigma^{-1}\mathbf J + \mathbf K^{-1}\big)^{-1},
\tag{10}
$$

whose diagonal $\sigma_{s,j}^2$ gives a multiplicative credibility interval on
the desmeared intensity,

$$
I_j\,\exp(-z\,\sigma_{s,j}) \;\le\; I_j \;\le\; I_j\,\exp(+z\,\sigma_{s,j}),
\tag{11}
$$

reported at $z=2$ (≈95%). The interval widens automatically in weak or
under-constrained regions (the low-$q$ edge, the extrapolated high-$q$ tail),
providing an honest, data-driven statement of where the reconstruction is
determined by the data versus by the prior — information the Lake method cannot
provide. The posterior 1σ is stored as the desmeared uncertainty.

### 3.5 Base trend and method variants

The base trend $\mu_0$ is a strictly positive, smooth envelope of the data,
obtained by a geometric-mean (log-space) moving average of the smeared
intensities after clipping to the noise floor of (7); log-space smoothing
preserves multi-decade power laws that a linear average would distort. Two
variants are exposed: the default (`gp`) uses this data envelope, while
`gp_lake_mean` uses a smoothed Lake reconstruction as $\mu_0$, so the GP reverts
to Lake's robust behaviour at the under-constrained low-$q$ edge while still
smoothing elsewhere.

### 3.6 Slit-resolution-aware regularization

A distinct failure mode appears on **noisy data with a flat (Guinier) plateau
when the slit length is comparable to or larger than q_min**. There the low-q
desmeared intensity lies in the near-null space of the operator: since a
measured point at $q$ integrates $I$ over $[q,\sqrt{q^2+L^2}]$, on a flat
plateau all low-q measurements return essentially the same value, so many
different low-q shapes reproduce the data equally well (numerically, the
low-q block of $\mathbf M$ is ill-conditioned, and a low-q oscillation of unit
amplitude is damped to $\sim\!25\%$ in the smeared data). A stationary GP fills
this null space with a coherent oscillation at $\sim\ell$ — a Gibbs-type ringing
artifact. Lake avoids coherent oscillation because it is a pointwise iteration
with no global basis; its noise stays incoherent and averages to the plateau.

The cure is to forbid structure the slit cannot resolve. The slit sets a
q-dependent resolution: in $\ln q$ its half-width is

$$
W(q) \;=\; \tfrac12\ln\!\big(1+(L/q)^2\big),
\tag{12}
$$

which is negligible at high $q$ and grows as $q\to0$. We therefore use a
**non-stationary** GP whose local correlation length is
$\ell(q)=\max(\ell_0,\,\alpha\,W(q))$ ($\alpha\simeq1$), implemented by warping
the input coordinate $u(x)=\int^{x}\!\mathrm{d}x'/\ell(x')$ and applying a
unit-length stationary kernel in $u$. This lengthens the prior correlation
exactly where the slit smears hardest, so the low-q null space is smoothed to
the physical plateau (no ringing), while high-q resolution is untouched.
Crucially it suppresses only short-wavelength *oscillations*: large, smooth,
data-constrained corrections — e.g. a genuine low-q upturn or correlation peak —
are unaffected, because there the likelihood dominates the prior. (On a test
sample with a real peak at $q\approx10^{-3}$ the recovered peak height is
unchanged with the resolution-aware prior on or off.)

A complementary safeguard operates directly on the GP output. Where the total
column weight of $\mathbf M$ is small (the under-constrained low-$q$ block on a
flat plateau), the desmeared intensity is not allowed to fall below the local
*smeared* plateau level — the correct value for a flat region — which removes the
artefactual low-$q$ ramp. A genuine rising low-$q$ power law (resolution-limited,
strong signal) is untouched, because there the desmeared curve exceeds the
smeared data and the floor never binds. This distinguishes the two physically
different low-$q$ regimes automatically: signal-limited (flat, needs smoothing)
versus resolution-limited (rising, well-determined).

For the noisiest flat-plateau data a simpler alternative is offered:
**smoothed Lake** — the Lake reconstruction followed by a geometric-mean
smoothing. It inherits Lake's null-space-safe shape and removes the incoherent
point noise, at the cost of the credibility band. No single method is optimal
for all data: the GP wins on weak / near-zero / over-subtracted data (where Lake
spikes), and smoothed Lake on noisy flat plateaus (where a regularized inversion
would ring).

## 4. Software implementation

The method is implemented in `matilda/desmearing_methods.py` (NumPy only). A
single dispatcher, `desmear_dispatch`, shares the signature and 4-tuple return of
the existing `desmearData`, so it is a drop-in replacement at the reduction call
site; `method='lake'` reproduces the legacy Lake path bit-for-bit. The slit
operator $\mathbf M$ is assembled once per curve on the extended log-$q$ grid;
the dense linear systems (9)–(10) are $O(N^3)$ in the number of points, which is
negligible for rebinned USAXS ($N\!\approx\!500$). The GP methods are selectable
from the reduction GUI (method, kernel, and $\ell_{\mathrm{dec}}$); Lake remains
the default for scripted and automated reduction.

## 5. Validation

The implementation was validated on data from the APS USAXS instrument
(Ilavsky *et al.*, 2018).

**Consistency.** On strongly scattering samples the GP reconstruction agrees
with Lake to ≈1% (median absolute log-intensity difference) across the full
$q$-range, confirming the method is unbiased where the data are informative.

**Robustness (batch).** A random sample of 80 datasets (of ~3100) spanning many
users and sample types was processed with Lake, GP-Matérn, and GP-RBF. Lake
produced spikes or negative intensities on **31/80 (39%)** of datasets —
verified to be present in the stored production data, not an artefact of
re-processing. After correcting two numerical issues (log-space base-trend
construction; the 1% relative-error floor, §3.3), the GP produced smooth,
strictly positive reconstructions on **all 80** datasets (maximum point-to-point
log-intensity jump 1.3, versus Lake's ≈690 where it reaches zero), with a median
credibility half-width of ±12%.

**Weak / over-subtracted data.** The characteristic mid-$q$ failure mode of the
instrument — where sample and background scattering are comparable and
subtraction drives the excess intensity to zero or slightly negative — was
reproduced by injecting an over-subtraction into a test dataset (38% of points
forced $\le 0$). Lake returned 192 negative points and spikes spanning >30
decades; the GP returned a smooth, positive curve tracking the true low-$q$
structure, with the credibility band widening to ≈±350% through the weak region
— the physically expected "small, smooth, magnitude-uncertain" result.

**Feature preservation.** On a weakly scattering sample with a real correlation
peak at $q\approx 10^{-3}\ \text{Å}^{-1}$, GP-Matérn reconstructed the peak,
whereas GP-RBF smoothed it away — consistent with the kernels' differentiability
(§3.2) and the basis for making Matérn-3/2 the default.

**Flat-plateau null space.** On noisy samples with a flat Guinier plateau and
$L\gtrsim q_{\min}$ (here $L/q_{\min}\approx 6$–$7$, relative noise 10–26%), the
stationary GP placed coherent low-$q$ oscillations on the plateau. This is a
genuine null-space artefact: the low-$q$ block of $\mathbf M$ is ill-conditioned
(condition number $\sim\!6\times10^4$) and a unit-amplitude low-$q$ oscillation
is damped to $\sim\!26\%$ in the smeared data, so the data cannot constrain the
low-$q$ shape and a regularized inversion rings there. The resolution-aware prior
(§3.6) removes the oscillation — the low-$q$ block behaves as a smooth plateau —
and, independently, smoothed Lake reproduces the plateau directly. Lake itself
does not oscillate (its noise is incoherent) but is visibly noisy; smoothing it
recovers the physical plateau.

## 6. Method selection in practice

The failure analysis yields a clear, physically-grounded guide — no single method
is optimal across all data, but the operator's conditioning predicts which to use:

| Data character | Recommended | Rationale |
|---|---|---|
| Strong scattering, well-defined features | GP (Matérn-3/2) | noise suppression + feature preservation + uncertainties |
| Weak / near-zero / over-subtracted intensity | GP (Matérn-3/2) | smooth and strictly positive where Lake spikes; honest band |
| Noisy flat Guinier plateau, $L\gtrsim q_{\min}$ | Smoothed Lake, or GP with the resolution-aware prior | avoids null-space ringing; recovers the physical plateau |
| Featureless and intrinsically smooth | GP (RBF) | maximal noise suppression (safe only without sharp features) |
| Automated / unattended reduction | Lake (current default) | conservative; no per-sample tuning |

The unifying principle is the slit operator's conditioning: where it is
well-posed the GP is unbiased and adds value (noise suppression, uncertainty);
where a flat plateau drives it into a near-null space, the honest choice is to
smooth rather than to invert.

## 7. Status and limitations

The GP methods are offered as GUI-selectable alternatives; Lake remains the
production default pending broader field testing. Current limitations: the
desmeared entry does not yet record the method used (provenance); the dense
solve is impractical for step scans retaining thousands of points; the discrete
slit operator differs from the legacy smearing routine by ≈20% at the high-$q$
tail and should be reconciled before cross-method $\chi^2$ is used
quantitatively; and automatic selection of kernel and correlation length by GP
marginal likelihood (evidence) is a natural but unimplemented extension.

## References

*(Verify exact pages/authors for the 2024–2026 references before submission.)*

- Bonse, U. & Hart, M. (1965). *Appl. Phys. Lett.* **7**, 238–240.
- Glatter, O. (1977). *J. Appl. Cryst.* **10**, 415–421.
- Hansen, S. (2000). *J. Appl. Cryst.* **33**, 1415–1421.
- Huang, G.-R., Tung, C.-H., Chen, M.-Z., Porcar, L., Shinohara, Y., Wildgruber,
  C. U., Do, C. & Chen, W.-R. (2023). *J. Appl. Cryst.* **56**, 1537–1543.
  (Model-free desmearing by central-moment expansion, 1D.)
- Huang, G.-R., Tung, C.-H., *et al.* (2025). *Desmearing two-dimensional SANS
  data by central moment expansions.* arXiv:2502.13488 (J. Appl. Cryst.,
  in press).
- Ilavsky, J. & Jemian, P. R. (2009). *J. Appl. Cryst.* **42**, 347–353. (Irena.)
- Ilavsky, J., Zhang, F., Andrews, R. N., Kuzmenko, I., Jemian, P. R., Levine,
  L. E. & Allen, A. J. (2018). *J. Appl. Cryst.* **51**, 867–882. (APS USAXS
  facility.)
- Lake, J. A. (1967). *Acta Cryst.* **23**, 191–194. (Iterative slit correction.)
- Rasmussen, C. E. & Williams, C. K. I. (2006). *Gaussian Processes for Machine
  Learning.* MIT Press.
- Tung, C.-H., Huang, G.-R., Do, C., *et al.* (2026a). *Desmearing Bonse–Hart
  USANS data using Bayesian Gaussian process regression.* *J. Chem. Phys.*
  **164**, 204102.
- Tung, C.-H., Huang, G.-R. & Do, C. (2026b). *A Bayesian desmearing algorithm
  for Bonse–Hart USANS with anisotropic scattering.* *J. Chem. Phys.* **164**,
  164108.
- Ryukhtin, V., Len, A., Almásy, L., Juszyńska-Gałązka, E., Zając, W. & Tomchuk,
  O. (2024). *J. Appl. Cryst.* **57**, 1551–1556. (Pinhole-SANS-based desmearing
  of slit USANS.)
