\defgroup roofit_likelihood_math Mathematics of the RooFit likelihood
\ingroup Roofitmain
\brief Exact definition of the negative log-likelihood created by RooAbsPdf::createNLL(), including all constant terms.

Constant terms of a likelihood do not affect fit results, uncertainties, or
likelihood *ratios*, but they do shift the absolute value that
`RooAbsReal::getVal()` returns for the object created by
RooAbsPdf::createNLL(). Some of these constants depend only on the data and are
therefore not part of a serialized model (an HS3 file, for example), which makes
comparing absolute negative log-likelihood (NLL) values between %RooFit and
other tools error-prone.

This page writes the terms out explicitly and states which constants are
included and which are not.

### Notation

| Symbol | Meaning |
|---|---|
| \f$ \theta \f$ | the model parameters |
| \f$ x_i, w_i \f$ | observable values and weight of dataset entry \f$ i \f$ (\f$ w_i = 1 \f$ for unweighted data) |
| \f$ N = \sum_i w_i \f$ | sum of weights, `RooAbsData::sumEntries()` |
| \f$ n_j \f$ | observed content of bin \f$ j \f$ of a binned dataset |
| \f$ p(x\vert\theta) \f$ | the pdf, normalized over the fit observables in the fit range |
| \f$ \nu(\theta) \f$ | expected number of events, `RooAbsPdf::expectedEvents()` |
| \f$ \mu_j(\theta) \f$ | expected number of events in bin \f$ j \f$ |

The result of RooAbsPdf::createNLL() is the sum of a term describing the data
and a term describing the parameter constraints:

\f[
   \mathrm{NLL}(\theta) = \mathrm{NLL}_\mathrm{data}(\theta) + \mathrm{NLL}_\mathrm{constr}(\theta) \; .
\f]

### The data term

#### Per-event likelihood (default)

\f[
   \mathrm{NLL}_\mathrm{data}(\theta) = -\sum_i w_i \log p(x_i\vert\theta)
\f]

The pdf is always normalized over the observables of the fit, restricted to the
fit range if `Range()` was used. A RooDataHist is treated as a set of weighted
entries located at the bin centres, and \f$ p \f$ is the probability
*density* at the bin centre: the bin volume is **not** multiplied in. Passing
`IntegrateBins()` replaces the density at the bin centre by the density averaged
over the bin (the pdf integrated over the bin, divided by the bin volume), so it
does not introduce a bin volume factor either.

#### Extended term

If the extended term is active, the following is added:

\f[
   +\; \nu(\theta) - N \log \nu(\theta)
\f]

The extended term defaults to being active for every pdf that provides an
expected number of events (`RooAbsPdf::extendMode()` returning `CanBeExtended`
or `MustBeExtended`); it can be forced on or off with `Extended()`.

\warning The Poisson normalization \f$ \log N! \f$ of the extended term is
**not** included. %RooFit's extended likelihood therefore differs from the full
Poisson log-probability of \f$ N \f$ by the data-only constant
\f$ \log \Gamma(N+1) \f$.

#### Weighted data and `SumW2Error()`

The likelihood returned by `createNLL()` is never a weight-squared likelihood:
the formulas above always use \f$ w_i \f$, not \f$ w_i^2 \f$. The
`SumW2Error(true)` option of `RooAbsPdf::fitTo()` does not change the minimized
function either. After the fit, it temporarily switches the same NLL object into
weight-squared mode (`RooAbsArg::applyWeightSquared()`, replacing \f$ w_i \f$
by \f$ w_i^2 \f$ in the data term and rescaling the extended term by
\f$ \sum_i w_i^2 / \sum_i w_i \f$), runs a second HESSE, and uses the two
covariance matrices to build the corrected one. Only the uncertainties are
affected; the NLL value that an external tool would compare against is the
plain, weight-linear one.

#### Binned likelihood

If the pdf carries the `BinnedLikelihood` attribute, %RooFit interprets its
values directly as per-bin yields and replaces both terms above by a per-bin
Poisson term:

\f[
   \mathrm{NLL}_\mathrm{data}(\theta) = \sum_j \left[ \mu_j(\theta) - n_j \log \mu_j(\theta) + \log \Gamma(n_j+1) \right]
\f]

Bins with \f$ \mu_j = n_j = 0 \f$ contribute nothing. This is the code path used
for the channels of a RooStats::HistFactory model, which are RooRealSumPdf
objects flagged with that attribute. It is also used for HS3 models, because the
attribute survives a JSON round trip through the ROOT-specific
`misc.ROOT_internal.attributes` section of the file.

\warning Unlike the extended term, this expression **does** include the
data-only constants \f$ \log n_j! = \log \Gamma(n_j+1) \f$. Reading the same
model without the `BinnedLikelihood` attribute evaluates the per-event
likelihood with the extended term instead, which for a constant bin width
\f$ v \f$ gives a result that is lower by exactly
\f$ \sum_j \log \Gamma(n_j+1) - N \log v \f$.

### Constraint terms

Constraint pdfs \f$ C_c \f$ (collected from the model and from
`ExternalConstraints()`) enter as

\f[
   \mathrm{NLL}_\mathrm{constr}(\theta) = -\sum_c \log C_c(g_c\vert\theta) \; ,
\f]

where \f$ g_c \f$ are the global observables. Each constraint pdf is evaluated
**normalized** over the global observables (or, if no global observables were
declared, over the constrained parameters), so all of its normalization
constants are part of the NLL:

| Constraint | \f$ -\log C \f$ |
|---|---|
| RooGaussian \f$ (g \vert \theta, \sigma) \f$ | \f$ \frac{(g-\theta)^2}{2\sigma^2} + \log\left(\sqrt{2\pi}\,\sigma\right) \f$ |
| RooPoisson \f$ (g \vert \lambda(\theta)) \f$ | \f$ \lambda(\theta) - n \log \lambda(\theta) + \log \Gamma(n+1) \f$, with \f$ n = \lfloor g \rfloor \f$ unless RooPoisson::setNoRounding() was called |
| RooLognormal \f$ (g \vert m_0(\theta), \kappa) \f$, median parameterization | \f$ \frac{1}{2}\left(\frac{\log(g/m_0(\theta))}{\log\kappa}\right)^2 + \log\left(g\,\sqrt{2\pi}\,\log\kappa\right) \f$ |

The normalization integral runs over the *declared range* of the global
observable, not over \f$ (-\infty, \infty) \f$, so the closed forms above are
reached only in the limit of a range that is wide compared to the width of the
constraint. Since the integral depends on the constrained parameter as well,
this piece is not strictly constant, but it is usually negligible: for a Gaussian
constraint of width \f$ \sigma = 1 \f$ whose global observable is declared on
\f$ [-5, 5] \f$, the term \f$ \log\left(\sqrt{2\pi}\,\sigma\right) \f$ is
reduced by \f$ 5.7 \cdot 10^{-7} \f$ at \f$ \theta = 0 \f$ and by
\f$ 3.4 \cdot 10^{-6} \f$ at \f$ \theta = 0.5 \f$.

### Simultaneous fits

In a simultaneous fit, the index category counts as an additional observable and
the likelihood is built from the joint pdf
\f$ p(c, x\vert\theta) = p_c(x\vert\theta) / N_\mathrm{ch} \f$: the channel pdf
\f$ p_c \f$ is normalized over the observables of channel \f$ c \f$ alone, and
each of the \f$ N_\mathrm{ch} \f$ index states that have a pdf attached gets the
same weight \f$ 1/N_\mathrm{ch} \f$. On top of the sum of the per-channel terms,
the NLL therefore contains

\f[
   +\; N \log N_\mathrm{ch} \; ,
\f]

with \f$ N \f$ the sum of weights over all channels. Note that the
\f$ 1/N_\mathrm{ch} \f$ factor is applied by the likelihood only:
RooSimultaneous::getVal() returns the value of the pdf of the current index
state, without it. This term is dropped when bin offsetting is used.

### Offsetting

`Offset("initial")` subtracts the NLL value of the first evaluation to improve
the numerical precision of the minimization. It does **not** change the value
returned by `getVal()`: the shift is hidden from the user and only the minimizer
sees it (see RooAbsReal::setHideOffset()).

`Offset("bin")` does change the value. It subtracts the saturated model bin by
bin, which for a binned likelihood leaves

\f[
   \sum_j \left[ \mu_j(\theta) - n_j - n_j \log\frac{\mu_j(\theta)}{n_j} \right] \; ,
\f]

so that \f$ 2\,\mathrm{NLL} \f$ is approximately \f$ \chi^2 \f$-distributed. The
\f$ \log \Gamma(n_j+1) \f$ constants and the \f$ N \log N_\mathrm{ch} \f$ term
of a simultaneous fit cancel in this expression. For unbinned data, a histogram
template pdf built from the dataset with the current binning of the observables
takes the role of the saturated model.

### Summary

| Constant | Included in the NLL? |
|---|---|
| Normalization integral of the pdf | yes, the pdf is always normalized over the fit observables and the fit range |
| Bin volume for a RooDataHist fitted with a density | **no**, also not with `IntegrateBins()` |
| \f$ \log N! \f$ of the extended term | **no** |
| \f$ \log n_j! \f$ of a binned likelihood | **yes** (dropped by `Offset("bin")`) |
| Normalization of the constraint pdfs | yes, over the declared range of the global observables |
| \f$ N \log N_\mathrm{ch} \f$ of a simultaneous fit | yes (dropped by `Offset("bin")`) |
| Offset of `Offset("initial")` | not visible in `getVal()` |
