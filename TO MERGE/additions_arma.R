### Forecasting AR(p) Models

One of the most interesting things in time series analysis is to predict the future unobserved values based on the values that have been observed up to now. However, this is not possible if the underlying (parametric) model is unknown, thus in this section we assume the time series $X_t$ is stationary and its model is known. In particular, we denote forecasts by $X^{T}_{T+m}$, where $n$ represents the data points collected (e.g. $\mathbf{X} = (X_{1}, X_{2}, \cdots , X_{T-1}, X_T)$) and $m$ represents the amount of points in the future we wish to predict. So, $X^{T}_{T+1}$ represents a one-step-ahead prediction $X_{T+1}$ given data $(X_{1}, X_{2}, \cdots, X_{T-1}, X_{T})$.

To obtain forecasts, we rely upon the best linear prediction (BLP). Recall that the BLP Definition \@ref(def:blp) is mentioned when we calculate the PACF of AR model. In general, the best approach to finding the BLP is to use the Theorem \@ref(thm:projtheo).

```{theorem, label="projtheo", name="Projection Theorem"}
Let $\mathcal{M} \subset \mathcal{L}_2$ be a closed linear subspace of Hibert space. 
For every $y \in \mathcal{L}_2$, there exists a unique element $\hat{y} \in \mathcal{M}$ that
minimizes $||y - l||^2$ over $l \in \mathcal{M}$. This element is uniquely 
determined by the requirements

1. $\hat{y} \in \mathcal{M}$ and 
2. $(y - \hat{y}) \perp \mathcal{M}$.
```

<!-- \mbox{(trivial for time series with zero mean)} -->
  
  The projection theorem naturally leads to an equivalent way to find the best linear 
predictor by solving the prediction equations,
\[ \mathbb{E}(X_{t+h} - \hat{X}_{t+h}) = 0,\]
and
\[ \mathbb{E} [(X_{t+h} - \hat{X}_{t+h})X_j ] = 0, \mbox{ for } i = 1, \dots, t.\]

If we denote $\mathbb{E}(X_{i}, X_{j})$ as $\gamma(|i - j|)$, the predition equations can
be represented in the following form

\begin{equation}
\begin{aligned}
\begin{pmatrix}
\gamma(0) & \gamma(1) & \cdots & \gamma(T-1) \\
\gamma(1) & \gamma(0) & \cdots & \gamma(T-2) \\
\vdots & \vdots & \ddots & \vdots \\
\gamma(T-1) & \gamma(T-2) & \cdots &\gamma(0)
\end{pmatrix}_{T \times T}
\begin{pmatrix}
\alpha_1 \\
\vdots \\
\alpha_T
\end{pmatrix}_{T \times 1}
&=
  \begin{pmatrix}
\gamma(T+h-1)  \\
\vdots \\
\gamma(h)
\end{pmatrix}_{T \times 1} \\
\Gamma_T \mathbf{\alpha}_T  &= \mathbf{\gamma}_T
\end{aligned}
\end{equation}.

Assuming that $\Gamma_T$ is non-singular, then the values of $\mathbf{\alpha}_T$ are
given as:
  
  $$\mathbf{\alpha}_T  = \Gamma^{-1}_T\mathbf{\gamma}_T$$
  
  
  ### Inference for AR(p) Models
  
  For all the above methods, it would be necessary to understand how "precise" 
their estimates are. To do so we would need to obtain confidence intervals for
these estimates and this can be done mainly in two manners:
  
  - using the asymptotic distribution of the parameter estimates;
- using parametric bootstrap.

The first approach consists in using the asymptotic distribution of the 
estimators presented earlier to deliver approximations of the confidence
intervals which get better as the length of the observed time series increases. 
Hence, if for example we wanted to find a 95% confidence interval for the 
parameter $\phi$, we would use the quantiles of the normal distribution 
(given that all methods presented earlier present this asymptotic distribution).
However, this approach can present some drawbacks, one of which its
behaviour when the parameters are close to the boundaries of the parameter space. 
Suppose we consider a realization of length $T = 100$ of the following AR(1) model:
  
  \[X_t = 0.96 X_{t-1} + W_t, \;\;\;\; W_t \sim \mathcal{N}(0,1),\]

which is represented in Figure \@ref(fig:simAR1ci)

``` {r simAR1ci, cache = TRUE, fig.height= 4.5, fig.width = 9, fig.cap = "AR(1) with $\\phi$ close to parameter bound"}
set.seed(7)
x = gen_gts(n = 100, AR1(0.96, 1))
plot(x)
```

It can be seen that the parameter $\phi = 0.96$ respects the condition for 
stationarity (i.e.$\left| \phi \right| < 1$) but is very close to its boundary.
Using the MLE and the GMWM estimators, we first estimate the parameters and
then compute confidence intervals for $\phi$ using the asymptotic 
normal distribution.

```{r asymAR1estci, cache = TRUE}
# Compute the parameter estimates using MLE
fit.ML = arima(x, order = c(1,0,0), include.mean = FALSE)
c("phi" = fit.ML$coef, "se" = sqrt(fit.ML$var.coef))

# Construct asymptotic confidence interval for phi
fit.ML$coef + c(-1,1)*1.96*as.numeric(sqrt(fit.ML$var.coef))

# Compute the parameter estimates with inference using GMWM
fit.gmwm = gmwm2::gmwm(AR1(), x)
summary(fit.gmwm, inference = TRUE)$estimate
```


From the estimation summary, we can notice that both confidence intervals contain 
values that would make the AR(1) non-stationary (i.e. values of $\phi$ larger
                                                 than 1). However, these confidence intervals are based on the (asymptotic) distributions of $\hat{\phi}$ which are shown in Figure
\@ref(fig:asymIC).

```{r asymIC, cache = TRUE, echo = FALSE, fig.height = 4.5, fig.width = 9, fig.cap = "Estimated asymptotic distribution of $\\hat{\\phi}$ for MLE and GMWM parameter estimates. The dashed vertical line represents the true value of $\\phi$, the solid line denotes the upper bound of the parameter space for $\\phi$ and the vertical ticks represent the limits of the 95% confidence intervals for both methods."}
# Define colors
gg_color_hue = function(n, alpha = 1) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
}

colors = gg_color_hue(6, alpha = 0.2)
colors2 = gg_color_hue(6, alpha = 1)
phi.sd.gmwm = summary(fit.gmwm)$estimate[,"SE"][1]

par(mfrow = c(1,2))
xx = seq(from = 0, to = 1.2, length.out = 10^3)
yy.ML = dnorm(xx, fit.ML$coef, sqrt(fit.ML$var.coef))
plot(NA, xlim = c(0.8,1.04), ylim = range(yy.ML), 
     xlab = expression(phi), ylab = "Density", main = "MLE")
grid()
abline(v = 1, lwd = 2, col = "grey60")
abline(v = 0.96, lwd = 2, lty = 2, col = "grey60")
polygon(c(xx,rev(xx)), c(rep(0,length(yy.ML)), rev(yy.ML)), border = NA, col = colors[1]) 
points(qnorm(0.975,fit.ML$coef, sqrt(fit.ML$var.coef)), 0, cex = 3, 
       col = colors2[1], pch = "|")
points(qnorm(0.025,fit.ML$coef, sqrt(fit.ML$var.coef)), 0, cex = 3, 
       col = colors2[1], pch = "|")

yy.gmwm = dnorm(xx, fit.gmwm$estimate[1], phi.sd.gmwm)
plot(NA, xlim = c(0.7,1.2), ylim = range(yy.gmwm), 
     xlab = expression(phi), ylab = "Density", main = "GMWM")
grid()
abline(v = 1, lwd = 2, col = "grey60")
abline(v = 0.96, lwd = 2, lty = 2, col = "grey60")
polygon(c(xx,rev(xx)), c(rep(0,length(yy.gmwm)), rev(yy.gmwm)), border = NA, col = colors[2])
points(qnorm(0.975, fit.gmwm$estimate[1], phi.sd.gmwm), 0, cex = 3, 
       col = colors2[2], pch = "|")
points(qnorm(0.025, fit.gmwm$estimate[1], phi.sd.gmwm), 0, cex = 3, 
       col = colors2[2], pch = "|")
```

Therefore, if we estimate a stationary AR(1) model, it would be convenient to have more "realistic" confidence intervals that give limits for a stationary AR(1) model. A viable solution for this purpose is to use parametric bootstrap, detailed in Theorem \@ref(thm:parabootstrap). Indeed, parametric bootstrap takes the estimated parameter values and uses them in order to simulate from an AR(1) process. For each simulation, the parameters are estimated again and stored.  Finally, the empirical quantiles at $\alpha/2$ and $1-\alpha/2$ are taken of the saved estimated parameter values to obtrain a confidence interval which does not suffer from boundary problems. The code below gives an example of how this confidence interval is built based on the same estimation procedure but using parametric bootstrap (using $B = 10000$ bootstrap replicates).

``` {r paraAR1ciest, cache = TRUE, warning = FALSE, fig.height= 6, fig.width= 9}
# Number of Iterations
B = 10000

# Set up storage for results
est.phi.gmwm = rep(NA,B)
est.phi.ML = rep(NA,B)

# Model generation statements
model.gmwm = AR1(fit.gmwm$estimate[1], fit.gmwm$estimate[2])
model.mle = AR1(fit.ML$coef, fit.ML$sigma2)

# Begin bootstrap
for(i in seq_len(B)){
  
  # Set seed for reproducibility
  set.seed(B + i)
  
  # Generate process under MLE parameter estimate
  x.star = gen_gts(100, model.mle)
  
  # Attempt to estimate phi by employing a try
  est.phi.ML[i] = tryCatch(arima(x.star, order = c(1,0,0), include.mean = FALSE)$coef, 
                           error = function(e) NA)
  
  # Generate process under GMWM parameter estimate
  x.star = gen_gts(100, model.gmwm)
  est.phi.gmwm[i] = gmwm2::gmwm(model.gmwm, x.star)$estimate[1]
}
```

``` {r paraAR1cigraphs, cache = TRUE, warning = FALSE, echo = FALSE, fig.height = 6, fig.width = 9, fig.cap = "Estimated parametric distributions of $\\hat{\\phi}$ for MLE and GMWM parameter estimates. The bar plots represent the results from the parametric bootstrap, the densities represent the estimate asymptotic distribution, the vertical line represents the true value of $\\phi$, the solid lines denotes the upper bound of the parameter space for $\\phi$ and the vertical ticks represent the limits of the 95% confidence intervals under both approaches for the estimation methods."}
par(mfrow = c(1,2))

hist.phi.ML = hist(est.phi.ML, plot = FALSE, breaks = 20)
hist.phi.ML$counts = hist.phi.ML$counts/sum(hist.phi.ML$counts)/diff(hist.phi.ML$mids)[1]
plot(NA, xlim = c(0.65,1.04), ylim = range(yy.ML), 
     xlab = expression(phi), ylab = "Density", main = "MLE")
grid()
abline(v = 1, lwd = 2, col = "grey60")
abline(v = 0.96, lwd = 2, lty = 2, col = "grey60")
plot(hist.phi.ML, xlim = c(0.6,1.1), ylim = range(yy.ML), 
     add = T, col = colors[3], border = "grey60")
polygon(c(xx,rev(xx)), c(rep(0,length(yy.ML)), rev(yy.ML)), border = NA, col = colors[1]) 

points(qnorm(0.975,fit.ML$coef, sqrt(fit.ML$var.coef)), 0, cex = 3, col = colors2[1], pch = "|")
points(qnorm(0.025,fit.ML$coef, sqrt(fit.ML$var.coef)), 0, cex = 3, col = colors2[1], pch = "|")

points(quantile(est.phi.ML, 0.975, na.rm = T), 0, cex = 3, col = colors2[3], pch = "|")
points(quantile(est.phi.ML, 0.025, na.rm = T), 0, cex = 3, col = colors2[3], pch = "|")


hist.phi.GMWM = hist(na.omit(est.phi.gmwm), plot = FALSE, breaks = 20)
hist.phi.GMWM$counts = hist.phi.GMWM$counts/sum(hist.phi.GMWM$counts)/diff(hist.phi.GMWM$mids)[1]
plot(NA, xlim = c(0.53,1.2), ylim = range(hist.phi.GMWM$counts), 
     xlab = expression(phi), ylab = "Density", main = "GMWM")
grid()
abline(v = 1, lwd = 2, col = "grey60")
abline(v = 0.96, lwd = 2, lty = 2, col = "grey60")
plot(hist.phi.ML, xlim = c(0.53,1.2), ylim = range(hist.phi.GMWM$counts), 
     add = T, col = colors[5], border = "grey60")
polygon(c(xx,rev(xx)), c(rep(0,length(yy.gmwm)), rev(yy.gmwm)), border = NA, col = colors[2]) 

points(qnorm(0.975,fit.gmwm$estimate[1], phi.sd.gmwm), 0, cex = 3, col = colors2[2], pch = "|")
points(qnorm(0.025,fit.gmwm$estimate[1], phi.sd.gmwm), 0, cex = 3, col = colors2[2], pch = "|")

points(quantile(est.phi.gmwm, 0.975, na.rm = T), 0, cex = 3, col = colors2[5], pch = "|")
points(quantile(est.phi.gmwm, 0.025, na.rm = T), 0, cex = 3, col = colors2[5], pch = "|")
```

In Figure \@ref(fig:paraAR1cigraphs), we compare the estimated densities
for $\hat{\phi}$ using asymptotic results and bootstrap techniques for the MLE 
and the GMWM estimator. It can be observed that the issue that arose previously 
with unrealistic confidence intervals have been resolved as the interval regions
lie entirely within the boundaries of the parameter space.

To emphasize that the effectiveness of parametric bootstrap, let us
consider one additonal example of an AR(2) process. For this example, discussion
will focus on observing the behavior exhibited solely by the MLE. Having said
this, the AR(2) process is defined to be:
  
  \begin{equation}
{X_t} = {1.98}{X_{t - 1}} - {0.99}{X_{t - 2}} + {W_t}
(\#eq:ciar2)
  \end{equation}
  
  where $W_t \sim WN(0, 1)$. The generated process is displayed in Figure \@ref(fig:CIAR2data). 
  
  ```{r CIAR2data, cache = TRUE, fig.cap="Generated AR(2) Process"}
  set.seed(432)
  Xt = gen_gts(500, AR(phi = c(1.98, -0.99), sigma2 = 1))
  plot(Xt)
  mod = arima(Xt, c(2,0,0), include.mean = FALSE)
  mod
  ```
  
  With the model's coefficients readily available, the parametric bootstrap
  is able to be constructed. As before, the bootstrap implementation mirrors prior
  versions with one cavaet regarding the storage of $\phi$ being two dimensional and,
  thus, requiring a `matrix` to store the result instead of a traditional atomic
  vector in _R_. 
  
  ```{r CIAR, cache = TRUE, warning = FALSE}
  B = 10000
  est.phi.ML = matrix(NA, B, 2)
  model = AR(phi = mod$coef, sigma2 = mod$sigma2)
  
  for(i in seq_len(B)) {
  set.seed(B + i)              # Set Seed for reproducibilty
  
  x.star = gen_gts(500, model) # Simulate the process 
  
  # Obtain parameter estimate with protection
  # that defaults to NA if unable to estimate
  est.phi.ML[i,] = tryCatch(arima(x.star, order = c(2,0,0),
  include.mean = FALSE)$coef,
  error = function(e) NA)
  }
  ```
  
  ```{r ciAR2d, cache = TRUE, echo = FALSE, fig.height= 4.5, fig.width= 9, warning = FALSE, fig.cap="Estimated distributions of $\\hat{\\phi_1}$ and $\\hat{\\phi_2}$ of the MLE using asymptotic and parametric bootstrap techniques. The colored contours represent the density of distribution and the dark grey lines represent the boundary constraints of $\\left|\\phi_2\\right|<1$ and $\\phi_2 = 1 - \\phi_1$."}
  # Graphing 
  
  # Requires use of MASS
  library("MASS")
  library("RColorBrewer")
  est.phi.ML = na.omit(est.phi.ML)
  z = kde2d(est.phi.ML[,1],est.phi.ML[,2], n=50)
  k = 11
  my.cols = rev(brewer.pal(k, "RdYlBu"))
  
  bivn = mvrnorm(100000, mu = mod$coef, Sigma = matrix(c(mod$var.coef), 2))
  bivn.kde = kde2d(bivn[,1], bivn[,2], n = 50)
  
  par(mfrow = c(1,2))
  plot(NA, xlim = c(1.96,2), ylim = c(-1.02,-0.97), xlab = expression(phi[1]),
  ylab = expression(phi[2]), cex.lab = 1.5, main = "Asymptotic")
  grid()
  
  
  # Adding boundary of constraint |phi_2| < 1
  abline(h = c(-1,1), lty = 2, col = "darkgrey")
  
  # Adding boundary of constraint phi_2 = 1 - phi_1 
  phi1 = seq(from = -2, to = 2, length.out = 10^3)
  phi2.c1 = 1 - phi1
  lines(phi1, phi2.c1, lty = 2, col = "darkgrey")
  
  # Adding boundary of constraint phi_2 = 1 + phi_1 
  phi1 = seq(from = -2, to = 2, length.out = 10^3)
  phi2.c2 = 1 + phi1
  lines(phi1, phi2.c2, lty = 2, col = "darkgrey")
  contour(bivn.kde, drawlabels=FALSE, nlevels=k, col=my.cols, add = TRUE)
  
  
  plot(NA, xlim = c(1.96,2), ylim = c(-1.02,-0.97), xlab = expression(phi[1]),
  ylab = expression(phi[2]), cex.lab = 1.5, main = "Bootstrap")
  grid()
  
  
  # Adding boundary of constraint |phi_2| < 1
  abline(h = c(-1,1), lty = 2, col = "darkgrey")
  
  # Adding boundary of constraint phi_2 = 1 - phi_1 
  phi1 = seq(from = -2, to = 2, length.out = 10^3)
  phi2.c1 = 1 - phi1
  lines(phi1, phi2.c1, lty = 2, col = "darkgrey")
  
  # Adding boundary of constraint phi_2 = 1 + phi_1 
  phi1 = seq(from = -2, to = 2, length.out = 10^3)
  phi2.c2 = 1 + phi1
  lines(phi1, phi2.c2, lty = 2, col = "darkgrey")
  contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add = TRUE)
  ```
  
  The cost of this approach is the assumption the model captures the dependency structure that is present within the time series. That is to say, we are able to successfully a regenerate new time series process that follows the appropriate distribution for each sampling phase. However, if we are not confident that the model we have selected is a valid estimation of the truth, then using parametric bootstrap with the estimated model is highly problematic as it does not represent the dependency between observations. To complicate matters further, the traditional bootstrap that  consists of simple random sampling with replacement, in the presence of dependency,  would also be highly inaccurate and downright reckless to use. Therefore, to  preserve the dependency structure of the original data one would have to use _block bootstrapping_, which is a form of non-parametric bootstrap. There  are many different kinds of block bootstraps for time series, which are descendents of the Moving Block Bootstrap (MBB) presented in \@ref(thm:mbb).
  
  ```{theorem, label="mbb", "Moving Block Bootstrap"}
  Suppose that $\left(X_t\right)$ is weakly stationary time series with $N$ 
  observations. 
  
  1. Divide time series into overlapping blocks $\left\{S_1, \ldots, S_M\right\}$ of
  length $\ell$, $1 \le \ell < N$, resulting in $M = N - \ell + 1$ blocks structured as: 
  
  $$\begin{aligned}   
  {S_1}& = & ({X_1}, &{X_2}, \cdots , {X_\ell}) & && \\
  {S_2}& = & &( {X_2}, {X_3}, \cdots , {X_{\ell + 1}}) & && \\
  & \cdots & & {} & \cdots && \\
  {S_M} & = & & & {} &&( {X_{N-\ell+1}}, {X_{N-\ell+2}}, \cdots , {X_{N}})
  \end{aligned}$$
  
  2. Draw $M = \left\lfloor {\frac{N}{\ell}} \right\rfloor$ blocks with replacement
  from these $\left\{S_1, \ldots, S_M\right\}$ blocks and place them in order to form
  $(X_t^*)$. 
  3. Compute the statistic of interest on the simulated 
  sample $(X_t^*)$.
  4. Repeat Steps 2 and 3 $B$ times where $B$ is sufficiently "large" 
  (typically $100 \leq B \leq 10000$).
  5. Compute the empirical mean and variance on the statistic of interest based on 
  the $B$ independent replications. 
  ```
  
  
  The approach taken by MBB ensures that within each block the dependency
  between observations is preserved. Though, one particular issue that now arises is that some inaccuracy is introduced as a result of successive blocks potentially being independent from each other. In reality, this is one of the trade-offs of the MBB approach that can be mitigated by selecting an optimal $\ell$. To lower the inaccuracy of the procedure the selection of $\ell = N^{1/3}$ as $N \to \infty$ is optimal (proof left as an exercise to the reader). An earlier variant of MBB was called nonoverlapping block bootstrap (NBB), which prohbited the sharing of data described in \@ref(thm:nbb).
  
  ```{theorem, label="nbb", "Nonoverlapping Block Bootstrap"}
  Suppose that $\left(X_t\right)$ is weakly stationary time series with $N$ 
  observations. 
  
  1. Divide time series into nonoverlapping blocks $\left\{S_1, \ldots, S_M\right\}$ of
  length $\ell$, $1 \le \ell < N$, resulting in $M = \left\lfloor {\frac{N}{\ell}} \right\rfloor$ blocks structured as: 
  
  $$\begin{aligned}   
  {S_1}& = & ({X_1}, {X_2}, \cdots , {X_\ell})& & && \\
  {S_2}& = & &( {X_{\ell+1}}, {X_{\ell+2}}, \cdots , {X_{2\ell}}) & && \\
  & \cdots & & {} & \cdots && \\
  {S_K} & = & & & {} &&( {X_{\left({}\right)}}, {X_{N-\ell+2}}, \cdots , {X_{N}})
  \end{aligned}$$
  
  2. Draw $M$ blocks with replacement
  from these $\left\{S_1, \ldots, S_M\right\}$ blocks and place them in order to form
  $(X_t^*)$. 
  3. Compute the statistic of interest on the simulated 
  sample $(X_t^*)$.
  4. Repeat Steps 2 and 3 $B$ times where $B$ is sufficiently "large" 
  (typically $100 \leq B \leq 10000$).
  5. Compute the empirical mean and variance on the statistic of interest based on 
  the $B$ independent replications. 
  
  ```
  
  Alternatively, depending on the case one can also use modification of MBB that seeks to change how the beginning and end of the time series is weighted such as a circular block-bootstrap (CBB), that seeks improve dependency by wrapping observations, or a stationary bootstrap (SB), that randomizes block length by under a geometric distribution of mean $\ell$. Regardless, the outcomes from using MBB on time series is considerably better than just resampling, which is also possible if we set $\ell = 1$ since the length of $(X_t^*)$ is $M\times\ell$.
  
  Having said this, the implementation of Theorem \@ref(thm:nbb) follows the
  similar mold of generating data, estimating, and returning a value. The most notably deviation from Theorem \@ref(thm:nbb) is the use of an indexing trick to indicate the start period of each block that allows for only one copy of the data to be held within memory instead of $m$ blocks of length $\ell$ in addition to the initial copy.
  
  ```{r funblockbootstrap}
  ar1_blockboot = function(Xt, block_len = 10, B = 500) {
  
  n = length(Xt)            # Length of Time Series
  res = rep(NA, B)          # Bootstrapped Statistics
  m = floor(n/block_len)    # Amount of Blocks
  
  for (i in seq_len(B)) {   # Begin MMB
  
  set.seed(i + 1199)      # Set seed for reproducibility
  x_star = rep(NA, n)     # Setup storage for new TS
  
  for (j in seq_len(m)) { # Simulate new time series 
  
  index = sample(m, 1)  # Randomize block starting points
  
  # Place block into time series
  x_star[(block_len*(j - 1) + 1):(block_len*j)] = 
  Xt[(block_len*(index - 1) + 1):(block_len*index)]
  }
  
  # Calculate parameter with protection
  res[i] = tryCatch(arima(x_star, order = c(1,0,0), include.mean = FALSE)$coef,
  error = function(e) NA)
  }
  
  na.omit(res)              # Release bootstrap result
  }
  ```
  
  With the block bootstrapping technique in hand, let us consider a scenario where the model's assumption that the residuals will be Gaussian is violated. Consider two AR(1) processes with the same coefficient but different noise generation procedures:
    
    $$
    \begin{aligned}
  \mathcal{M}_1:&{}& X_t &=0.5 X_{t-1} + W_t,&{}& W_t\sim \mathcal{N} (0,1) \\
  \mathcal{M}_2:&{}& X_t &=0.5 X_{t-1} + V_t,&{}& V_t\sim t_4
  \end{aligned}
  $$
    
    The generation procedure for $\mathcal{M}_1$ is straightforward, use: `gen_gts()`. However, the underlying generating mechanism for `gen_gts()` relies upon the noise being Gaussian. Therefore, to generate $\mathcal{M}_2$, one must use _R_'s `arima.sim()` function with a custom noise generator defined to sample from a $t$ distribution.
  
  ```{r mbbdata, cache = TRUE}
  set.seed(1)                                        # Set seed for reproducibilty
  xt_m1 = gen_gts(500, AR1(phi = 0.5, sigma2 = 1))   # Gaussian noise only
  xt_m2 = gts(arima.sim(n = 500, list(ar = 0.5, ma = 0), # Multiple noises supported
  rand.gen = function(n, ...) rt(n, df = 4)))
  ```
  
  ```{r mbbdatavis, echo = F, cache = TRUE, fig.cap="AR(1) processes generated under different noise conditions.", dependson="mbbdata"}
  par(mfrow = c(1,2))
  plot(xt_m1, main = "Model 1")
  plot(xt_m2, main = "Model 2")
  ```
  
  From Figure \@ref(fig:mbbdatavis), the time series look remarkably different as a result of the noise process being altered slightly even though the $\phi$ coefficient was held constant. Principally, the difference is the due related to the variance of the $t$ distribution is defined to be $\frac{\nu}{\nu-2}$ for $\nu > 2$ and $\infty$ for $\nu \le 2$. Therefore, the i.i.d white noise processes generated are $\sigma^2_{\mathcal{M}_1} = 1$ when compared to the alternative $\sigma^2_{\mathcal{M}_2} = 2$. With this being said, both the underlying nonparametric and parametric bootstrap estimation procedures betweenthe models will be the same. The implementation for blocking approach regarding an AR(1) was discussed previously and the parametric approach is as follows:
  
  ```{r ar1bootstrapper, cache = TRUE}
  ar1_paraboot = function(model, B = 10000) {
  est.phi = rep(NA,B)    # Define a storage vector
  
  for(i in seq_len(B)) { # Perform bootstrap
  
  set.seed(B + i)      # Set seed for reproducibility
  
  # Simulate time series underneath the estimated model
  x.star = arima.sim(n = 500, list(ar = model$coef, ma = 0),
  sd = sqrt(model$sigma2))
  
  # Attempt to estimate parameters with recovery
  est.phi[i] = tryCatch(arima(x.star, order = c(1,0,0),
  include.mean = FALSE)$coef, 
  error = function(e) NA)
  
  }
  
  na.omit(est.phi)       # Return estimated phis.
  }
  ```
  
  With both functions in hand, the procedure is regulated to calling them
  on the different sets of data and models. 
  
  ```{r blockbootsim, cache = TRUE}
  B = 10000 
  
  # Model 1
  fit_m1_mle = arima(xt_m1, order = c(1,0,0), include.mean = FALSE)
  para_m1_phi  = ar1_paraboot(fit_m1_mle, B = B)
  block_m1_phi = ar1_blockboot(xt_m1, block_len = 25, B = B)
  
  # Model 2
  fit_m2_mle = arima(xt_m2, order = c(1,0,0), include.mean = FALSE)
  para_m2_phi  = ar1_paraboot(fit_m2_mle, B = B)
  block_m2_phi = ar1_blockboot(xt_m2, block_len = 25, B = B) 
  ```
  
  ```{r blockbootmodels, echo = FALSE, cache = T, warning = FALSE, fig.height= 6, fig.width= 9, fig.cap="Estimated parametric and non-parametric block  bootstrap distributions of $\\hat{\\phi}$ for MLE parameter estimates. The histogram bars represent the empirical results from the bootstraps with the green representing parametric bootstrap and the red represent the block bootstrap approach. The dashed vertical line represents the true value of $\\phi$ and the vertical ticks correspond to the limits of the 95% confidence intervals for both estimation techniques."}
  
  ## Rewrite as a function later... 
  # Has a dependency on previously written code.
  
  par(mfrow = c(1,2))
  
  ## Model 1
  hist.phi.ML = hist(para_m1_phi, plot = FALSE, breaks = 25)
  hist.phi.ML$counts = hist.phi.ML$counts/sum(hist.phi.ML$counts)/diff(hist.phi.ML$mids)[1]
  
  hist.phi.ML.bboot = hist(block_m1_phi, plot = FALSE, breaks = 25)
  hist.phi.ML.bboot$counts = hist.phi.ML.bboot$counts/sum(hist.phi.ML.bboot$counts)/diff(hist.phi.ML.bboot$mids)[1]
  
  plot(NA, xlim = c(0.28,0.6), ylim = range(yy.ML), 
  xlab = expression(phi), ylab = "Density", main = "Model 1")
  grid()
  abline(v = 1, lwd = 2, col = "grey60")
  abline(v = 0.95, lwd = 2, lty = 2, col = "grey60")
  plot(hist.phi.ML, xlim = c(0.28,0.6), ylim = range(yy.ML), 
  add = T, col = colors[3], border = "grey60")
  
  plot(hist.phi.ML.bboot, xlim = c(0.28,0.6), ylim = range(yy.ML), 
  add = T, col = colors[1], border = "grey60")
  
  points(quantile(para_m1_phi, 0.975), 0, cex = 3, col = colors[3], pch = "|")
  points(quantile(para_m1_phi, 0.025), 0, cex = 3, col = colors[3], pch = "|")
  
  points(quantile(block_m1_phi, 0.975), 0, cex = 3, col = colors[1], pch = "|")
  points(quantile(block_m1_phi, 0.025), 0, cex = 3, col = colors[1], pch = "|")
  abline(v = 0.5, lwd = 2, lty = 2, col = "grey60")
  
  ## Model 2
  hist.phi.ML = hist(para_m2_phi, plot = FALSE, breaks = 25)
  hist.phi.ML$counts = hist.phi.ML$counts/sum(hist.phi.ML$counts)/diff(hist.phi.ML$mids)[1]
  
  hist.phi.ML.bboot = hist(block_m2_phi, plot = FALSE, breaks = 25)
  hist.phi.ML.bboot$counts = hist.phi.ML.bboot$counts/sum(hist.phi.ML.bboot$counts)/diff(hist.phi.ML.bboot$mids)[1]
  
  plot(NA, xlim = c(0.28,0.6), ylim = range(yy.ML), 
  xlab = expression(phi), ylab = "Density", main = "Model 2")
  grid()
  abline(v = 1, lwd = 2, col = "grey60")
  abline(v = 0.95, lwd = 2, lty = 2, col = "grey60")
  plot(hist.phi.ML, xlim = c(0.28,0.6), ylim = range(yy.ML), 
  add = T, col = colors[3], border = "grey60")
  
  plot(hist.phi.ML.bboot, xlim = c(0.28,0.6), ylim = range(yy.ML), 
  add = T, col = colors[1], border = "grey60")
  
  points(quantile(para_m2_phi, 0.975), 0, cex = 3, col = colors[3], pch = "|")
  points(quantile(para_m2_phi, 0.025), 0, cex = 3, col = colors[3], pch = "|")
  
  points(quantile(block_m2_phi, 0.975), 0, cex = 3, col = colors[1], pch = "|")
  points(quantile(block_m2_phi, 0.025), 0, cex = 3, col = colors[1], pch = "|")
  
  abline(v = 0.5, lwd = 2, lty = 2, col = "grey60")
  
  ```
  
  The results from the parametric and nonparametric block bootstrapping techniques
  are displayed in Figure \@ref(fig:blockbootmodels). In both instances, we can 
  see that the block bootstrap had a narrower distribution. However, there was a
  considerable bias that lead the distribution to be off-centered. The origins of
  bias come from the independence associated between blocks as they are placed
  into $(X^*_t)$. 