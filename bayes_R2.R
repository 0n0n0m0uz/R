#' ---
#' title: "Bayesian R2 and LOO-R2"
#' author: "Aki Vehtari, Andrew Gelman, Ben Goodrich, Jonah Gabry"
#' date: "`r format(Sys.Date())`"
#' output:
#'   html_document:
#'     fig_caption: yes
#'     toc: TRUE
#'     toc_depth: 2
#'     number_sections: TRUE
#'     toc_float:
#'       smooth_scroll: FALSE
#' bibliography: bayes_R2.bib
#' csl: harvard-cite-them-right.csl
#' ---

#+ knitr-setup, include=FALSE
knitr::opts_chunk$set(
  comment=NA, 
  cache=FALSE,
  fig.width = 10
)

#' # Introduction
#' 
#' This notebook demonstrates Bayesian posterior distributions of
#' model based R-squared and LOO-R2 measures for linear and logistic
#' regression. This notebook is a supplement for the paper
#'
#' - Andrew Gelman, Ben Goodrich, Jonah Gabry, and Aki Vehtari (2018). R-squared for Bayesian regression models. The American Statistician, doi:10.1080/00031305.2018.1549100. [Online](https://doi.org/10.1080/00031305.2018.1549100) [Preprint](http://www.stat.columbia.edu/~gelman/research/unpublished/bayes_R2_v3.pdf).
#'
#' We specifically define R-squared as
#' $$
#' R^2 = \frac{{\mathrm Var}_{\mu}}{{\mathrm Var}_{\mu}+{\mathrm Var}_{\mathrm res}},
#' $$
#' where ${\mathrm Var}_{\mu}$ is variance of modelled predictive means,
#' and ${\mathrm Var}_{\mathrm res}$ is the modelled residual variance.
#' Specifically both of these are computed only using posterior
#' quantities from the fitted model.
#'
#' The residual based R2 uses draws from the residual distribution and
#' is defined as
#' $$
#' {\mathrm Var}_{\mu}^s = V_{n=1}^N \hat{y}_n^s\\
#' {\mathrm Var}_{\mathrm res}^s = V_{n=1}^N \hat{e}_n^s,
#' $$
#' where $\hat{e}_n^s=y_n-\hat{y}_n^s$.
#' 
#' The model based R2 uses draws from the modeled residual variances.
#' For linear regression we define
#' $$
#' {\mathrm Var}_{\mathrm res}^s = (\sigma^2)^s,
#' $$
#' and for logistic regression, following @Tjur2009, we define
#' $$
#' {\mathrm Var}_{\mathrm res}^s = \frac{1}{N}\sum_{n=1}^N (\pi_n^s(1-\pi_n^s)),
#' $$
#' where $\pi_n^s$ are predicted probabilities.
#' 
#' The LOO-R2 uses LOO residuals and is defined
#' $1-{\mathrm Var}_{\mathrm loo-res}/{\mathrm Var}_{y}$ and
#' $$
#' {\mathrm Var}_{y} = V_{n=1}^N {y}_n\\
#' {\mathrm Var}_{\mathrm loo-res} = V_{n=1}^N \hat{e}_{{\mathrm loo},n},
#' $$
#' where $\hat{e}_{{\mathrm loo},n}=y_n-\hat{y}_{{\mathrm loo},n}$.
#' The uncertainty in LOO-R2 comes from not knowing the future data
#' distribution [@Vehtari+Ojanen:2012]. We can approximate draws from
#' the corresponding distribution using Bayesian bootstrap
#' [@Vehtari+Lampinen:2002b].
#' 

#' File contents:</br>
#' 
#'   - bayes_R2_res function (using residuals)
#'   - bayes_R2 function     (using sigma and latent space for binary)
#'   - loo_R2 function       (using LOO residuals and Bayesian bootstrap for uncertainty)
#'   - code for fitting the models for the example in the paper
#'   - code for producing the plots in the paper (both base graphics and ggplot2 code)
#'   - several examples
#'
#' [Corresponding Bayes_R2.R file](bayes_R2.R)
#' 
#' -------------
#' 

#+ include=FALSE
# switch this to TRUE to save figures in separate files
savefigs <- FALSE

#' **Load packages**
#+ setup, message=FALSE, error=FALSE, warning=FALSE
library("rprojroot")
root<-has_dirname("RAOS-Examples")$make_fix_file()
library("rstanarm")
options(mc.cores = parallel::detectCores())
library("loo")
library("ggplot2")
library("bayesplot")
theme_set(bayesplot::theme_default(base_family = "sans"))
library("latex2exp")
library("foreign")
library("bayesboot")
SEED <- 1800
set.seed(SEED)

#' # Functions for Bayesian R-squared for stan_glm models. 
#' 
#' A more complicated version of this function that is compatible with
#' stan_glmer models and specifying newdata will be available in the
#' rstanarm package.
#'
#' \@param fit A fitted model object returned by stan_glm.</br>
#' \@return A vector of R-squared values with length equal to</br>
#'      the number of posterior draws.
#'
#' -------------
#' 

#' **Bayes-R2 function using residuals**
bayes_R2_res <- function(fit) {
  y <- rstanarm::get_y(fit)
  ypred <- rstanarm::posterior_linpred(fit, transform = TRUE)
  if (family(fit)$family == "binomial" && NCOL(y) == 2) {
    trials <- rowSums(y)
    y <- y[, 1]
    ypred <- ypred %*% diag(trials)
  }
  e <- -1 * sweep(ypred, 2, y)
  var_ypred <- apply(ypred, 1, var)
  var_e <- apply(e, 1, var)
  var_ypred / (var_ypred + var_e)
}

#' **Bayes-R2 function using modelled (approximate) residual variance**
bayes_R2 <- function(fit) {
  mupred <- rstanarm::posterior_linpred(fit, transform = TRUE)
  var_mupred <- apply(mupred, 1, var)
  if (family(fit)$family == "binomial" && NCOL(y) == 1) {
      sigma2 <- apply(mupred*(1-mupred), 1, mean)
  } else {
      sigma2 <- as.matrix(fit, pars = c("sigma"))^2
  }
  var_mupred / (var_mupred + sigma2)
}

#' **LOO-R2 function using LOO-residuals and Bayesian bootstrap**
loo_R2 <- function(fit) {
    y <- get_y(fit)
    ypred <- posterior_linpred(fit, transform = TRUE)
    ll <- log_lik(fit)
    M <- length(fit$stanfit@sim$n_save)
    N <- fit$stanfit@sim$n_save[[1]]
    r_eff <- relative_eff(exp(ll), chain_id = rep(1:M, each = N))
    psis_object <- psis(log_ratios = -ll, r_eff = r_eff)
    ypredloo <- E_loo(ypred, psis_object, log_ratios = -ll)$value
    eloo <- ypredloo-y
    n <- length(y)
    rd<-rudirichlet(4000, n)
    vary <- (rowSums(sweep(rd, 2, y^2, FUN = "*")) -
             rowSums(sweep(rd, 2, y, FUN = "*"))^2)*(n/(n-1))
    vareloo <- (rowSums(sweep(rd, 2, eloo^2, FUN = "*")) -
                rowSums(sweep(rd, 2, eloo, FUN = "*")^2))*(n/(n-1))
    looR2 <- 1-vareloo/vary
    looR2[looR2 < -1] <- -1
    looR2[looR2 > 1] <- 1
    return(looR2)
}

#' # Experiments

#' ## Toy data with n=5
x <- 1:5 - 3
y <- c(1.7, 2.6, 2.5, 4.4, 3.8) - 3
xy <- data.frame(x,y)

#' **Lsq fit**
fit <- lm(y ~ x, data = xy)
ols_coef <- coef(fit)
yhat <- ols_coef[1] + ols_coef[2] * x
r <- y - yhat
rsq_1 <- var(yhat)/(var(y))
rsq_2 <- var(yhat)/(var(yhat) + var(r))
round(c(rsq_1, rsq_2), 3)

#' **Bayes fit**
#+ results='hide'
fit_bayes <- stan_glm(y ~ x, data = xy,
  prior_intercept = normal(0, 0.2, autoscale = FALSE),
  prior = normal(1, 0.2, autoscale = FALSE),
  prior_aux = NULL,
  seed = SEED
)
posterior <- as.matrix(fit_bayes, pars = c("(Intercept)", "x"))
post_means <- colMeans(posterior)

#' **Median Bayesian R2**
round(median(bayesR2res<-bayes_R2_res(fit_bayes)), 2)
round(median(bayesR2<-bayes_R2(fit_bayes)), 2)

#' **LOO R2**
looR2<-loo_R2(fit_bayes)
round(median(looR2), 2)

#' With small n LOO-R2 is much lower than median of Bayesian R2
#' 

#' **Figures**
#'
#' The first section of code below creates plots using base R graphics. </br>
#' Below that there is code to produce the plots using ggplot2.

# take a sample of 20 posterior draws
keep <- sample(nrow(posterior), 20)
samp_20_draws <- posterior[keep, ]

#' **Base graphics version**

#+ eval=FALSE, include=FALSE
if (savefigs) pdf("fig/rsquared1a.pdf", height=4, width=5)
#+
par(mar=c(3,3,1,1), mgp=c(1.7,.5,0), tck=-.01)
plot(
  x, y,
  ylim = range(x),
  xlab = "x",
  ylab = "y",
  main = "Least squares and Bayes fits",
  bty = "l",
  pch = 20
)
abline(coef(fit)[1], coef(fit)[2], col = "black")
text(-1.6,-.7, "Least-squares\nfit", cex = .9)
abline(0, 1, col = "blue", lty = 2)
text(-1, -1.8, "(Prior regression line)", col = "blue", cex = .9)
abline(coef(fit_bayes)[1], coef(fit_bayes)[2], col = "blue")
text(1.4, 1.2, "Posterior mean fit", col = "blue", cex = .9)
points(
  x,
  coef(fit_bayes)[1] + coef(fit_bayes)[2] * x,
  pch = 20,
  col = "blue"
)
#+ eval=FALSE, include=FALSE
if (savefigs) dev.off()

#+ eval=FALSE, include=FALSE
if (savefigs) pdf("fig/rsquared1b.pdf", height=4, width=5)
#+
par(mar=c(3,3,1,1), mgp=c(1.7,.5,0), tck=-.01)
plot(
  x, y,
  ylim = range(x),
  xlab = "x",
  ylab = "y",
  bty = "l",
  pch = 20,
  main = "Bayes posterior simulations"
)
for (s in 1:nrow(samp_20_draws)) {
  abline(samp_20_draws[s, 1], samp_20_draws[s, 2], col = "#9497eb")
}
abline(
  coef(fit_bayes)[1],
  coef(fit_bayes)[2],
  col = "#1c35c4",
  lwd = 2
)
points(x, y, pch = 20, col = "black")
#+ eval=FALSE, include=FALSE
if (savefigs) dev.off()

#' **ggplot version**
theme_update(
  plot.title = element_text(face = "bold", hjust = 0.5), 
  axis.text = element_text(size = rel(1.1))
)
fig_1a <-
  ggplot(xy, aes(x, y)) +
  geom_point() +
  geom_abline( # ols regression line
    intercept = ols_coef[1],
    slope = ols_coef[2],
    size = 0.4
  ) +
  geom_abline( # prior regression line
    intercept = 0,
    slope = 1,
    color = "blue",
    linetype = 2,
    size = 0.3
  ) +
  geom_abline( # posterior mean regression line
    intercept = post_means[1],
    slope = post_means[2],
    color = "blue",
    size = 0.4
  ) +
  geom_point(
    aes(y = post_means[1] + post_means[2] * x),
    color = "blue"
  ) +
  annotate(
    geom = "text",
    x = c(-1.6, -1, 1.4),
    y = c(-0.7, -1.8, 1.2),
    label = c(
      "Least-squares\nfit",
      "(Prior regression line)",
      "Posterior mean fit"
    ),
    color = c("black", "blue", "blue"),
    size = 3.8
  ) +
  ylim(range(x)) + 
  ggtitle("Least squares and Bayes fits")
plot(fig_1a)
#+ eval=FALSE, include=FALSE
ggsave("fig/rsquared1a-gg.pdf", width = 5, height = 4)
#+
fig_1b <-
  ggplot(xy, aes(x, y)) +
  geom_abline( # 20 posterior draws of the regression line
    intercept = samp_20_draws[, 1],
    slope = samp_20_draws[, 2],
    color = "#9497eb",
    size = 0.25
  ) +
  geom_abline( # posterior mean regression line
    intercept = post_means[1],
    slope = post_means[2],
    color = "#1c35c4",
    size = 1
  ) +
  geom_point() +
  ylim(range(x)) + 
  ggtitle("Bayes posterior simulations")
plot(fig_1b)
#+ eval=FALSE, include=FALSE
ggsave("fig/rsquared1b-gg.pdf", width = 5, height = 4)

#' **Bayesian R^2 Posterior and median**
mcmc_hist(data.frame(bayesR2), binwidth=0.02)+xlim(c(0,1))+
  xlab(TeX('$R^2$'))+
  geom_vline(xintercept=median(bayesR2))+
  ggtitle('Bayesian R squared posterior and median')
#+ eval=FALSE, include=FALSE
ggsave("fig/bayesr2post.pdf", width = 5, height = 4)

#' **Compare draws from residual variances and sigma**
fit<-fit_bayes
y <- rstanarm::get_y(fit)
ypred <- rstanarm::posterior_linpred(fit, transform = TRUE)
e <- -1 * sweep(ypred, 2, y)
var_e <- apply(e, 1, var)
sigma <- as.matrix(fit, pars = c("sigma"))
qplot(sqrt(var_e),sigma)+geom_abline()+
    ggtitle('Toy data with n=5')+
    xlab('Draws from residual sd')+
    ylab('Draws from sigma')

#' Draws from residual sd are truncated and not properly using the
#' model information, which is reflected in R2 distribution below.
qplot(bayes_R2_res(fit_bayes),bayes_R2(fit_bayes))+geom_abline()+
    ggtitle('Toy data with n=5')+
    xlab('Bayesian R2 with draws from residual variance')+
    ylab('Bayesian R2 with draws from sigma')

#' **Compare residual and model based Bayesian R2, and Bayesian LOO-R2**
#+ message=FALSE, error=FALSE, warning=FALSE
pxl<-xlim(-0.5,1)
p1<-mcmc_hist(data.frame(bayesR2res), binwidth=0.02)+pxl+
    ggtitle('Toy data with n=5')+
    xlab('Bayesian R2 with draws from residuals')+
    geom_vline(xintercept=median(bayesR2res))
p2<-mcmc_hist(data.frame(bayesR2), binwidth=0.02)+pxl+
    xlab('Bayesian R2 with draws from sigma')+
    geom_vline(xintercept=median(bayesR2))
p3<-mcmc_hist(data.frame(looR2), binwidth=0.02)+pxl+
  xlab('Bayesian LOO-R2')+
  geom_vline(xintercept=median(looR2))
bayesplot_grid(p1,p2,p3)

#' With small n, draws from the residual variance are not a good
#' representation of the uncertainty.
#' 

#' ## Toy logistic regression example, n=20
set.seed(20)
y<-rbinom(n=20,size=1,prob=(1:20-0.5)/20)
data <- data.frame(rvote=y, income=1:20)
#+ results='hide'
fit_logit <- stan_glm(rvote ~ income, family=binomial(link="logit"), data=data)

#' **Median Bayesian R2**
round(median(bayesR2res<-bayes_R2_res(fit_logit)), 2)
round(median(bayesR2<-bayes_R2(fit_logit)), 2)

#' **LOO R2**
looR2<-loo_R2(fit_logit)
round(median(looR2), 2)

#' With small n LOO-R2 is much lower than median of Bayesian R2
#'

#' **Compare residual and model based Bayesian R2, and Bayesian LOO-R2**
#+ message=FALSE, error=FALSE, warning=FALSE
pxl<-xlim(-.5,1)
p1<-mcmc_hist(data.frame(bayesR2res), binwidth=0.05)+pxl+
    ggtitle('Toy logistic data with n=20')+
    xlab('Bayesian R2 with draws from residuals')+
    geom_vline(xintercept=median(bayesR2res))
p2<-mcmc_hist(data.frame(bayesR2), binwidth=0.05)+pxl+
    xlab('Bayesian R2 with draws from sigma')+
    geom_vline(xintercept=median(bayesR2))
p3<-mcmc_hist(data.frame(looR2), binwidth=0.05)+pxl+
  xlab('Bayesian LOO-R2')+
  geom_vline(xintercept=median(looR2))
bayesplot_grid(p1,p2,p3)

#' With small n, draws from the residual variance are not a good
#' representation of the uncertainty.
#' 

#' ## Mesquite - linear regression
#' 
#' Predicting the yields of mesquite bushes, n=46
#' 

#' **Load data**
mesquite <- read.table(root("Mesquite/data","mesquite.dat"), header=TRUE)
mesquite$canopy_volume <- mesquite$diam1 * mesquite$diam2 * mesquite$canopy_height
mesquite$canopy_area <- mesquite$diam1 * mesquite$diam2
mesquite$canopy_shape <- mesquite$diam1 / mesquite$diam2
(n <- nrow(mesquite))

#' **Predict log weight model with log canopy volume, log canopy shape, and group**
#+ results='hide'
fit_5 <- stan_glm(log(weight) ~ log(canopy_volume) + log(canopy_shape) +
    group, data=mesquite)

#' **Median Bayesian R2**
round(median(bayesR2res<-bayes_R2_res(fit_5)), 2)
round(median(bayesR2<-bayes_R2(fit_5)), 2)

#' **LOO R2**
looR2<-loo_R2(fit_5)
round(median(looR2), 2)

#' LOO-R2 is slightly lower than median of Bayesian R2
#' 

#' **Compare residual and model based Bayesian R2**
#+ message=FALSE, error=FALSE, warning=FALSE
pxl<-xlim(0.6, 0.95)
p1<-mcmc_hist(data.frame(bayesR2res), binwidth=0.01)+pxl+
    ggtitle('Mesquite data with n=46')+
    xlab('Bayesian R2 with draws from residuals')+
    geom_vline(xintercept=median(bayesR2res))
p2<-mcmc_hist(data.frame(bayesR2), binwidth=0.01)+pxl+
    xlab('Bayesian R2 with draws from sigma')+
    geom_vline(xintercept=median(bayesR2))
p3<-mcmc_hist(data.frame(looR2), binwidth=0.01)+pxl+
  xlab('Bayesian LOO-R2')+
  geom_vline(xintercept=median(looR2))
bayesplot_grid(p1,p2,p3)

#' With small n, draws from the residual variance are not a good
#' representation of the uncertainty and R2 distribution is more
#' narrow (50%) than when using draws from the model sigma.
#' 

#' ## LowBwt -- logistic regression
#' 
#' Predict low birth weight, n=189, from @HosmerLemeshow2000 </br>
#' This data was also used by @Tjur2009
#' 

#' **Load data**
lowbwt <- read.table(root("LowBwt/data","lowbwt.dat"), header=TRUE)
lowbwt$race <- factor(lowbwt$race)
(n <- nrow(lowbwt))

#' **Predict low birth weight**
#+ results='hide'
fit <- stan_glm(low ~ age + lwt + race + smoke,
                family=binomial(link="logit"), data=lowbwt)

#' **Median Bayesian R2**
round(median(bayesR2res<-bayes_R2_res(fit)), 2)
round(median(bayesR2<-bayes_R2(fit)), 2)

#' **LOO R2**
looR2<-loo_R2(fit)
round(median(looR2), 2)

#' LOO-R2 is much lower than median of Bayesian R2
#' 

#' **Compare residual and model based Bayesian R2**
#+ message=FALSE, error=FALSE, warning=FALSE
pxl<-xlim(-0.15, 0.3)
p1<-mcmc_hist(data.frame(bayesR2res), binwidth=0.01)+pxl+
    ggtitle('Low birth weight data with n=189')+
    xlab('Bayesian R2 with draws from residuals')+
    geom_vline(xintercept=median(bayesR2res))
p2<-mcmc_hist(data.frame(bayesR2), binwidth=0.01)+pxl+
    xlab('Bayesian R2 with draws from sigma')+
    geom_vline(xintercept=median(bayesR2))
p3<-mcmc_hist(data.frame(looR2), binwidth=0.01)+pxl+
  xlab('Bayesian LOO-R2')+
  geom_vline(xintercept=median(looR2))
bayesplot_grid(p1,p2,p3)

#' With moderate n, the distribution obtained using draws from the
#' residual variances is still more narrow than the distribution
#' obtained using draws from the model based variance.
#' 

#' ## LowBwt -- linear regression
#' 
#' Predict birth weight, n=189, from @HosmerLemeshow2000 </br>
#' @Tjur2009 used logistic regression for dichotomized birth weight.</br>
#' Below we use the continuos valued birth weight.
#' 

#' **Predict birth weight**
#+ results='hide'
fit <- stan_glm(bwt ~ age + lwt + race + smoke, data=lowbwt)

#' **Median Bayesian R2**
round(median(bayesR2res<-bayes_R2_res(fit)), 2)
round(median(bayesR2<-bayes_R2(fit)), 2)

#' **LOO R2**
looR2<-loo_R2(fit)
round(median(looR2), 2)

#' LOO-R2 is much lower than median of Bayesian R2
#' 

#' **Compare residual and model based Bayesian R2**
#+ message=FALSE, error=FALSE, warning=FALSE
pxl<-xlim(-0.12, 0.36)
p1<-mcmc_hist(data.frame(bayesR2res), binwidth=0.01)+pxl+
    ggtitle('Birth weight data with n=189')+
    xlab('Bayesian R2 with draws from residuals')+
    geom_vline(xintercept=median(bayesR2res))
p2<-mcmc_hist(data.frame(bayesR2), binwidth=0.01)+pxl+
    xlab('Bayesian R2 with draws from sigma')+
    geom_vline(xintercept=median(bayesR2))
p3<-mcmc_hist(data.frame(looR2), binwidth=0.01)+pxl+
  xlab('Bayesian LOO-R2')+
  geom_vline(xintercept=median(looR2))
bayesplot_grid(p1,p2,p3)

#' With moderate n, the distribution obtained using draws from the
#' residual variances is still more narrow than the distribution
#' obtained using draws from the model based variance.
#' 

#' ## KidIQ - linear regression
#' 
#' Children's test scores data, n=434
#' 

#' **Load children's test scores data**

kidiq <- read.dta(file=root("KidIQ/data","kidiq.dta"))
(n <- nrow(kidiq))

#' **Predict test score**
#+ results='hide'
fit_3 <- stan_glm(kid_score ~ mom_hs + mom_iq, data=kidiq, seed=1507)

#' **Median Bayesian R2**
round(median(bayesR2res<-bayes_R2_res(fit_3)), 2)
round(median(bayesR2<-bayes_R2(fit_3)), 2)

#' **LOO R2**
looR2<-loo_R2(fit_3)
round(median(looR2), 2)

#' LOO-R2 is slightly lower than median of Bayesian R2
#' 

#' **Compare residual and model based Bayesian R2**
#+ message=FALSE, error=FALSE, warning=FALSE
pxl<-xlim(0.05, 0.35)
p1<-mcmc_hist(data.frame(bayesR2res), binwidth=0.01)+pxl+
    ggtitle('KidIQ data with n=434')+
    xlab('Bayesian R2 with draws from residuals')+
    geom_vline(xintercept=median(bayesR2res))
p2<-mcmc_hist(data.frame(bayesR2), binwidth=0.01)+pxl+
    xlab('Bayesian R2 with draws from sigma')+
    geom_vline(xintercept=median(bayesR2))
p3<-mcmc_hist(data.frame(looR2), binwidth=0.01)+pxl+
  xlab('Bayesian LOO-R2')+
  geom_vline(xintercept=median(looR2))
bayesplot_grid(p1,p2,p3)

#' ## Earnings - logistic and linear regression
#' 
#' Predict respondents' yearly earnings using survey data from 1990.</br>
#' logistic regression n=1374, linear regression n=1187
#' 

#' **Load data**
earnings_all <- read.csv(root("Earnings/data","earnings.csv")) 
earnings_all$positive <- earnings_all$earn > 0
(n_all <- nrow(earnings_all))
# only non-zero earnings
earnings <- earnings_all[earnings_all$positive, ]
(n <- nrow(earnings))
earnings$log_earn <- log(earnings$earn)

#' **Bayesian logistic regression on non-zero earnings**</br>
#' Null model
#+ results='hide'
fit_0 <- stan_glm(positive ~ 1,
                   family = binomial(link = "logit"),
                  data = earnings_all)
#+
loo0<-loo(fit_0)
#' Predict using height and sex
#+ results='hide'
fit_1a <- stan_glm(positive ~ height + male,
                   family = binomial(link = "logit"),
                   data = earnings_all)
#+
loo1a<-loo(fit_1a)
compare_models(loo0, loo1a)

#' There is a clear difference in predictive performance.
#' 

#' **Median Bayesian R2**
round(median(bayesR2res<-bayes_R2_res(fit_1a)), 3)
round(median(bayesR2<-bayes_R2(fit_1a)), 3)

#' **LOO R2**
looR2<-loo_R2(fit_1a)
round(median(looR2), 2)

#' LOO-R2 is slightly lower than median of Bayesian R2
#' 

#' **Compare residual and model based Bayesian R2**
#+ message=FALSE, error=FALSE, warning=FALSE
pxl<-xlim(0.02, 0.11)
p1<-mcmc_hist(data.frame(bayesR2res), binwidth=0.002)+pxl+
    ggtitle('Earnings data with n=1374')+
    xlab('Bayesian R2 with draws from residuals')+
    geom_vline(xintercept=median(bayesR2res))
p2<-mcmc_hist(data.frame(bayesR2), binwidth=0.002)+pxl+
    xlab('Bayesian R2 with draws from sigma')+
    geom_vline(xintercept=median(bayesR2))
p3<-mcmc_hist(data.frame(looR2), binwidth=0.002)+pxl+
  xlab('Bayesian LOO-R2')+
  geom_vline(xintercept=median(looR2))
bayesplot_grid(p1,p2,p3)

#' With plenty of data, there is not much difference between using
#' draws from residuals or from sigma.
#' 

#' **Bayesian probit regression on non-zero earnings**</br>
#+ results='hide'
fit_1p <- stan_glm(positive ~ height + male,
                   family = binomial(link = "probit"),
                   data = earnings_all)
#+
loo1p<-loo(fit_1p)
compare_models(loo1a, loo1p)

#' There is no practical difference in predictive performance between
#' logit and probit.
#' 
    
#' **Median Bayesian R2**
round(median(bayesR2), 3)
round(median(bayesR2p<-bayes_R2(fit_1p)), 3)

#' **LOO R2**
looR2p<-loo_R2(fit_1p)
round(median(looR2), 2)
round(median(looR2p), 2)

#' LOO-R2 is slightly lower than median of Bayesian R2
#' 

#' **Compare logistic and probit models using new Bayesian R2**
#+ message=FALSE, error=FALSE, warning=FALSE
pxl<-xlim(0.02, 0.11)
p1<-mcmc_hist(data.frame(bayesR2), binwidth=0.002)+pxl+
    ggtitle('Earnings data with n=1374')+
    xlab('Bayesian R2 for logistic model')+
    geom_vline(xintercept=median(bayesR2))
p2<-mcmc_hist(data.frame(bayesR2p), binwidth=0.002)+pxl+
    xlab('Bayesian R2 for probit model')+
    geom_vline(xintercept=median(bayesR2p))
p3<-mcmc_hist(data.frame(looR2), binwidth=0.002)+pxl+
  xlab('Bayesian LOO-R2 for logistic model')+
  geom_vline(xintercept=median(looR2))
p4<-mcmc_hist(data.frame(looR2p), binwidth=0.002)+pxl+
  xlab('Bayesian LOO-R2 for probit model')+
  geom_vline(xintercept=median(looR2p))
bayesplot_grid(p1,p3,p2,p4)

#' There is no practical difference in predictive performance between
#' logit and probit.
#' 

#' **Bayesian model on positive earnings on log scale**
#+ results='hide'
fit_1b <- stan_glm(log_earn ~ height + male, data = earnings)

#' **Median Bayesian R2**
round(median(bayesR2res<-bayes_R2_res(fit_1b)), 3)
round(median(bayesR2<-bayes_R2(fit_1b)), 3)

#' **LOO R2**
looR2<-loo_R2(fit_1b)
round(median(looR2), 2)

#' LOO-R2 is slightly lower than median of Bayesian R2
#' 

#' **Compare residual and model based Bayesian R2**
#+ message=FALSE, error=FALSE, warning=FALSE
pxl<-xlim(0.02, 0.15)
p1<-mcmc_hist(data.frame(bayesR2res), binwidth=0.002)+pxl+
    ggtitle('Positive earnings data with n=1187')+
    xlab('Bayesian R2 with draws from residuals')+
    geom_vline(xintercept=median(bayesR2res))
p2<-mcmc_hist(data.frame(bayesR2), binwidth=0.002)+pxl+
    xlab('Bayesian R2 with draws from sigma')+
    geom_vline(xintercept=median(bayesR2))
p3<-mcmc_hist(data.frame(looR2), binwidth=0.002)+pxl+
  xlab('Bayesian LOO-R2')+
  geom_vline(xintercept=median(looR2))
bayesplot_grid(p1,p2,p3)

#' With plenty of data, there is not much difference between using
#' draws from residuals or from sigma.
#' 

#' # References {.unnumbered}
#'
#' <div id="refs"></div>
#' 

