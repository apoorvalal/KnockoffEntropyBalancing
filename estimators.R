# %% ####################################################
library(LalRUtils)
LalRUtils::libreq(data.table, hdm, ebal, glue,
                  lfe, janitor, caret, kosel)

# %% prepare X matrix only
prep_matrix = function(
    dummies, # dummy vars
    continuouses, # continuous vars to modify fn form
    dat, # must be data.table,
    corr_cut = 0.98, # cutoff for correlation threshold to drop one of the vars
    k = 2, # default pairwise interactions
    m = 2  # default quadratic functions
  ) {
  controls = c(dummies, continuouses)
  #############################################
  # functional form changes for continuous vars
  #############################################
  polynomial_dfs = list()
  # log(x+1) for all continuous columns
  polynomial_dfs[[1]] = log1p(dat[, ..continuouses])
  names(polynomial_dfs[[1]]) = paste0("log_", continuouses)
  # construct polynomials
  for (i in 2:m){
    polynomial_dfs[[i]] = dat[, ..continuouses]^i
    names(polynomial_dfs[[i]]) = paste0(continuouses, glue::glue("_{i}"))
  }
  polynomials = list.cbind(polynomial_dfs) %>% as.matrix()
  data = cbind(dat[, ..controls], polynomials)
  # all n-way interactions with polynomials
  X = model.matrix(as.formula(glue::glue("~.^{k} - 1")), data)
  # drop non-varying Xs
  varyvar = apply(X, 2, function(col) nunique(col) > 1)
  X = X[, varyvar]
  # X <- X[, -nearZeroVar(X)]
  # drop highly correlated Xs - this drops polynomials of binary variables
  corm = cor(X)
  hc = findCorrelation(corm, cutoff=corr_cut) # put any value as a "cutoff"
  hc = sort(hc)
  X = X[,-c(hc)]
  return(X)
}


# %% # data prep function for all x, y, w vars
prep_matrices = function(
    outcome, treatment,  # characters
    controls,  # character vector
    dat, # must be data.table,
    corr_cut = 0.95, # cutoff for correlation threshold to drop one of the vars
    k = 2, # default pairwise interactions
    m = 2,  # default quadratic functions
    log_cols = F
  ) {
  polynomial_dfs = list()
  if (log_cols){ # log(x+1) for all columns
    # subset to nonnegative columns
    log_allowed = dat[, lapply(.SD, function(x) all(x >= 0)), .SDcols = controls] %>%
      as.logical
    logcols = controls[log_allowed]
    polynomial_dfs[[1]] = log1p(dat[, ..logcols])
    names(polynomial_dfs[[1]]) = paste0("log_", logcols)
    j = 2
  } else{
    j = 1
  }
  # construct polynomials
  for (i in j:m){
    polynomial_dfs[[i]] = dat[, ..controls]^i
    names(polynomial_dfs[[i]]) = paste0(controls, glue::glue("_{i}"))
  }
  polynomials = list.cbind(polynomial_dfs) %>% as.matrix()
  data = cbind(dat[, ..controls], polynomials)
  # all n-way interactions with polynomials
  X = model.matrix(as.formula(glue::glue("~.^{k} - 1")), data)
  # drop non-varying Xs
  X <- X[, -nearZeroVar(X)]
  # drop highly correlated Xs - this drops polynomials of binary variables
  corm = cor(X)
  hc = findCorrelation(corm, cutoff=corr_cut) # put any value as a "cutoff"
  hc = sort(hc)
  X = X[,-c(hc)]
  # center X matrix
  X = scale(X, center = TRUE, scale = T)
  # drop missings
  df = cbind(dat[[outcome]], dat[[treatment]], X) %>% na.omit
  df = df[order(df[, 2]), ] # move untreated up top
  y = df[, 1]; w = df[, 2]; X = df[, -(1:2)]
  return(list(y, w, X))
}
# %% vanilla difference in means (no controls)

############################################################
# vanilla estimation
############################################################

diffmeans = function(df, y, w){
  m0 = felm(as.formula(glue("{y} ~ {w}")), data = df) %>% robustify
  ols_est0 = summary(m0)$coefficients[2,1:2]
  return(ols_est0)
}

# %% covariate adjustment (lin correction)
covaradjust = function(df, y, w, ctrls){
  fmla = as.formula(glue('{y} ~ {w} * (', glue_collapse(ctrls, '+'), ")"))
  df_mod_centered = data.frame(scale(df, center = TRUE, scale = FALSE))
  m1 = felm(fmla, data = df_mod_centered) %>% robustify
  ols_est1 = summary(m1)$coefficients[2,1:2]
  return(ols_est1)
}

# %% ipw
ipw_reg = function(df, y, w, ctrls){
  fmla = as.formula(glue('{w} ~ ', glue_collapse(ctrls, '+')))
  m0 = glm(fmla, df, family = binomial())
  p = predict(m0, type = 'response')
  a = df[[w]]
  wt <- (a / p) + ((1 - a) / (1 - p))
  ipw_mod = felm(as.formula(glue("{y} ~ {w}")), data = df, weights = wt) %>% robustify
  ipw_est = summary(ipw_mod)$coefficients[2,1:2]
  return(ipw_est)
}

# %% ebal reg
eb_reg = function(df, y, w, ctrls){
  setorderv(df, w)
  X = df[, ..ctrls]; w = df[[w]]; y = df[[y]]
  out_eb <- ebalance(Treatment= w, X= X, print.level = -1)
  wt1 = c(out_eb$w, w[w==1])
  eff = felm(y ~ w, weights = wt1) %>% robustify()
  summary(eff)$coefficients[2,1:2]
}

# %%
############################################################
# ML-based estimation
############################################################

# main function to estimate effects using augmented data in  5 different ways
estimate_hdc = function(y, w, X, bak_X = NULL){
  # estimate 1, 2 - lassoEffect
  eff1 = rlassoEffect(X, y, w, method = "partialling out")
  eff2 = rlassoEffect(X, y, w, method = "double selection")
  # estimator 4 - knockoff selector
  y_ko = ko.glm(X, scale(y)); y_selected = ko.sel(y_ko)$estimation
  w_ko = ko.glm(X, w);        w_selected = ko.sel(w_ko)$estimation
  union_list = union(which(y_selected != 0), which(w_selected != 0))
  Xnew = X[, union_list]; eff4 = felm(y ~ cbind(w, Xnew)) %>% robustify
  # estimator 3- entropy balance
  # balance only on vars selected by the two Knockoff selectors
  try(out_eb <- ebalance(Treatment=w, X=Xnew, print.level = -1), silent=FALSE)
  if(exists("out_eb")){
    wt1 = c(out_eb$w, w[w==1])
    eff3 = felm(y ~ w, weights = wt1) %>% robustify
    eb_est = summary(eff3)$coefficients[2, 1:2] # Entropy
  } else {
    message("original EB did not converge, using backup")
    out_eb <- ebalance(Treatment=w, X= bak_X, print.level = -1)
    wt1 = c(out_eb$w, w[w==1])
    eff3 = felm(y ~ w, weights = wt1) %>% robustify
    eb_est = summary(eff3)$coefficients[2, 1:2] # Entropy
  }
    # %% collect estimate and SE
  res = rbind(
    summary(eff1)$coefficients[,  1:2], # PO
    summary(eff2)$coefficients[,  1:2], # DS
    summary(eff4)$coefficients[2, 1:2], # KOSEL
    eb_est
  )
  rownames(res) = c("Partial-Out", "Double Selection", "Knockoff Selector", "Knockoff EB")
  return(res)
}
