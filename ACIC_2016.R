rm(list = ls())
library(LalRUtils)
libreq(data.table, aciccomp2016, assertthat, fastDummies, glmnet, hdm, ebal,
    kosel, tictoc)
set.seed(42)
notif = \(x) RPushbullet::pbPost("note", x)

# %% parameter grid
npar = 77  # 1:77
nsim = 100 # 1:100
ijs = expand.grid(1:npar, 1:nsim)

# %% covariate matrix (basic and interacted)
k = 2
X = dummy_cols(input_2016, remove_first_dummy = T,
        remove_selected_columns = T) %>% as.matrix()
# matrix with interactions
XX = model.matrix(as.formula(glue::glue("~.^{k} - 1")), data.frame(scale(X, scale = F)))
keepX = apply(XX, 2, var) %>% .[.>1] %>% names()
XX = XX[, keepX]
hc = findCorrelation(cor(XX), cutoff=0.9) %>% sort
XX = XX[,-hc]

# %% knockoff selector fn
koSelect = function(X, y, w){
  y_ko       = ko.glm(X, scale(y));
  y_selected = ko.sel(y_ko)$estimation
  w_ko       = ko.glm(X, w);
  w_selected = ko.sel(w_ko)$estimation
  union(which(y_selected != 0), which(w_selected != 0))
}

# %% entropy balancing function
eb_run = \(ij, X, kosel = F){
  # simulation
  sim = dgp_2016(X, ij[1], ij[2])
  # outcome and treatment
  y = sim$y; z = sim$z
  # true SATT
  true_satt = (mean(sim$y.1[z == 1]) - mean(sim$y.0[z == 1]))
  if(kosel){ Xmat = X[, koSelect(X, y, z)]
  } else{ Xmat = X}
  # run ebal - fails often
  status = try({
      time = system.time({
        ebal.out = ebalance(Treatment=z, X=Xmat)
      })
  })
  # failure
  if (is.error(status)) { satt_hat = NA; bias = NA
  } else {
      # get weights
      weights0 = ebal.out$w / sum(ebal.out$w)
      weights1 = rep(1 / sum(z == 1), sum(z == 1))
      # estimate satt and standardized satt
      satt_hat = (weighted.mean(y[z == 1], weights1) -
                  weighted.mean(y[z == 0], weights0))
      bias = satt_hat - true_satt
  }
  c(satt_hat, bias)
}
# %% dry run
ijs %>% .[sample(nrow(.), 1), ] %>% as.numeric %>% eb_run(X)
# %%

# %% parallel loop for vanilla eb
libreq(foreach, doParallel, tictoc)
cl = makeCluster(detectCores() - 1) #not to overload your computer
registerDoParallel(cl)
nr = nrow(ijs)

tic()
ebalRes = foreach(k=1:nr, .combine=rbind,
    .packages = c("ebal", "kosel", "aciccomp2016", "assertthat")
  ) %dopar% {
  ij = ijs[k, ] %>% as.numeric()
  eb_run(ij, X)
}
toc()
#stop cluster
stopCluster(cl)

notif("basic sims done")
save(ebalRes, file = "eb1.rds")

# %% EB with big matrix
cl = makeCluster(detectCores() - 1) #not to overload your computer
registerDoParallel(cl)
nr = nrow(ijs)

tic()
ebalRes = foreach(k=1:nr, .combine=rbind,
    .packages = c("ebal", "kosel", "aciccomp2016", "assertthat")
  ) %dopar% {
  ij = ijs[k, ] %>% as.numeric()
  eb_run(ij, XX)
}
toc()
#stop cluster
stopCluster(cl)

notif("eb bigmat sims done")
save(ebalRes, file = "eb2.rds")

# %% EB + kosel
cl = makeCluster(detectCores() - 1) #not to overload your computer
registerDoParallel(cl)
nr = nrow(ijs)

tic()
ebalRes = foreach(k=1:nr, .combine=rbind,
    .packages = c("ebal", "kosel", "aciccomp2016", "assertthat")
  ) %dopar% {
  ij = ijs[k, ] %>% as.numeric()
  eb_run(ij, XX, kosel = T)
}
toc()
#stop cluster
stopCluster(cl)

notif("eb 3 sims done")
save(ebalRes, file = "eb3.rds")

# %%
