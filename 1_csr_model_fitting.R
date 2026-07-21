###################################################################
## Script for fitting multi-species models of plant CSR strategy ##
## Tom Mason, July 2026 ########################################### 
###################################################################

#############################
## 1 ## Set up environment ##
#############################

### Install and load all required packages
list.packages <- c("rstan","brms","cmdstanr","bayesplot","loo","ade4","usdm",
                   "terra","Ternary","ggplot2","dplyr","patchwork") 
# Function to install packages
check.install.packages <- function(list.packages){
  for(req.lib in list.packages){
    is.installed <- is.element(req.lib, installed.packages()[,1])
    if(is.installed == FALSE){install.packages(req.lib)}
    require(req.lib,character.only=TRUE)
  }
}
check.install.packages(list.packages)

### Set up directories - Change as appropriate
data_dir <- "C://Users/tomhe/Dropbox/Work/Jobs/WCMC/Work/Plant trait SDMs/Data"
output_dir <- "C://Users/tomhe/Dropbox/Work/Jobs/WCMC/Work/Plant trait SDMs/Output"
#data_dir <- "C://Users/tomm/Dropbox/Work/Jobs/WCMC/Work/Plant trait SDMs/Data"
#output_dir <- "C://Users/tomm/Dropbox/Work/Jobs/WCMC/Work/Plant trait SDMs/Output"


###############################
## 2 ## Load and manage data ##
###############################

### Load CSR data, for model fitting
CSRCompo_raw <- data.table::fread(paste0(data_dir,"/StanData.csv"))#,
sapply(CSRCompo_raw, function(x) sum(is.infinite(x)))
sapply(CSRCompo_raw, function(x) sum(is.na(x)))

# Remove selected columns
CSRCompo_raw <- CSRCompo_raw |>
  dplyr::select(
    -c("Life form","Growth form","Dispersal syndrome","Deciduousness")
  )

# Remove NAs and Inf values
CSRCompo_raw <- CSRCompo_raw |>
  dplyr::filter(
    dplyr::if_all(
      dplyr::everything(),
      ~ !is.na(.x) & !is.infinite(.x)
    )
  )
dim(CSRCompo_raw)

# Reformat data
CSRCompo <- CSRCompo_raw |>
  # rename component
  dplyr::rename(C=C,
                S=S,
                R=R) |>
  # Divide by 100 to rescale between 0 and 1, and add small offset
  dplyr::mutate(dplyr::across(C:R, ~ ./100 + 1e-10)) |>
  # Renormalize to ensure no exact zeros
  dplyr::mutate(dplyr::across(C:R, ~ ./rowSums(dplyr::across(C:R))))

# Z-transform predictors - only for numeric environmental predictors
plant_data_transformed <- CSRCompo |>
  dplyr::mutate(
    dplyr::across(
      -(1:9) & where(is.numeric),
      ~ as.numeric(scale(.x))
    )
  )

####################################################################
## 3 ## Variable selection using PCA and Correlation coefficients ##
####################################################################

### Run PCA
ndim = 4 # Retain first 4 PC
pca_data <- plant_data_transformed |>
  dplyr::select(-c(1:9))
pca_data <- plant_data_transformed |>
  dplyr::select(
    -c(1:9),
    -dplyr::where(is.character)
  )
pca <- ade4::dudi.pca(df= pca_data,center= F, scale = F, 
                     scannf= F,nf= ndim) 

### Inertia explained by each dimension
inertia = pca$eig/sum(pca$eig)*100
m <- data.frame(Comp=paste0('Dim',1:length(inertia)), inertia=inertia)%>%
  dplyr::slice(1:ndim)

### Evaluate absolute contribution of a variable to an axis
cont = ade4::inertia.dudi(pca,col.inertia=TRUE)$col.abs
ctr = cont[order(cont[,1],decreasing=FALSE),] # Sort in decreasing order
mm <- ctr |>
  tibble::rownames_to_column(var="Variable") |>
  tidyr::pivot_longer(cols=dplyr::contains("Axis"),
                      names_to="Axis",values_to="Contribution") |>
  dplyr::group_by(Axis) |>
  dplyr::arrange(dplyr::desc(Contribution),.by_group=TRUE) |>
  dplyr::ungroup() |>
  dplyr::mutate(Variable = factor(paste(Variable, Axis, sep = "__"),
                                  levels = rev(paste(Variable, Axis, sep = "__"))))

### Refine predictor set, based on their importance to PCs
preds <- c("pH_mean") ## Optionally, keep standard variables, e.g. most plant species have strict soil pH requirements
freq = 1
top = .75
preds <- ctr |>
  tibble::rownames_to_column(var="Variable") |>
  # for each dimension (axis) select variables in the top 75%
  dplyr::mutate(keep1 = Axis1 >= quantile(Axis1, probs=top), # Add predictors ranked at least once in the top quarter (75%) of contributing variables
                keep2= Axis2  >= max(Axis2)*.5, # Or where a few variables are dominant 
                keep3= Axis3  >= max(Axis3)*.5,
                keep4 = Axis4  >= max(Axis4)*.5) %>%
  dplyr::filter(rowSums(dplyr::across(keep1:keep4, ~ . %in% TRUE))>=freq) %>%
  dplyr::pull(Variable) %>%
  # add to standard variables
  append(preds,.)
preds

### Classify as variables of mean strategy and environmental filtering
# Mean strategy
patterns_mu<-c("pH","OCS","Nit","bio_1","bio_2","bio_3","bio_11","bio_12","bio_13","bio_14","bio_16",
                 "SBIO1","SBI10","SBIO11")
preds_mu <- preds[grep(paste(patterns_mu, collapse = "|"), preds)]
# Envtl filtering
patterns_phi<-c("slope","rough","TRI","def","bio_6","SBIO4","SBIO5","SBIO6","SBIO7")
preds_phi <- preds[grep(paste(patterns_phi, collapse = "|"), preds)]

### Now refine further, by assessing collinearity within each group of variables
# Mean strategy
cor_mu <- as.data.frame(plant_data_transformed)[,preds_mu]
cor_mu_mat <- cor(cor_mu,use = "pairwise.complete.obs",method = "pearson")
# Envtl filtering
cor_phi <- as.data.frame(plant_data_transformed)[,preds_phi]
cor_phi_mat <- cor(cor_phi,use = "pairwise.complete.obs",method = "pearson")

# Within each group of strongly correlated variables (Pearson's r > .7), choose most important according to PCA
### This step was performed manually ###
# Mean strategy
patterns_mu <-c("pH","OCS","Nit","bio_2","bio_14","SBIO11")
preds_mu <- preds[grep(paste(patterns_mu, collapse = "|"), preds)]
# Envtl filtering
patterns_phi <- c("def_mean","slope","SBIO5","SBIO6")
preds_phi <- preds[grep(paste(patterns_phi, collapse = "|"), preds)]


########################
## 4 ## Model fitting ##
########################

##################################################################
### OPTION a: Fit different combinations and select best model ###
##################################################################

### Create all combinations of variables

# Correlated variables
forbidden_groups <- list(
  # Envtl filtering variables
  c("slope","roughness","TRIrmsd"),
  c("def_max_mean","def_cumul_mean","def_mean","def_proba"),
  c("wc2.1_2.5m_bio_6","SBIO4_Temperature_Seasonality",
    "SBIO6_Min_Temperature_of_Coldest_Month",
    "SBIO7_Temperature_Annual_Range)"),
  # Main variables
  c("SBIO1_Annual_Mean_Temperature","SBIO10_Mean_Temperature_of_Warmest_Quarter",
    "SBIO11_Mean_Temperature_of_Coldest_Quarter","wc2.1_2.5m_bio_1",
    "wc2.1_2.5m_bio_3","wc2.1_2.5m_bio_11","wc2.1_2.5m_bio_12",
    "wc2.1_2.5m_bio_13","wc2.1_2.5m_bio_16")
)

# Function to identify valid combinations
is_valid_set <- function(set, forbidden_groups) {
  for (grp in forbidden_groups) {
    if (sum(set %in% grp) > 1) return(FALSE)
  }
  TRUE
}

# List of variable sets and their max lengths
var_sets <- list(
  mu    = list(vars = preds_mu,    max_k = 4),
  phi    = list(vars = preds_phi, max_k = 2)
)

# Generate all valid combinations for each set
valid_sets <- lapply(var_sets, function(vs) {
  all_sets <- unlist(
    lapply(1:vs$max_k, function(k) {
      combn(vs$vars, k, simplify = FALSE)
    }),
    recursive = FALSE
  )
  Filter(function(s) is_valid_set(s, forbidden_groups), all_sets)
})


### Convert into formulae
# For mean strategy
mu_formula <- lapply(valid_sets$mu, function(vars) {
  rhs <- paste0(vars, collapse = " + ")
  as.formula(paste("cbind(C,S,R) ~", rhs, "+ (1|AccSpeciesName)"))
})
# For envtl filtering
phi_formula <- lapply(valid_sets$phi, function(vars) {
  rhs <- paste0(vars,collapse = " + ")
  as.formula(paste("phi ~", rhs))
})

### Set up model
nburn=1000 # No. samples discarded from warm-up
nchains=3 # No. MCMC chains
niter=2500 # No. of iterations (length of chains)
ncores =nchains

### Fit models, running through each combination of variables
# This is done first for mu, with phi~1
# Then for phi, using best mu formulation

# Generate subsample to estimate loo from
set.seed(123)
sub_idx <- sample(dim(plant_data_transformed)[1], size = 1500, replace = FALSE)

###a Model fitting for mu
csr_mods_mu <- lapply(1:length(mu_formula),function(x) {

  #Model form
  brms_form = brms::bf(
    mu_formula[[x]],
    phi ~ 1
  )

  # Check structure
  str(plant_data_transformed[, c("C", "S", "R")])

  # Check each row sums to 1
  rowSums(plant_data_transformed[, c("C", "S", "R")]) |> range()

  #fit model
  fit <- brms::brm(brms_form,
                 data= plant_data_transformed,
                 family= brms::dirichlet(link= 'logit',link_phi='log',refcat= "R"), # Refcat is Ruderal (coefficients set to 0)
                 backend= "rstan",
                 warmup= nburn,
                 iter= niter,
                 chains= nchains,
                 cores= ncores,
                 save_pars= brms::save_pars(all= TRUE)
  )

  # Add LOO (cached inside model)
  ##fit<-add_criterion(fit,"loo",moment_match=TRUE,overwrite=TRUE)
  loo_ss <- loo_subsample(fit,cores = 4,seed = 123,subsample = sub_idx)
  fit$criteria$loo_subsample <- loo_ss

  # Save model WITH LOO attached
  saveRDS(fit,file=paste0(output_dir,"/CSR_model_mu_",x,".rds"))
  
  # Print progress
  message(sprintf("[%s] Model %d / %d",format(Sys.time(),
                                              "%Y-%m-%d %H:%M:%S"),x,length(x)))
  
  # Return model
  return(fit)
  
})

### For models with different formulation for mu, 
###  compare leave-one-out(loo) expected log predictive accuracy (elpd)
loo_list <- lapply(csr_mods_mu, function(x) x$criteria$loo_subsample)
which_best_mu <- which.max(sapply(loo_list, function(x) x$estimates["elpd_loo", "Estimate"]))
csr_best<-csr_mods_mu[which_best_mu][[1]]
saveRDS(csr_best,file=paste0(output_dir,"/CSR_model_mu_best.rds"))

###b Model fitting for phi
csr_mods_mu_phi <- lapply(1:length(phi_formula),function(x) {
  
  #Model form
  brms_form = brms::bf(
    mu_formula[[which_best_mu]],
    phi_formula[[x]]
  )
  
  # Check structure
  str(plant_data_transformed[, c("C", "S", "R")])
  
  # Check each row sums to 1
  rowSums(plant_data_transformed[, c("C", "S", "R")]) |> range()
  
  #fit model
  fit <- brms::brm(brms_form,
                 data= plant_data_transformed,
                 family= brms::dirichlet(link= 'logit',link_phi='log',refcat= "R"), # Refcat is Ruderal (coefficients set to 0)
                 backend= "rstan",
                 warmup= nburn,
                 iter= niter,
                 chains= nchains,
                 cores= ncores,
                 save_pars= brms::save_pars(all= TRUE)
  )
  
  # Add LOO (cached inside model)
  ##fit<-add_criterion(fit,"loo",moment_match=TRUE,overwrite=TRUE)
  loo_ss <- loo_subsample(fit,cores = 4,seed = 123,subsample = sub_idx)
  fit$criteria$loo_subsample <- loo_ss
  
  # Save model WITH LOO attached
  #saveRDS(fit,file=paste0(output_dir,"/CSR_model_mu_",x,".rds"))
  
  # Print progress
  message(sprintf("[%s] Model %d / %d",format(Sys.time(),
                                              "%Y-%m-%d %H:%M:%S"),x,length(x)))
  
  # Return model
  return(fit)
  
})

### For models with different formulation for mu, 
###  compare leave-one-out(loo) expected log predictive accuracy (elpd)
loo_list <- lapply(csr_mods_mu_phi, function(x) x$criteria$loo_subsample)
which_best_mu_phi <- which.max(sapply(loo_list, function(x) x$estimates["elpd_loo", "Estimate"]))
csr_best<-csr_mods_mu_phi[which_best_mu_phi][[1]]
saveRDS(csr_best,file=paste0(output_dir,"/CSR_model_mu_phi_best.rds"))

### Plot of fitted effects
csr_best <- readRDS(file=paste0(output_dir,"/CSR_model_mu_phi_best.rds"))
p1a <- plot(conditional_effects(csr_best,categorical= T,
                             effects= "wc2.1_2.5m_bio_2"),plot = F)[[1]]
p1b <- plot(conditional_effects(csr_best,categorical= T,
                               effects= "Nitrogem_mean"),plot = F)[[1]]
p1c <- plot(conditional_effects(csr_best,categorical= T,
                               effects= "OCS_mean"),plot = F)[[1]]
p1d <- plot(conditional_effects(csr_best,categorical= T,
                               effects= "pH_mean"),plot = F)[[1]]
p1 <- wrap_plots(p1a, p1b, p1c, p1d, ncol = 2)  # stacked vertically
p1


####################################
### OPTION b: Maximal model only ###
####################################

### Convert into formula
# For mean strategy
mu_rhs <- paste0("s(", preds_mu, ")", collapse = " + ")
mu_formula <- as.formula(paste("cbind(C,S,R) ~", mu_rhs, "+ (1|AccSpeciesName)"))
# For envtl filtering
phi_rhs <- paste0("s(", preds_phi, ")", collapse = " + ")
phi_formula <-  as.formula(paste("phi ~", phi_rhs))
# Best model
best_formula <- brms::bf(
  mu_formula,
  phi_formula
)

### Set up model
nburn = 1000 # No. samples discarded from warm-up
nchains = 3 # No. MCMC chains
niter = 4000 # No. of iterations (length of chains)
ncores = nchains

### Generate sub-sample to estimate loo from
# set.seed(123)
# sub_idx <- sample(dim(plant_data_transformed)[1], size = 1500, replace = FALSE)

### Checks 
str(plant_data_transformed[, c("C", "S", "R")]) # Correct structure
rowSums(plant_data_transformed[, c("C", "S", "R")]) |> range() # Check each row sums to 1

### Fit model
fit<-brms::brm(best_formula,
               data= plant_data_transformed,
               family= brms::dirichlet(link= 'logit',link_phi='log',refcat= "R"), # Refcat is Ruderal (coefficients set to 0)
               backend= "rstan",
               warmup= nburn,
               iter= niter,
               chains= nchains,
               cores= ncores,
               save_pars= brms::save_pars(all= TRUE)
)

# Add LOO (cached inside model)
fit<-add_criterion(fit,"loo",moment_match=TRUE,overwrite=TRUE)
#loo_ss <- loo_subsample(fit,cores = 4,seed = 123)#,subsample = sub_idx)
#fit$criteria$loo_subsample <- loo_ss
# Save model WITH LOO attached
saveRDS(fit,file=paste0(output_dir,"/CSR_model_maximal.rds"))

### Compare fit to null model
# Model formula
null_formula<-brms::bf(
  cbind(C, S, R) ~ 1 + (1 | AccSpeciesName),
  phi ~ 1 
)
# Fit null model
fit_null<-brms::brm(null_formula,
                    data= plant_data_transformed,
                    family= brms::dirichlet(link= 'logit',link_phi='log',refcat= "R"), # Refcat is Ruderal (coefficients set to 0)
                    backend= "rstan",
                    warmup= nburn,
                    iter= niter,
                    chains= nchains,
                    cores= ncores,
                    save_pars= brms::save_pars(all= TRUE)
)
# Add LOO (cached inside model)
fit_null<-add_criterion(fit_null,"loo",moment_match=TRUE,overwrite=TRUE)
# Save model WITH LOO attached
saveRDS(fit_null,file=paste0(output_dir,"/CSR_model_null.rds"))

### Compare leave-one-out(loo) expected log predictive accuracy (elpd)
csr_mods<-list(fit,fit_null)
loo_list<-lapply(csr_mods, function(x) x$criteria$loo)
csr_loos_comp <-loo::loo_compare(loo_list)



# ### PCA plots ###
# 
# #a Plot eigen values
# bplt <-ggplot2::ggplot(data=m,ggplot2::aes(x=Comp,y=inertia)) +
#   ggplot2::geom_bar(stat="identity",ggplot2::aes(fill=as.factor(inertia)),
#                     show.legend= FALSE) +
#   ggplot2::geom_text(ggplot2::aes(label=paste0(round(inertia,2),"%"),
#                                   y=inertia+0.1, fontface='bold'),
#                      vjust=0,
#                      color="black",
#                      position = ggplot2::position_dodge(1), size=4.5) +
#   ggplot2::theme_classic() +
#   ggplot2::labs(title="Barplot of Eigenvalues",
#                 subtitle=paste0("Using the first ",ndim," dimensions"),
#                 caption=paste0("The ",ndim," first coordinates (dimensions) explain ",
#                                round(sum(inertia[1:ndim]),2),"% of the variability")) +
#   ggplot2::theme(
#     text=ggplot2::element_text(family="serif"),
#     plot.title = ggplot2::element_text(face="bold", size=14, hjust = 0),
#     axis.text=ggplot2::element_text(colour="black", size=12),
#     axis.title=ggplot2::element_text(colour="black", size=16,face="bold")
#   ) +
#   ggplot2::xlab("Dimensions") + ggplot2::ylab("Percentage of explained variances") +
#   ggplot2::scale_fill_brewer(palette='Blues') +
#   ggplot2::scale_x_discrete(limits=paste0('Dim',1:ndim))
# bplt
# 
# ###b# Variable contributions
# 
# # Set up envt
# my_font_family = "serif"
# colr <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,'Blues'))
# cli::cli_progress_step("visualize contribution of each predictor")
# 
# # Plot
# bplt_axis <- ggplot2::ggplot(data=mm, ggplot2::aes(x=Variable,
#                                                    y=Contribution,
#                                                    group=Axis)) +
#   ggplot2::geom_bar(stat ="identity",
#                     ggplot2::aes(fill=as.factor(Contribution)),
#                     position = ggplot2::position_dodge(width=5),
#                     show.legend = FALSE) +
#   ggplot2::geom_text(ggplot2::aes(label=scales::percent(round(Contribution,3),
#                                                         scale =1,
#                                                         accuracy=0.01),
#                                   group=Axis),
#                      fontface='bold',
#                      hjust=-0.15,
#                      color="black",
#                      position = ggplot2::position_dodge(width=0.9),
#                      size=2.5) +
#   ggplot2::theme_classic() +
#   ggplot2::labs(title=paste0("Barplot of variables' contribution to each Axis"),
#                 subtitle=paste0("Using the first ",ndim," dimensions"),
#                 caption=paste0("The ",ndim," first coordinates (dimensions) explain ",
#                                round(sum(inertia[1:ndim]),2),"% of the variability"),
#                 x="Variable") +
#   ggplot2::theme(text=ggplot2::element_text(family=my_font_family, size=14),
#                  plot.title = ggplot2::element_text(hjust = 0, size=24),
#                  strip.background =ggplot2::element_rect(fill="gray"),
#                  strip.text = ggplot2::element_text(colour = 'black', size=14),
#                  axis.ticks = ggplot2::element_line(colour = "gray"),
#                  panel.grid = ggplot2::element_line(colour="gray", linewidth=0.1),
#                  panel.background = ggplot2::element_rect(colour="gray", fill="transparent"),
#                  #plot.background = ggplot2::element_rect(fill="transparent",color=NA),
#                  legend.position="bottom"
#   ) +
#   ggplot2::coord_flip() +
#   ggplot2::xlab("Variables") + ggplot2::ylab("% of contribution") +
#   ggplot2::scale_fill_manual(values=colr(nlevels(as.factor(mm$Contribution)))) +
#   ggplot2::scale_x_discrete(labels= ~ gsub("__.+$", "", .x)) +
#   ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult=c(0.001,0.1))) +
#   ggplot2::facet_wrap(~Axis, ncol=2, scales="free")
# bplt_axis


