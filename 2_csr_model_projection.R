###############################################################################
## Script for projecting plant species CSR strategy over environmental space ##
## and comparing with independent species' occurrences ########################
## Tom Mason, July 2026 ####################################################### 
###############################################################################

#############################
## 1 ## Set up environment ##
#############################

### Install and load CRAN packages
list.packages <- c("devtools","rstan","brms","cmdstanr","bayesplot","loo","ade4","usdm",
                   "terra","Ternary","ggplot2","dplyr","patchwork","sf","sfheaders",
                   "effsize","ggtern") 
# Function to install packages
check.install.packages <- function(list.packages){
  for(req.lib in list.packages){
    is.installed <- is.element(req.lib, installed.packages()[,1])
    if(is.installed == FALSE){install.packages(req.lib)}
    require(req.lib,character.only=TRUE)
  }
}
check.install.packages(list.packages)

### Install and load github packages
#pak::pak("matildabrown/rWCVP")
#pak::pak("matildabrown/rWCVPdata")
#pak::pak("IanOndo/TDWG")
#pak::pak("IanOndo/RGeodata")
list.packages <- c("rWCVP","rWCVPdata","TDWG","RGeodata")
check.install.packages(list.packages)

### Set up directories - Change as appropriate
data_dir <- "O://f01_projects_active/Global/p09623_CreatingFoundationSystems4EnvironmentalAI/work_in_progress/Progress July 2026/Data"
output_dir <- "O://f01_projects_active/Global/p09623_CreatingFoundationSystems4EnvironmentalAI/work_in_progress/Progress July 2026/Output"

### Load additional functions

#a Utility functions
meanAndsd <-function(x){
  cbind(mean=colMeans(x, na.rm=T), sd=apply(x,2,sd, na.rm=T))
}
obtain <-`[`

#b Auxiliary function to reconciliate taxonomy when creating species' range map
reconcile_names_withWCVP <- function(taxon, taxon_rank = c("species", "genus", "family"), wcvp_names=NULL){
  
  if(is.null(wcvp_names)) wcvp_names <- rWCVPdata::wcvp_names
  
  taxon_rank = match.arg(taxon_rank)
  
  taxon_name <- gsub("_"," ",taxon)
  
  taxon_name_accepted_id <- wcvp_names |>
    dplyr::filter(taxon_name==!!taxon_name) |>
    dplyr::pull(accepted_plant_name_id)
  
  taxon_name <- wcvp_names |>
    dplyr::filter(accepted_plant_name_id %in% taxon_name_accepted_id,
                  taxon_status=="Accepted") |>
    dplyr::pull(taxon_name)
  
  taxon <- switch(taxon_rank,
                  species= taxon_name[1],
                  genus  = wcvp_names |> dplyr::filter(taxon_name==!!taxon_name[1]) |> dplyr::pull(genus),
                  family = wcvp_names |> dplyr::filter(taxon_name==!!taxon_name[1]) |> dplyr::pull(family))
  return(taxon)
}

#c Main function to generate species' range map
make_projection_domain_fromWCVP <- function(taxon,
                                            taxon_rank = c("species", "genus", "family", "order", "higher"),
                                            code_level = c("LEVEL3_COD","LEVEL2_COD","LEVEL1_COD"),
                                            output_dir=NULL,
                                            dissolve=FALSE,
                                            native = TRUE,
                                            introduced = TRUE,
                                            location_doubtful= FALSE,
                                            reconcile = TRUE,
                                            #use_ecoregions=TRUE,
                                            verbose=TRUE){
  
  if(missing(taxon))
    stop("Please provide the name of the species you would like the distribution area from.")
  
  code_level = match.arg(code_level)
  
  # update WCVP data if out of date
  if(!rWCVPdata::wcvp_check_version()){
    message("updating package `rWCVPdata`")
    unloadNamespace("rWCVPdata") #to be sure
    if(!"devtools" %in% names(installed.packages()[,3])) install.packages("devtools")
    devtools::install_github("matildabrown/rWCVPdata")
  }
  
  taxon_name <- gsub("_"," ",taxon)
  
  if(reconcile){
    taxon_name <- reconcile_names_withWCVP(taxon)
  }
  
  wcvp_dist <- try(rWCVP::wcvp_distribution(taxon = taxon_name,
                                            taxon_rank = taxon_rank,
                                            native = native,
                                            introduced = introduced,
                                            location_doubtful = location_doubtful))
  
  if(inherits(wcvp_dist,"try-error")){
    warnings("Cannot retrieve species distribution from WCVP")
    return(NULL)
  }
  wcvp_code <- unique(wcvp_dist$LEVEL3_COD)
  y = switch(code_level,
             "LEVEL1_COD"= TDWG:::tdwg_level1 |> dplyr::filter(LEVEL1_COD %in% wcvp_code),
             "LEVEL2_COD"= TDWG:::tdwg_level2 |> dplyr::filter(LEVEL2_COD %in% wcvp_code),
             "LEVEL3_COD"= TDWG:::tdwg_level3 |> dplyr::filter(LEVEL3_COD %in% wcvp_code))
  
  projection_domain <- sf::st_filter(x = RGeodata:::ecoregions_split,  ## Added an :
                                     y = y,
                                     .predicate=sf::st_overlaps)
  if(dissolve){
    tryCatch({
      projection_domain %<>%
        sf::st_set_crs(4326) %>%
        dplyr::summarise()
    },error=function(err){
      sf::sf_use_s2(FALSE)
      projection_domain %<>%
        sf::st_set_crs(4326) %>%
        dplyr::summarise()
    },finally = {sf::sf_use_s2(TRUE)})
  }
  
  projection_domain %<>%
    sf::st_make_valid() %>%
    sfheaders::sf_remove_holes()  ## Added sfheaders
  
  # ensure that the output directory provided exists
  if(!is.null(output_dir)){
    if(!dir.exists(output_dir)){
      warning(paste("Output directory:", output_dir,"does not exist or could not be found."));flush.console()
    }else
      saveRDS(projection_domain, file.path(output_dir,paste0(species_name,".rds")))
  }
  
  return(projection_domain)
  
}


#####################################
## 2 ## Load fitted model and data ##
#####################################

### Load best CSR model - for spatial projection
csr_best <- readRDS(file=paste0(output_dir,"/CSR_model_mu_phi_best.rds"))
preds <- c(all.vars(csr_best$formula$form)[-(1:3)],
         c(all.vars(csr_best$formula$pforms$phi))[-1])

### Load model data 
CSRCompo_raw <- CSRCompo_raw_orig <- data.table::fread(paste0(data_dir,"/StanData.csv"))#,
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
  # Rename component
  dplyr::rename(C=C,
                S=S,
                R=R) |>
  # Divide by 100 to rescale between 0 and 1, and add small offset
  dplyr::mutate(dplyr::across(C:R, ~ ./100 + 1e-10)) |>
  # Renormalize to ensure no exact zeros
  dplyr::mutate(dplyr::across(C:R, ~ ./rowSums(dplyr::across(C:R))))

### Load environmental spatial data - for spatial projection
env <- rast(paste0(data_dir,"/explvars.tif"))

# Fix names in spatial data
replacement_vars <- c("nitrogen_0-5cm_mean" = "Nitrogem_mean",
                     "ocs_0-30cm_mean" = "OCS_mean",
                     "phh2o_0-5cm_mean" = "pH_mean")
names(env)<-stringr::str_replace_all(names(env),replacement_vars)

# Subset spatial data
env_mod <- terra::subset(env,subset=intersect(preds,names(env)))
#plot(env_mod)

# z-transform environmental parameters (as in model)
env_rescale <- CSRCompo |>
  dplyr::select(dplyr::where(is.numeric))|>
  dplyr::select(dplyr::any_of(preds)) |>
  meanAndsd()|>
  obtain(names(env_mod),)
env_z <- terra::scale(env_mod,center= unname(env_rescale[,"mean"]),
                    scale= unname(env_rescale[,"sd"]))

### Load independent occurrence records- for validation
utilised_spp <- list.files(paste0(data_dir,"/Utilised plants/Occurrences/")) |> tools::file_path_sans_ext() |> gsub("_", " ", x = _)
spp.names <- unique(CSRCompo$AccSpeciesName)
utilised_spp_sub <- spp.names[which(spp.names %in% utilised_spp)]


#############################
## 3 ## Spatial projection ##
#############################

### Define prediction function 
predfun<-function(cat) {
  function(model, data, ...) {
    mu  <- brms::posterior_epred(model, newdata = data, cores = 10, ...)
    phi <- brms::posterior_epred(model, newdata = data, cores = 10, dpar = "phi", ...)
    
    # Full summary matrix: nobs x (n_cats + phi)
    out <- apply(mu, c(2, 3), mean) |>
      cbind(matrix(colMeans(phi), dimnames = list(NULL, "phi")))
    
    # Select single column for terra
    out[, cat]
  }
}

### Select species for projection
species_mod <- utilised_spp_sub

# Correcting a few species names to match names in species' data
species_range <- gsub("(\\w+) caffra", "\\1 afra", species_mod)
species_range[species_range == "Pentanema britannicum"] <- "Pentanema britannica"

### Create dataframe to save output
results <- vector("list", length(species_mod))

### Run through selected species
for (x in 1:length(species_mod)){
 
  #a# Load species occurrence and range data 
  
  ## Select one species
  one_sp_mod <- species_mod[x]
  one_sp_range <- species_range[x]
  print(paste0("Starting ",x,". ",one_sp_mod))
  
  ## Load species point records
  # Calibrated
  sp_xy_calib <- CSRCompo[CSRCompo$AccSpeciesName==one_sp_mod,c("Longitude","Latitude")]
  sp_xy_calib <- terra::vect(sp_xy_calib, geom = c("Longitude", "Latitude"), 
                       crs = "EPSG:4326")
  # Independent
  sp_xy <- read.csv(paste0(data_dir,"/Utilised plants/Occurrences/",gsub(" ","_",x=one_sp_mod),".csv"))
  sp_xy <- sp_xy[!duplicated(sp_xy[, c("decimalLatitude", "decimalLongitude")]), ]
  sp_xy <- terra::vect(sp_xy, geom = c("decimalLongitude","decimalLatitude"), 
                       crs = "EPSG:4326")
  
  ## Only proceed if there are sufficient occurrence records for validation
  sp_xy_n <- dim(sp_xy)[1]
  min_n <- 10
  if (sp_xy_n < min_n) {
    print(paste0("Insufficient species records- skipping ",one_sp_range))
    results[[x]] <- NULL
    next
  }
  
  ## Compute species range
  one_range <- make_projection_domain_fromWCVP(taxon=one_sp_range,
                                      taxon_rank = "species",
                                      code_level = "LEVEL3_COD",
                                      dissolve=T,
                                      native = T,
                                      introduced = T)
  
  ## Don't proceed if range is too large (due to memory limitations)
  max_area <- 5e13
  range_area <- as.numeric(st_area(one_range)) # in m² if CRS is projected
  if (range_area > max_area) {
    print(paste0("Range too large- skipping ",one_sp_range))
    results[[x]] <- NULL
    next
  }
  
  ## Only proceed if occurrences and range sufficiently overlap 
  sp_xy <- project(sp_xy, crs(vect(one_range)))
  overlap <- relate(sp_xy, vect(one_range), relation = "intersects")
  pct_overlap <- mean(overlap) * 100
  if (pct_overlap < 80) {
    print(sprintf("Only %.1f%% of points overlap with range - skipping %s", pct_overlap, one_sp_range))
    results[[x]] <- NULL
    next 
  }
  
  ## Crop spatial layers to extent of species' ranges
  env_crop <- terra::crop(env_z,one_range,mask=F)
  env_crop <- aggregate(env_crop,fact=10,fun="median",na.rm=T) # Aggregate for quicker result
  env_noInf <- terra::subst(env_crop, c(-Inf, Inf),NA) # Remove inf values

  
  #b# Run projection
  
  ## Get all category names + phi dynamically
  cats <- c(dimnames(brms::posterior_epred(csr_best, 
                                           newdata = as.data.frame(env_noInf, na.rm=TRUE)[1,],
                                           allow_new_levels = TRUE))[[3]], "phi")
  ## Project
  csr_pred <- lapply(cats, function(cat) {
    terra::predict(env_noInf, csr_best, predfun(cat),
                   allow_new_levels = TRUE, ndraws = 100, na.rm = TRUE,
                   const = data.frame(AccSpeciesName = one_sp_mod),
                   cpkgs = "brms")
  }) |> rast()
  names(csr_pred) <- cats
  
  ## Define Dirichlet density probability function in Stan, for faster computation
  fun_code <-'
  functions {
  /* dirichlet-logit log-PDF
  * Args:
  * y: vector of real response values
  * mu: vector of category logit probabilities
  * phi: precision parameter
  * Returns:
  * a scalar to be added to the log posterior
  */
  real dirichlet_logit_lpdf(vector y, vector mu, real phi) {
  for (i in 1:num_elements(mu)) { if(is_nan(mu[i])) return not_a_number(); };
  return dirichlet_lpdf(y | mu * phi);
  }
  }
  '

  # Expose the function to R
  rstan::expose_stan_functions(rstan::stanc(model_code = fun_code))
  
  ## Define model projections
  mu <- csr_pred[[1:3]] # Mean CSR strategy
  phi <- csr_pred[[4]] # Precision parameter
  
  ## Calculate mean CSR strategy of target species, from CSR data (rescaled between 0 and 1)
  target_csr <- as.numeric(apply(CSRCompo[CSRCompo$AccSpeciesName==one_sp_mod,c("C","S","R")],2,mean))
  
  ## Create a vectorised version of our density probability function
  ddirch <- Vectorize(
    function(C,S,R,phi) return(dirichlet_logit_lpdf(target_csr, c(C,S,R), phi))
  )
  
  ## Compute habitat suitability, here using percentile-based rescaling to deal with outliers
  suitab <- terra::lapp(c(mu, phi),fun=ddirch,usenames=TRUE)
  q_low  <- global(suitab, quantile, probs = 0.01, na.rm = TRUE)[[1]]
  q_high <- global(suitab, quantile, probs = 0.99, na.rm = TRUE)[[1]]
  suitab_clamp <- clamp(suitab, lower = q_low, upper = q_high)
  suitab_centered <- suitab_clamp - global(suitab_clamp, max, na.rm = TRUE)[[1]]
  proba_sc <- exp(suitab_centered)
  proba_sc_mask <- mask(proba_sc,one_range) # Mask
  
  #c# Output results
  
  ## Habitat suitability raster
  writeRaster(proba_sc_mask, paste0(output_dir,"/Habitat suitability/",one_sp_mod,".tif"),
              overwrite=T)
  
  ## Plot 
  pdf(paste0(output_dir,"/Plots/",one_sp_mod,".pdf"),
       width=10,height=8)
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE))
  par(oma = c(0, 0, 3, 0), bty = "n")  # outer margin at top for overall title
  
  #i Habitat suitability
  plot(proba_sc_mask,axes=F)
  mtext("a) Trait-informed habitat suitability", side = 3, line = 1, cex = 1, font = 2,padj=-2)
  
  #ii Species occurrences
  plot(vect(one_range),xaxs = "i", yaxs = "i",axes=F)
  points(sp_xy, pch = 21, col = "firebrick", cex = 1.2, lwd = 2)
  points(sp_xy_calib, pch = 21, col = "steelblue", cex = 1.2, lwd = 2)
  mtext("b) Independent species occurrences", side = 3, line = 1, cex = 1, font = 2,padj=-2)
  legend("bottomright",c("Calibration","Validation"), pch=21,
         col=c("steelblue","firebrick"),bty="n",xpd=NA)
  
  #iii Histograms
  
  # xy
  proba_xy_dat <- extract(proba_sc_mask, sp_xy, na.rm = TRUE)
  proba_xy <- proba_xy_dat[!is.na(proba_xy_dat$lyr1), "lyr1"]
  hist_xy <- density(proba_xy, bw = .1, from = 0, to = 1)
  hist_xy$y <- (hist_xy$y - min(hist_xy$y)) / (max(hist_xy$y) - min(hist_xy$y))
  # range
  proba_range_na <- as.numeric(values(proba_sc_mask))
  proba_range <- proba_range_na[!is.na(proba_range_na)]
  hist_range <- density(proba_range, bw = .1, from = 0, to = 1)
  hist_range$y <- (hist_range$y - min(hist_range$y)) / (max(hist_range$y) - min(hist_range$y))
  # now plot
  plot(NULL, xlim = c(0, 1), ylim = c(0, 1.01), xaxs = "i", yaxs = "i",bty="l",
       ylab = "Scaled probability density", xlab = "Trait-informed habitat suitability")
  polygon(x = c(hist_range$x, rev(hist_range$x)),
          y = c(hist_range$y, rep(0, length(hist_range$y))),
          col = adjustcolor("gray50", alpha.f = .5), border = NA)
  lines(x = hist_xy$x, y = hist_xy$y, col = "firebrick", lwd = 2)
  legend("topright", c("Occurrences", "Species range"),
         lty = c(1, NA),pch = c(NA, 15), lwd=c(2,NA),bty="n",
         col = c("firebrick", adjustcolor("gray50", alpha.f = 0.5)))
  mtext("c) Habitat suitability of occurrences", side = 3, line = 1, cex = 1, font = 2)
  
  # Calculate Cohen's d metric
  #ks.test(proba_xy, proba_range)
  proba_range_sub <- sample(proba_range, min(length(proba_range), length(proba_xy))) # subsample 
  cohen_d<-cohen.d(proba_xy, proba_range_sub)
  mtext(paste0("Cohen's d = ",round(cohen_d$estimate,3),side=3))
  
  # Alternative metric- difference between peaks
  peak_xy <- hist_xy$x[which.max(hist_xy$y)]
  peak_range <- hist_range$x[which.max(hist_range$y)]
  peak_diff <- peak_xy - peak_range
  
  # Overall title
  mtext(bquote(italic(.(one_sp_range))), side = 3, line = 1, cex = 1.4, font = 2, outer = TRUE)
  
  dev.off()
  
  ## Validation results
  # Add row to dataframe 
  results[[x]] <- data.frame(
    species = one_sp_mod,
    life_form = CSRCompo_raw_orig[CSRCompo_raw_orig$AccSpeciesName==one_sp_mod,'Life form'][1],
    growth_form = CSRCompo_raw_orig[CSRCompo_raw_orig$AccSpeciesName==one_sp_mod,'Growth form'][1],
    dispersal_synd = CSRCompo_raw_orig[CSRCompo_raw_orig$AccSpeciesName==one_sp_mod,'Dispersal syndrome'][1],  
    decid = CSRCompo_raw_orig[CSRCompo_raw_orig$AccSpeciesName==one_sp_mod,'Deciduousness'][1],
    model_n = dim(CSRCompo[CSRCompo$AccSpeciesName==one_sp_mod,])[1],
    c_mean = target_csr[1],
    s_mean = target_csr[2],
    r_mean = target_csr[3],
    cohens_d = cohen_d$estimate,
    cohens_d_min = cohen_d$conf.int[1],
    cohens_d_max = cohen_d$conf.int[2],
    peak_diff,
    n_points = length(proba_xy),
    n_range = length(proba_range),
    row.names=NULL
  )
  # Save
  write.csv(results[[x]],paste0(output_dir,"/Validation/",one_sp_mod,".csv"))
  
  
  ## End
  print(paste0("Completed ",x,". ",one_sp_mod))
  
}

#dev.off()

### Combine validation data
validation <- list.files(paste0(output_dir, "/Validation/"), full.names = TRUE)
validation_dat <- do.call(rbind, lapply(validation, read.csv))


### Plotting ###

# ### Density plot of predictive performance
# pdf(paste0(output_dir,"/Density plot Cohen's d.pdf"),
#     width=8,height=6)
# 
# dens <- with(validation_dat, density(cohens_d))
# n <- 1000
# breaks <- seq(min(dens$x), max(dens$x), length.out = n + 1)
# # Find which slice index corresponds to x=0
# zero_idx <- which.min(abs(breaks - 0))
# # Build palette with white exactly at zero_idx
# cols <- c(
#   colorRampPalette(c("firebrick", "white"))(zero_idx),
#   colorRampPalette(c("white", "steelblue"))(n - zero_idx + 1)[-1]
# )
# plot(dens, type = "n", main = "", xlab = "Cohen's d")
# for (i in 1:n) {
#   xi <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
#   yi <- c(0, 0,
#           approx(dens$x, dens$y, breaks[i+1])$y,
#           approx(dens$x, dens$y, breaks[i])$y)
#   polygon(xi, yi, col = cols[i], border = NA)
# }
# lines(dens, col = "gray50")
# segments(0, 0, 0, approx(dens$x, dens$y, 0)$y, lty = 2)
# 
# dev.off()


# # Ternary plot of predictive performance
# pdf(paste0(output_dir,"/CSR Cohen's d.pdf"),
#     width=8,height=6)
# 
# eps <- 0.018
# validation_dat$c_mean2 <- validation_dat$c_mean * (1 - 3*eps) + eps
# validation_dat$s_mean2 <- validation_dat$s_mean * (1 - 3*eps) + eps
# validation_dat$r_mean2 <- validation_dat$r_mean * (1 - 3*eps) + eps
# validation_dat$cohens_d_clamped <- pmax(pmin(validation_dat$cohens_d, 1), -1)
# ggtern(data = validation_dat,
#        aes(x = c_mean2,y = s_mean2,z = r_mean2,fill = cohens_d_clamped)) +
#   geom_point(shape = 21,size = 3,colour = "gray50",stroke = 0.5) +
#   scale_fill_gradient2(name = "Cohen's d",low = "firebrick",mid = "white", 
#                        high = "steelblue",midpoint = 0,limits = c(-1, 1)) +
#   labs(T = "C", L = "S", R = "R") +
#   theme_bw() +
#   theme(tern.axis.arrow.show = TRUE,tern.panel.background = element_blank(),
#         tern.axis.line = element_line(linewidth = 0.3))
# 
# dev.off()

# ## Reorganise plots according to validation performance
# # Create subfolders
# plot_dir <- paste0(output_dir, "/Validation/")
# dirs <- c("very_poor", "poor", "good", "very_good")
# lapply(paste0(plot_dir, dirs), dir.create, showWarnings = FALSE)
#
# # Classify and move each species
# for (i in 1:nrow(validation_dat)) {
#   species <- validation_dat$species[i]  # adjust column name as needed
#   d <- validation_dat$cohens_d[i]
#   
#   subfolder <- ifelse(d < -0.5, "very_poor",
#                       ifelse(d < 0,    "poor",
#                              ifelse(d < 0.5,  "good",
#                                     "very_good")))
#   
#   file.copy(
#     from = paste0(output_dir,"/Plots/", species, ".pdf"),
#     to   = paste0(output_dir,"/Plots/", subfolder, "/", species, ".pdf")
#   )
# }



# ## Global projection
# 
# env_glob<-aggregate(env_z,fact=20,fun="median",na.rm=T)
# env_glob <-terra::subst(env_glob, c(-Inf, Inf),NA) # Remove inf values
# plot(env_glob[[1]])
# 
# # Get all category names + phi dynamically
# cats <- c(dimnames(brms::posterior_epred(csr_best, 
#                                          newdata = as.data.frame(env_glob, na.rm=TRUE)[1,],
#                                          allow_new_levels = TRUE))[[3]], "phi")
# # Project
# csr_pred <- lapply(cats, function(cat) {
#   terra::predict(env_glob, csr_best, predfun(cat),
#                  allow_new_levels = TRUE, ndraws = 100, na.rm = TRUE,
#                  #const = data.frame(AccSpeciesName = one_sp_mod),
#                  cpkgs = "brms")
# }) |> rast()
# names(csr_pred) <- cats
# #plot(csr_pred,axes=F)
# 
# # Defining Dirichlet density probability function in Stan
# # (for fast computation)
# fun_code <-'
#   functions {
#   /* dirichlet-logit log-PDF
#   * Args:
#   * y: vector of real response values
#   * mu: vector of category logit probabilities
#   * phi: precision parameter
#   * Returns:
#   * a scalar to be added to the log posterior
#   */
#   real dirichlet_logit_lpdf(vector y, vector mu, real phi) {
#   for (i in 1:num_elements(mu)) { if(is_nan(mu[i])) return not_a_number(); };
#   return dirichlet_lpdf(y | mu * phi);
#   }
#   }
#   '
# 
# # Expose the function to R
# rstan::expose_stan_functions(rstan::stanc(model_code = fun_code))
# 
# # Define model projections
# mu <- csr_pred[[1:3]] # Mean CSR strategy
# phi <- csr_pred[[4]] # Precision parameter
# 
# #iii# Output results
# 
# #b# Plot 
# pdf(paste0(output_dir,"/Global CSR.pdf"),
#     width=15,height=8)
# 
# mu_clamped <- clamp(mu, lower = 0, upper = .8)
# plot(mu_clamped,range=c(0,.8))
# 
# dev.off()






