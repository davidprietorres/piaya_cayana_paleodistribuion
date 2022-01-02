library(dplyr) 
library(kuenm)
library(raster)

##########################Separación de los datos de presencia################
rm(list=ls()) #elimina TODOS los objetos del ambiente de trabajo

setwd("C:/project_sigs/Piaya/modelos/piaya_mexicana/")
dir() # check what is in your working directory

############################################
######### The model creation ###############
############################################
# Preparing variables to be used in arguments
file_name <- "piaya_mexicana"
kuenm_start(file.name = file_name)

# Variables with information to be used as arguments. Change "YOUR/DIRECTORY" by your actual directory.
occ_joint <- "piaya_mexicana_joint.csv"
occ_tra <- "piaya_mexicana_train.csv"
M_var_dir <- "M"
batch_cal <- "Candidate_models"
out_dir <- "Candidate_Models"
reg_mult <- c(seq(0.4, 2, 0.2), seq(2, 6, 0.5),8)
f_clas <- "all"
args <- NULL # e.g., "maximumbackground=20000" for increasing the number of pixels in the bacground or


# note that some arguments are fixed in the function and should not be changed
maxent_path <- "C:/Users/dprie/Documents/R/win-library/4.0/dismo/java"
wait <- FALSE
run <- TRUE
kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
          out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, args = args,
          maxent.path = maxent_path, wait = wait, run = run)
warnings()

############################################
######### The model evaluation #############
############################################
occ_test <- "piaya_mexicana_test.csv"
out_eval <- "Calibration_results"
threshold <- 5
rand_percent <- 50
iterations <- 500
kept <- TRUE
selection <- "OR_AICc"
paral_proc <- FALSE # make this true to perform pROC calculations in parallel, recommended
# only if a powerfull computer is used (see function's help)
# Note, some of the variables used here as arguments were already created for previous function
cal_eval <- kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, batch = batch_cal,
                        out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, iterations = iterations,
                        kept = kept, selection = selection, parallel.proc = paral_proc)


############################################
######### The model creation ###############
############################################
batch_fin <- "Final_models"
mod_dir <- "Final_Models"
rep_n <- 10
rep_type <- "Bootstrap"
jackknife <- FALSE
out_format <- "logistic"
project <- TRUE
G_var_dir <- "G"
ext_type <- "all"
write_mess <- TRUE
write_clamp <- TRUE
wait1 <- FALSE
run1 <- TRUE
args <- NULL # e.g., "maximumbackground=20000" for increasing the number of pixels in the bacground or
# "outputgrids=false" which avoids writing grids of replicated models and only writes the 
# summary of them (e.g., average, median, etc.) when rep.n > 1
# note that some arguments are fixed in the function and should not be changed
# Again, some of the variables used here as arguments were already created for previous functions
kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, batch = batch_fin, rep.n = rep_n,
          rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir, out.format = out_format, project = project,
          G.var.dir = G_var_dir, ext.type = ext_type, write.mess = write_mess, write.clamp = write_clamp, 
          maxent.path = maxent_path, args = args, wait = wait1, run = run1)


############################################
####### The model evaluation final #########
############################################
occ_ind <- "piaya_mexicana_ind.csv"
replicates <- TRUE
out_feval <- "Final_Models_evaluation"
# Most of the variables used here as arguments were already created for previous functions
fin_eval <- kuenm_feval(path = mod_dir, occ.joint = occ_joint, occ.ind = occ_ind, replicates = replicates,
                        out.eval = out_feval, threshold = threshold, rand.percent = rand_percent,
                        iterations = iterations, parallel.proc = paral_proc)

best <- read.csv("Calibration_results/best_candidate_models_OR_AICc.csv")
knitr::kable(best, caption = "Models selected based on significance, omission rates, and AICc, in that order.")



##########################Separación de los datos de presencia################
rm(list=ls()) #elimina TODOS los objetos del ambiente de trabajo

setwd("C:/project_sigs/Piaya/modelos/piaya_therma/")
dir() # check what is in your working directory

############################################
######### The model creation ###############
############################################
# Preparing variables to be used in arguments
file_name <- "piaya_therma"
kuenm_start(file.name = file_name)

# Variables with information to be used as arguments. Change "YOUR/DIRECTORY" by your actual directory.
occ_joint <- "piaya_therma_joint.csv"
occ_tra <- "piaya_therma_train.csv"
M_var_dir <- "M"
batch_cal <- "Candidate_models"
out_dir <- "Candidate_Models"
reg_mult <- c(seq(0.4, 2, 0.2), seq(2, 6, 0.5),8)
f_clas <- "all"
args <- NULL # e.g., "maximumbackground=20000" for increasing the number of pixels in the bacground or


# note that some arguments are fixed in the function and should not be changed
maxent_path <- "C:/Users/dprie/Documents/R/win-library/4.0/dismo/java"
wait <- FALSE
run <- TRUE
kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
          out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, args = args,
          maxent.path = maxent_path, wait = wait, run = run)
warnings()

############################################
######### The model evaluation #############
############################################
occ_test <- "piaya_therma_test.csv"
out_eval <- "Calibration_results"
threshold <- 5
rand_percent <- 50
iterations <- 500
kept <- TRUE
selection <- "OR_AICc"
paral_proc <- FALSE # make this true to perform pROC calculations in parallel, recommended
# only if a powerfull computer is used (see function's help)
# Note, some of the variables used here as arguments were already created for previous function
cal_eval <- kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, batch = batch_cal,
                        out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, iterations = iterations,
                        kept = kept, selection = selection, parallel.proc = paral_proc)


############################################
######### The model creation ###############
############################################
batch_fin <- "Final_models"
mod_dir <- "Final_Models"
rep_n <- 10
rep_type <- "Bootstrap"
jackknife <- FALSE
out_format <- "logistic"
project <- TRUE
G_var_dir <- "G"
ext_type <- "all"
write_mess <- TRUE
write_clamp <- TRUE
wait1 <- FALSE
run1 <- TRUE
args <- NULL # e.g., "maximumbackground=20000" for increasing the number of pixels in the bacground or
# "outputgrids=false" which avoids writing grids of replicated models and only writes the 
# summary of them (e.g., average, median, etc.) when rep.n > 1
# note that some arguments are fixed in the function and should not be changed
# Again, some of the variables used here as arguments were already created for previous functions
kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, batch = batch_fin, rep.n = rep_n,
          rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir, out.format = out_format, project = project,
          G.var.dir = G_var_dir, ext.type = ext_type, write.mess = write_mess, write.clamp = write_clamp, 
          maxent.path = maxent_path, args = args, wait = wait1, run = run1)


############################################
####### The model evaluation final #########
############################################
occ_ind <- "piaya_therma_ind.csv"
replicates <- TRUE
out_feval <- "Final_Models_evaluation"
# Most of the variables used here as arguments were already created for previous functions
fin_eval <- kuenm_feval(path = mod_dir, occ.joint = occ_joint, occ.ind = occ_ind, replicates = replicates,
                        out.eval = out_feval, threshold = threshold, rand.percent = rand_percent,
                        iterations = iterations, parallel.proc = paral_proc)

best <- read.csv("Calibration_results/best_candidate_models_OR_AICc.csv")
knitr::kable(best, caption = "Models selected based on significance, omission rates, and AICc, in that order.")



