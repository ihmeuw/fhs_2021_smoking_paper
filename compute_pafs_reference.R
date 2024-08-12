###############################################################################################################
## Purpose: Calculate Smoking PAFs
###############################################################################################################

set.seed(34)

args <- commandArgs(trailingOnly = TRUE)

version        <- args[1]
gbd_round_id   <- as.integer(args[2])

scenario_id    <- 0
draws          <- c(0:999)
end_year       <- 2100

prev_path         <- paste0("FILEPATH")
paf_path          <- paste0("FILEPATH")
dir.create(paf_path)

estimation_yrs <- c(2019, 2020, 2022, 2023, seq(2027, end_year, by = 5), end_year)
  
current_prev_path <- paste0("FILEPATH")
former_prev_path  <- paste0("FILEPATH")
  
r_id = 6
  

###### 1. Source Libraries #########################################################################################################
message(Sys.time())
message("Preparing Environment")

library(tidyr)
library(dplyr)
library(data.table)
library(magrittr)


###### 2. Source Functions #########################################################################################################
setwd("FILEPATH")
source(paste0("FILEPATH", "exposure_simulation_functions.R"))
source(paste0("FILEPATH", "get_location_metadata.R"))
source(paste0("FILEPATH", "get_age_metadata.R"))
source(paste0("FILEPATH", "save.R"))
source(paste0("FILEPATH", "utils.R"))

###### 3. Setting up for array jobs #########################################################################################################
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste0('task_id: ', task_id))

locs <- get_location_metadata(39, release_id = r_id)[level >= 3]$location_id

fhs_path <- "FILEPATH"
ref <- fread("FILEPATH")
temp <- fread("FILEPATH") %>%
  .[location_id %in% locs]

parameters <- NULL
for (c in unique(ref$cause_id)) {
  add_df <- temp[, cause_id:=c]
  parameters <- rbind(parameters, add_df, fill = T)
}
parameters <- parameters[!((cause_id %in% c(429, 432) & sex_id == 1) | (cause_id == 438 & sex_id == 2))]

l           <- parameters[task_id, location_id]
s           <- parameters[task_id, sex_id]
c           <- parameters[task_id, cause_id]
output_path <- parameters[task_id, output_path]


###### 4. Define some useful objects #########################################################################################################
idvars   <- c("location_id", "year_id", "age_group_id", "sex_id")
drawvars <- c(paste0("draw_", draws))
ages_og <- c(11:20, 30, 31, 32, 235)

age_map <- get_age_metadata(age_group_set_id = 19)

ref[, lag:=0] 
if (s == 1) { ref<-ref[!(cause_id%in%c(429, 432))]}
if (s == 2) { ref<-ref[!(cause_id%in%c(438))]}

###### 5. Import exposure distributions #########################################################################################################
# Read in the simulated histories
message("Loading simulated histories")

attach("FILEPATH", name = "exposure")

for (y in c(2019, 2020, 2022)){
  for (a in c(9:20, 30:32, 235)){
    print(paste0("loading ", "amt_",a,"_",s,"_", y))
    assign(paste0("py_",a,"_",s,"_", y), get(paste0("py_",a,"_",s,"_", y)))
    assign(paste0("amt_",a,"_",s,"_", y), get(paste0("amt_",a,"_",s,"_", y)))
    assign(paste0("cess_",a,"_",s,"_", y), get(paste0("cess_",a,"_",s,"_", y)))
  }}

detach(exposure)


###### 6. Import prevalence #########################################################################################################
current_prev <- fread(paste0(current_prev_path, l, "_", s, ".csv")) %>%
  .[scenario == scenario_id] %>%
  .[,c("scenario", "V1") := NULL] %>%
  setnames(old = "value", new = "current_prev")

former_prev <- fread(paste0(former_prev_path, l, "_", s, ".csv"))%>%
  .[scenario == scenario_id] %>%
  .[,c("scenario", "V1") := NULL] %>%
  setnames(old = "value", new = "former_prev")

future_prev <- merge(current_prev, former_prev, by = c("location_id", "age_group_id", "sex_id", "year_id", "draw"), all = T) %>%
  .[year_id %in% estimation_yrs] %>%
  .[age_group_id %in% c(7,8), former_prev := 0] %>%
  .[(current_prev + former_prev) > 1, former_prev := 1 - current_prev] %>% 
  setnames(old = "draw", new = "variable") %>%
  .[,variable := paste0("draw_", variable)]  %>%
  .[variable %in% drawvars]

rm(current_prev, former_prev)

#################################
## Calculate the PAFs
#################################

# Set some cause-specific parameters
ref_use <- ref[cause_id == c]

if (!c %in% c(878, 923)){
  cause_c <- ref_use$cause
} else {
  cause_c <- "fractures"
}

path_current <- "FILEPATH"
path_former  <- "FILEPATH"
exp_def      <- ref_use$current_exp_def
exp_max      <- ref_use$max_exp

# read in RRs
if (exp_def == "prevalence") {
  rr_current<-fread(paste0(path_current)) %>%
    .[cause_id==c] %>%
    .[, draw := paste0("draw_", draw)] %>%
    setnames("draw", "variable") %>%
    .[variable %in% drawvars]
  
} else {
  # Read in the RR draws for that cause (current and former), age, sex
  rr_current <- fread(paste0(path_current)) %>%
    .[draw %in% draws]
  
  # Cut at RR max
  if (!(is.na(exp_max))) {
    rr_current<-rr_current[exposure<=exp_max]
  } 
  
  rr_former <- fread(paste0(path_former)) %>%
    .[exposure < 60]  %>% 
    .[draw %in% draws]
  
  temp_rr_former <- rr_former[exposure == 0] %>%
    .[, exposure := 60] %>%
    .[, rr := 1]
  rr_former <- rbind(rr_former, temp_rr_former)
  
  if (c != 544) {rr_former[rr<1, rr:=1]} 
  
}

# Initialize the outframe
out<-NULL

# Loop through the causes to calculate PAFs
for (future_year in estimation_yrs) {
  print(future_year)
  
  if (future_year %in% c(2019:2022)){
    years <- future_year
  } else {
    years <- 2022
  }
  
  for (cohort_a in ages_og) {

      if (exp_def == "prevalence") {
        
        pafs_temp<-merge(rr_current, future_prev, by = c("variable", "age_group_id", "sex_id")) %>%
          .[age_group_id == cohort_a & year_id == future_year] %>%
          .[, variable:=as.numeric(gsub(variable, pattern = "draw_", replacement = ""))] %>%
          setnames("variable", "draw") %>%
          .[, paf:=((1-current_prev)+(current_prev*rr)-1)/((1-current_prev)+(current_prev*rr))] %>%
          .[,.(draw, sex_id, cause_id, location_id, year_id, age_group_id, paf)]
        
      } else {
        
        # Flatten current smoker RR into a list
        rr_current_list<-split(rr_current[age_group_id==cohort_a & sex_id==s], by = "draw")
        # Interpolate and define a RR function
        rr_current_list<-lapply(rr_current_list, function(x) approxfun(x$exposure, x$rr, method = "linear", rule =2))
        
        # Calculate the current smoker exposure integral (exposure-weighted RR)
        current_smoker_integral<-mapply(calc_integral, get(paste0(exp_def, "_", cohort_a, "_", s, "_", years))[draws+1], rr_current_list, 0, cap_exposure(exp_def))
        
        # Calculate the former smoking risk curve based on the current smoking integral
        new_max_temp<-data.table(new_max=current_smoker_integral, draw=draws)
        rr_former_temp<-rr_former[age_group_id==cohort_a & sex_id==s] %>%
          .[, max:=max(rr), by = c("draw", "age_group_id", "sex_id")] %>%
          merge(new_max_temp, by = "draw")
        
        if (c==544) {
          rr_former_temp$rr_new<-mapply(transform_former, orig_val = rr_former_temp$rr, new_max = 1, new_min = rr_former_temp$new_max)
        } else {
          rr_former_temp$rr_new<-mapply(transform_former, orig_val = rr_former_temp$rr, new_max = rr_former_temp$new_max, orig_max = rr_former_temp$max)
        }
        rr_former_list<-split(rr_former_temp, by = "draw")
        rr_former_list<-lapply(rr_former_list, function(x) approxfun(x$exposure, x$rr_new, method = "linear", rule =2))
        
        # Calculate the former smoking intergral (exposure-weighted RR)
        former_smoker_integral<-mapply(calc_integral, get(paste0("cess_", cohort_a, "_", s, "_", years))[draws+1], rr_former_list, 0, age_map[age_group_id==cohort_a, age_group_years_end])
        
        # Calculate the PAF
        future_prev <- future_prev[order(year_id, age_group_id, variable)] 
        
        current_prev_vec<-as.vector(future_prev$current_prev[future_prev$age_group_id==cohort_a & future_prev$sex_id==s & future_prev$year_id == future_year])
        former_prev_vec<-as.vector(future_prev$former_prev[future_prev$age_group_id==cohort_a & future_prev$sex_id==s & future_prev$year_id == future_year])
        
        never_smoker_burden      <- (1-(current_prev_vec+former_prev_vec))
        current_smoker_burden    <- current_prev_vec*current_smoker_integral
        former_smoker_burden     <- former_prev_vec*former_smoker_integral
        
        pafs_temp<-data.table(cause_id = c, location_id = l, year_id = future_year, age_group_id = cohort_a, 
                              sex_id = s, draw = draws, 
                              former_prev = former_prev_vec,
                              former_integral = former_smoker_integral,
                              current_prev = current_prev_vec,
                              current_integral = current_smoker_integral,
                              never_smoker_burden = never_smoker_burden,
                              current_smoker_burden = current_smoker_burden,
                              former_smoker_burden = former_smoker_burden)
        
        pafs_temp[, paf := (never_smoker_burden + current_smoker_burden + former_smoker_burden - 1) /
                    (never_smoker_burden + current_smoker_burden + former_smoker_burden)]
        
      }
      out<-rbind(out, pafs_temp, fill = T)  
    }
}

out <- expand_causes_parent(out)

# Format for save
out<-out[, c("mortality", "morbidity"):=1]
out<-out[, rei_id:=99]
sex_id<-s
location_id<-l
cause_id <- c

out <- merge(out, age_map[,.(age_group_id, age_group_years_start)], by = "age_group_id") %>%
  .[, cohort_start_age:=year_id-age_group_years_start] %>%
  .[, age_group_years_start:=NULL]

dir.create(paste0(paf_path,"/pafs_ss_annual/"))
dir.create(paste0(paf_path,"/paf_summaries_ss_annual/"))

save_paf(dt = out, rid=99, rei="smoking_direct", n_draws=length(draws), out_dir=paste0(paf_path,"/pafs_ss_annual/"))

# Save summaries
out <- out %>%
  .[,cohort_start_age := NULL] %>%
  .[, mean:=mean(paf), by = c(idvars, "cause_id")] %>%
  .[, lower:=quantile(paf, 0.025), by = c(idvars, "cause_id")] %>%
  .[, upper:=quantile(paf, 0.975), by = c(idvars, "cause_id")]
out<-unique(out[,.(location_id, year_id, age_group_id, sex_id, cause_id, mean, lower, upper)])

write.csv(out, paste0(paf_path,"/paf_summaries_ss_annual/", l, "_", s, "_", c,".csv"), na = "", row.names = F)

