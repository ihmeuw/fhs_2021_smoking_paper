###################################################################################################################################
## Purpose: Calculate Future Smoking PAFs
###################################################################################################################################

###### 1. Source Libraries ########################################################################################################
library(data.table)
library(tidyr)
library(dplyr)
library(magrittr)


###### 2. Source Functions ########################################################################################################
setwd("FILEPATH")
source(paste0("FILEPATH", "exposure_simulation_functions.R"))
source(paste0("FILEPATH", "get_location_metadata.R"))
source(paste0("FILEPATH", "get_age_metadata.R"))
source(paste0("FILEPATH", "save.R"))
source(paste0("FILEPATH", "utils.R"))


###### 3. Define Variables ########################################################################################################

set.seed(34)
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
args    <- commandArgs(trailingOnly = TRUE)

version        <- args[1]
diff_mortality <- as.logical(args[2])
target_yr      <- as.numeric(args[3])
reduction      <- as.numeric(args[4])
gbd_round_id   <- as.numeric(args[5])

if (gbd_round_id == 6){
  r_id = 6
  
} else if (gbd_round_id == 7){
  r_id = 9
}

locs <- get_location_metadata(39, release_id = r_id)[level >= 3]$location_id
age_map <- get_age_metadata(age_group_set_id = 19)

fhs_path <- "FILEPATH"
ref <- fread("FILEPATH")
temp <- fread("FILEPATH") %>%
  .[location_id %in% locs]

## Paths where differential mortality and PAFs saved, respectively
mort_path <- paste0("FILEPATH")
paf_path  <-  paste0("FILEPATH")

parameters <- NULL
for (c in unique(ref$cause_id)) {
  add_df <- temp[, cause_id:=c]
  parameters <- rbind(parameters, add_df, fill = T)
}
parameters <- parameters[!((cause_id %in% c(429, 432) & sex_id == 1) | (cause_id == 438 & sex_id == 2))]

l               <- parameters[task_id, location_id]
s               <- parameters[task_id, sex_id]
c               <- parameters[task_id, cause_id]
exp_path        <- parameters[task_id, output_path]


###### 4. Define some useful objects #########################################################################################################

if (diff_mortality){
  prev_path <- "FILEPATH"
} else {
  prev_path <- "FILEPATH"
}

idvars   <-c("location_id", "year_id", "age_group_id", "sex_id")
drawvars <-c(paste0("draw_", 0:999))
ages_og  <-c(11:20, 30, 31, 32, 235)
years    <- 2022

ref[, lag:=0] 
if (s == 1) { ref<-ref[!(cause_id%in%c(429, 432))]}
if (s == 2) { ref<-ref[!(cause_id%in%c(438))]}

###### 5. Import exposure distributions #########################################################################################################
# Read in the simulated histories
message("Loading simulated histories")

attach("FILEPATH", name = "exposure")

for (y in c(2022)){
  for (a in c(9:20, 30:32, 235)){
    assign(paste0("py_",a,"_",s,"_", y), get(paste0("py_",a,"_",s,"_", y)))
    assign(paste0("amt_",a,"_",s,"_", y), get(paste0("amt_",a,"_",s,"_", y)))
    assign(paste0("cess_",a,"_",s,"_", y), get(paste0("cess_",a,"_",s,"_", y)))
  }}

detach(exposure)


###### 6. Import prevalence #########################################################################################################

future_prev <- fread(paste0(prev_path, l, "_", s,".csv"))

## Compute "new" new former smokers over some time period
future_prev_pastyr <- copy(future_prev) %>%
  .[,year_id := year_id - 5] %>%
  .[,c("current_prev", "former_prev", "annual_reduction", "current_age_group_id") := NULL] %>%
  setnames("former_new_prev", "former_new_prev_futureyr") 

future_prev <- future_prev %>%
  merge(future_prev_pastyr, by = c("location_id", "year_id", "age_group_id", "sex_id", "variable"), all.x = T) %>%
  .[!is.na(former_new_prev_futureyr)] %>% 
  .[,new_newformer_prev := former_new_prev_futureyr - former_new_prev]


###### 7. Import RRs #########################################################################################################

# Set some cause-specific parameters
cvd_outcomes <- ref[cause %in% c("stroke", "afib_and_flutter", "aortic_aneurism", "ihd", "peripheral_artery_disease")]$cause_id
ref_use <- ref[cause_id == c]

if (!c %in% c(878, 923)){
  cause_c <- ref_use$cause
} else {
  cause_c <- "fractures"
}

path_current <- "FILEPATH"
path_former  <- "FILEPATH"
exp_def      <- ref_use$current_exp_def ## pack-year or cig-day cause
exp_max      <- ref_use$max_exp ## what is the max exposure we want to evaluate RR at for the given cause

if (!(is.na(exp_max))) {
  exp_max <- Inf
}

if (exp_def == "prevalence") {
  
  rr_current<-fread(path_current) %>%
    .[cause_id==c] %>%
    .[, draw := paste0("draw_", draw)] %>%
    setnames("draw", "variable")
  
} else {
  # Read in the RR draws for that cause (current and former), age, sex
  rr_current <- fread(paste0(path_current)) %>%
    .[exposure<=exp_max & sex_id == s]
  
  rr_former  <- fread(paste0(path_former)) %>%
    .[exposure < 60  & sex_id == s] 
  
  temp_rr_former <- copy(rr_former[exposure==0]) %>%
    .[, exposure:=60] %>%
    .[, rr:=1]
  
  rr_former <- rbind(rr_former, temp_rr_former)
  
  if (c != 544) {rr_former[rr<1, rr:=1]} 
  
  # Create age-specific former smoker scaling factors 
  age_cessation <- c(9:10, ages_og)
  
  if (c %in% cvd_outcomes){
    
    for (a in age_cessation){
      current_age <- age_cessation[age_cessation >= a] #possible future ages that a quitter at age "a" can be
      
      for (ca in current_age){
        assign(paste0("current_smoker_integral_",a, "_", ca), compute_exposure_weighted_rr(distribution_age = a, rr_age = ca))
        assign(paste0("rr_former_list_", a, "_", ca), compute_scaled_former_rr(scaling_factor_name = paste0("current_smoker_integral_", a, "_", ca),
                                                                               rr_age = ca))
      }
    }
    
  } else {
    
    for (a in age_cessation){
      assign(paste0("current_smoker_integral_",a), compute_exposure_weighted_rr(rr_age = a, distribution_age = a))
      assign(paste0("rr_former_list_", a), compute_scaled_former_rr(scaling_factor_name = paste0("current_smoker_integral_",a),
                                                                    rr_age = a))
      
    }
  }
}


####################################################################################################################################
## 8. Calculate the PAFs
####################################################################################################################################
# Initialize 
out                                   <- NULL #This will append all of the PAFs together

current_smoker_integral_save_total    <- NULL #These three objects will append all the inputs for differential mortality together
former_smoker_integral_save_total     <- NULL
new_former_smoker_integral_save_total <- NULL

## Start PAF computation: we want to observe the same cohorts into the future, so we estimate PAFs at 5 year intervals
# Therefore, we loop over every 5th future year and every cohort, based on their age in 2022
cohort_start_year <- 2022
cohort_end_year <- 2052

for (future_year in seq(cohort_start_year, cohort_end_year, by = 5)) {
  print(future_year)
  
  for (cohort_a in unique(future_prev$age_group_id)) {
    print(paste0("Cohort: ", cohort_a))
    
    ## Find the birth-cohort's current age in the current "future_year"
    current_age_raw <- age_map[age_group_id==cohort_a]$age_group_years_start + (future_year-cohort_start_year) 
    current_age <- age_map[age_group_years_start==current_age_raw]$age_group_id
    if (current_age_raw >= 100) {current_age <- 300} 
    
    if(current_age>=11 & current_age_raw < 100) {
      message(paste0("Calculating PAFs for Cause ID: ", c))
      
      if (exp_def == "prevalence") {
        
        pafs_temp <- merge(rr_current, future_prev, by = c("variable", "age_group_id", "sex_id")) %>%
          .[age_group_id == cohort_a & year_id == future_year] %>%
          .[, variable:=as.numeric(gsub(variable, pattern = "draw_", replacement = ""))] %>%
          setnames("variable", "draw") %>%
          .[, paf:=((1-current_prev)+(current_prev*rr)-1)/((1-current_prev)+(current_prev*rr))]
        
        current_smoker_integral    <- pafs_temp$rr
        former_smoker_integral     <- 0
        new_former_smoker_integral <- 0
        
        current_smoker_integral_save <- format_weighted_rr(current_smoker_integral)
        former_smoker_integral_save <- format_weighted_rr(former_smoker_integral)
        new_former_smoker_integral_save <- format_weighted_rr(new_former_smoker_integral)
        
        current_smoker_integral_save_total <- rbind(current_smoker_integral_save_total, current_smoker_integral_save)
        former_smoker_integral_save_total <- rbind(former_smoker_integral_save_total, former_smoker_integral_save)
        new_former_smoker_integral_save_total <- rbind(new_former_smoker_integral_save_total, new_former_smoker_integral_save)
        
        pafs_temp <- pafs_temp[,.(draw, sex_id, cause_id, location_id, year_id, paf)]
        
        pafs_temp <- pafs_temp %>%
          .[, year_id:=future_year] %>%
          .[, age_group_id:=current_age] %>%
          .[, cohort_start_age:=cohort_a]
        
      } else {
        
        ### 1. Calculate burden for current smokers
        if (c %in% cvd_outcomes){
          current_smoker_integral <- get(paste0("current_smoker_integral_", current_age, "_", current_age))
        } else {
          current_smoker_integral <- get(paste0("current_smoker_integral_", current_age))
        }
        
        ### 2. Calculate burden for "existing" former smokers
        
        # Calculate the former smoking integral (exposure-weighted RR)
        if (cohort_a >= 9) {
          
          if (c %in% cvd_outcomes){
            rr_former_list <- get(paste0("rr_former_list_", cohort_a, "_", current_age))
          } else {
            rr_former_list <- get(paste0("rr_former_list_", cohort_a))
          }
          
          # take the cessation distribution this cohort had in the last year they could quit (2022), and shift it as the future_year changes
          shift_cess <- lapply(get(paste0("cess_", cohort_a, "_", s, "_", years)), FUN = new_cess, t = future_year-cohort_start_year) 
          former_smoker_integral<-mapply(calc_integral, exp = shift_cess, rr = rr_former_list, minval = 0, maxval = age_map[age_group_id==current_age, age_group_years_end])
        } else {
          former_smoker_integral <- 1
        }
        
        ### 3. Calculate burden for "new" former smokers
        if (future_year==2022 | (cohort_a < 9 & target_yr == 2023) | (current_age < 9 & target_yr != 2023)) { 
          new_former_smoker_integral <- 1
          
        } else if (target_yr == 2023 & reduction == 1){

          new_former_smoker_integral<-mapply(calc_integral_new_former, rr = rr_former_list, 
                                             unif_start = future_year-target_yr, unif_max = future_year-target_yr+1,  
                                             minval = future_year-target_yr, maxval = future_year-target_yr+1)
          
        } else {
          
          ## First, need to determine the first year that smoking prevalence hits 0
          temp <- copy(future_prev[age_group_id == cohort_a] )%>%
            .[,min_current_prev := min(current_prev), by = variable] %>%
            .[min_current_prev == 0 & current_prev == 0, min_year := min(year_id), by = variable] %>%
            .[min_current_prev > 0, min_year := 2052, by = variable]
          
          min_prev_year <- temp[year_id == max(temp$year_id)]$min_year 
          
          #track sub-cohorts of quitters
          need_yrs_total <- seq(2022, (future_year), by = 5)
          
          new_former_smoker_integral <- NULL
          
          for (d in c(1:1000)){
            
            ## In case differential mortality causes smoking prevalence to be eliminated before the year currently being estimated, need 0 beyond
            last_highest <- min_prev_year[d]
            
            if (last_highest == 2022){ 
              
              new_former_smoker_integral_d <- 1
              
            } else {
              
              need_yrs <- seq(2022, min(future_year-5, last_highest-1), by = 5) 
              zero_yrs <- setdiff(need_yrs_total, need_yrs) 
              
              for (ny in need_yrs_total){
                age_2022 <- cohort_a
                int_increase <- (ny-2022)/5
                
                ages_list <- c(7:20, 30:32, 235)
                position <- match(age_2022,ages_list)
                
                age_during_ny <- ages_list[position + int_increase]
                
                
                if (ny %in% zero_yrs | age_during_ny < 9) {
                  assign(paste0("new_former_smoker_integral_",ny), 1)
                } else {
                  
                  if (ny+5 <= last_highest) {
                    start_range <- future_year - (ny + 5)
                    end_range <- future_year - (ny)
                  } else { 
                    start_range <- future_year - last_highest 
                    end_range <- future_year - ny
                  }
                  
                  if (c %in% cvd_outcomes){
                    assign(paste0("new_former_smoker_integral_", ny),
                           mapply(calc_integral_new_former, unif_start = start_range, unif_max = end_range, rr = get(paste0("rr_former_list_", age_during_ny,"_",current_age))[d], minval = start_range, maxval = end_range) 
                    )  
                  } else {
                    assign(paste0("new_former_smoker_integral_", ny),
                           mapply(calc_integral_new_former, unif_start = start_range, unif_max = end_range, rr = get(paste0("rr_former_list_", age_during_ny))[d], minval = start_range, maxval = end_range) 
                    ) 
                  }
                }
              }
              
              draw_cols <- unique(future_prev$variable)
              
              demo_sub <- copy(future_prev[age_group_id==cohort_a & sex_id==s & year_id %in% need_yrs & variable == draw_cols[d]]) %>% 
                .[new_newformer_prev <= 0, new_newformer_prev := 1e-6] %>% 
                .[,total_new_former := sum(new_newformer_prev), by = c("location_id", "sex_id", "age_group_id", "variable")] %>%
                .[,prop_new_former := new_newformer_prev/total_new_former]
              
              for (ny in need_yrs){
                
                prev_vec <- as.vector(demo_sub[year_id == ny]$prop_new_former)
                integral_vec <- get(paste0("new_former_smoker_integral_", ny))
                
                assign(paste0("weighted_new_former_smoker_integral_",ny),
                       prev_vec * integral_vec)
                
              }
              
              weighted_interals_list <- paste0("weighted_new_former_smoker_integral_",need_yrs)
              
              new_former_smoker_integral_d <- 0
              for (w in weighted_interals_list){
                new_former_smoker_integral_d <- new_former_smoker_integral_d + get(w)
              }
            }
            
            new_former_smoker_integral <- c(new_former_smoker_integral, new_former_smoker_integral_d)
            
          }
        } 
        
        if (!diff_mortality){
          current_smoker_integral_save <- format_weighted_rr(current_smoker_integral)
          former_smoker_integral_save <- format_weighted_rr(former_smoker_integral)
          new_former_smoker_integral_save <- format_weighted_rr(new_former_smoker_integral)
          
          current_smoker_integral_save_total <- rbind(current_smoker_integral_save_total, current_smoker_integral_save)
          former_smoker_integral_save_total <- rbind(former_smoker_integral_save_total, former_smoker_integral_save)
          new_former_smoker_integral_save_total <- rbind(new_former_smoker_integral_save_total, new_former_smoker_integral_save)
        }
        
        future_prev <- future_prev[order(year_id, age_group_id, variable)]
        
        # Calculate the PAF
        demo_sub <- copy(future_prev[age_group_id==cohort_a & sex_id==s & year_id == future_year])
        demo_sub[,former_prev := as.numeric(former_prev)]
        
        current_prev_vec<-as.vector(demo_sub$current_prev)
        former_prev_vec<-as.vector(demo_sub$former_prev)
        former_prev_new_vec<-as.vector(demo_sub$former_new_prev)
        
        never_smoker_burden      <- (1-(current_prev_vec+former_prev_vec+former_prev_new_vec))
        current_smoker_burden    <- current_prev_vec*current_smoker_integral
        former_smoker_burden     <- former_prev_vec*former_smoker_integral
        new_former_smoker_burden <- former_prev_new_vec*new_former_smoker_integral
        
        pafs_temp<-data.table(cause_id = c, location_id = l, year_id = future_year, age_group_id = current_age, 
                              cohort_start_age = cohort_a, sex_id = s, draw = 0:999, 
                              never_smoker_burden = never_smoker_burden,
                              current_smoker_burden = current_smoker_burden,
                              former_smoker_burden = former_smoker_burden,
                              new_former_smoker_burden = new_former_smoker_burden)
        
        pafs_temp[, paf := (never_smoker_burden + current_smoker_burden + former_smoker_burden + new_former_smoker_burden - 1) /
                    (never_smoker_burden + current_smoker_burden + former_smoker_burden + new_former_smoker_burden)]
        
      }
      
      out<-rbind(out, pafs_temp, fill = T)  
    }
  }
}


if (!diff_mortality){
  dir.create(paste0(mort_path,"/current/"), recursive = T)
  dir.create(paste0(mort_path,"/former/"), recursive = T)
  dir.create(paste0(mort_path,"/new_former/"), recursive = T)
  fwrite(current_smoker_integral_save_total, paste0(mort_path,"/current/",l,"_",s,"_",c,".csv"))
  fwrite(former_smoker_integral_save_total, paste0(mort_path,"/former/",l,"_",s,"_",c,".csv"))
  fwrite(new_former_smoker_integral_save_total, paste0(mort_path,"/new_former/",l,"_",s,"_",c,".csv"))
  
} else {

  out <- expand_causes_parent(out)
  
  # Format for save
  out <- out %>%
    .[, c("mortality", "morbidity"):=1] %>%
    .[, rei_id:=99]
  
  sex_id<-s
  location_id<-l
  cause_id <- c
  
  dir.create(paste0(paf_path,"/pafs_ss_annual/"), recursive = T)
  dir.create(paste0(paf_path,"/paf_summaries_ss_annual/"))
  
  save_paf(dt = out, rid=99, rei="smoking_direct", n_draws=1000, out_dir=paste0(paf_path,"/pafs_ss_annual/"))
  
  # Save summaries
  out <- out %>%
    .[, mean:=mean(paf), by = c(idvars, "cause_id", "cohort_start_age")] %>%
    .[, lower:=quantile(paf, 0.025), by = c(idvars, "cause_id", "cohort_start_age")] %>%
    .[, upper:=quantile(paf, 0.975), by = c(idvars, "cause_id", "cohort_start_age")] %>%
    .[,.(location_id, year_id, age_group_id, cohort_start_age, sex_id, cause_id, mean, lower, upper)] %>%
    unique()
  
  write.csv(out, paste0(paf_path,"/paf_summaries_ss_annual/", l, "_", s, "_", c,".csv"), na = "", row.names = F)
  
}

