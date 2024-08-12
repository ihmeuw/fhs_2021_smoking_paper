###############################################################################################################
## Purpose: Calculate Future Smoking Prevalence
###############################################################################################################

args <- commandArgs(trailingOnly = TRUE)

version        <- args[1]
diff_mortality <- as.logical(args[2])
target_yr      <- as.numeric(args[3])
reduction      <- as.numeric(args[4])
gbd_round_id   <- as.numeric(args[5])
reduction_type <- args[6]

print(version)
print(diff_mortality)
print(target_yr)
print(reduction)

fhs_path       <- "FILEPATH"
past_prev_path <- paste0("FILEPATH")
prev_path      <- paste0("FILEPATH")
mort_path      <- paste0("FILEPATH")

if (diff_mortality){
  prev_path <- "FILEPATH"
} else {
  prev_path <- "FILEPATH"
}

set.seed(34)

###### 1. Source Libraries #########################################################################################################
message(Sys.time())
message("Preparing Environment")
library(tidyr)
library(data.table)

###### 2. Source Functions #########################################################################################################
source(paste0("FILEPATH", "get_draws.R"))
source(paste0("FILEPATH", "get_population.R"))
source(paste0("FILEPATH", "get_location_metadata.R"))
source(paste0("FILEPATH", "utils.R"))

###### 2. Define some useful objects #########################################################################################################
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste0('task_id: ', task_id))

current_prev_path <- "FILEPATH"
former_prev_path  <- "FILEPATH"
  
r_id = 6
  
locs <- get_location_metadata(39, release_id = r_id)[level >= 3]$location_id

parameters <- expand_grid(location_id = locs,
                          sex_id = c(1,2)) %>%
  as.data.table()

temp <- fread(paste0("FILEPATH")) %>%
  .[location_id %in% locs]

l               <- parameters[task_id, location_id]
s               <- parameters[task_id, sex_id]
output_path     <- temp$output_path %>% unique()
current_prev_id <- temp$current_prev_id %>% unique()
former_prev_id  <- temp$former_prev_id %>% unique()
years           <- 2022
idvars          <- c("location_id", "year_id", "age_group_id", "sex_id")
drawvars        <- c(paste0("draw_", 0:999))


###### 6. Import prevalence #########################################################################################################

# Pull current smoking prevalence
message("Reading in Prevalence Draws")

current_prev <- fread("FILEPATH") %>%
  .[,version := "current_prev"] 
former_prev <- fread("FILEPATH") %>%
  .[,version := "former_prev"]

future_prev <- rbind(current_prev, former_prev) %>%
  .[scenario == 0] %>%
  .[year_id %in% years] %>%
  dcast(location_id + year_id + age_group_id + sex_id + draw ~ version, value.var = "value") %>%
  .[age_group_id %in% c(7,8), former_prev := 0] %>%
  .[(current_prev + former_prev) > 1, former_prev := 1 - current_prev] %>%
  .[,variable := paste0("draw_", draw)] %>%
  .[,c("draw") := NULL] %>%
  .[, current_age_group_id:=age_group_id]  

###### 7. Create specifications for future scenarios #########################################################################################################

if (reduction_type == "relative"){
  future_prev <- future_prev %>%
    .[, former_new_prev:=0] %>%
    .[, annual_reduction:=(current_prev * reduction)/(target_yr-2022)]
  
} else if (reduction_type == "absolute"){
  
  future_prev <- future_prev %>%
    .[, former_new_prev:=0] %>%
    .[,reduction_goal := current_prev - reduction] %>%
    .[reduction_goal <= 0, annual_reduction := 0] %>%
    .[reduction_goal > 0, annual_reduction := reduction_goal/(target_yr-2022)]
  
}

## Pull inputs for differential mortality
if (diff_mortality){
  message("Computing prevalence with differential mortality")
  
  ## Pull mortality
  pops <- get_population(location_id = l, sex_id = s, age_group_id = c(11:20, 30, 31, 32, 235), year_id = 2019,
                             gbd_round_id = 6, decomp_step = "step5") %>%
    .[,run_id := NULL]
  
  mort <- get_draws("cause_id", gbd_id=c(294), location_id= l, sex_id=s, age_group_id=c(11:20, 30, 31, 32, 235), year_id = 2019,
                    num_workers=5, gbd_round_id=6, decomp_step = "step5", measure_id=1, metric_id=1, source="codcorrect", version_id = 135) %>%
    .[,c("version_id", "metric_id", "measure_id") := NULL] %>%
    melt(id.vars = c("location_id", "year_id", "sex_id", "age_group_id", "cause_id")) %>%
    merge(pops, by = c("location_id", "year_id", "sex_id", "age_group_id"), all.x = T) %>%
    .[,val := value/population] %>%
    .[,.(age_group_id, variable, val)]
  
  mort[val > 1, val := 1]
  
} else {
  message("Computing prevalence without differential mortality")
}

keep_track <- NULL
for (t in 2022:2056) {
  print(t)
  
  current_yr <- t+1
  
  temp <- copy(future_prev[year_id == t]) %>%
    .[, year_id:=current_yr]
  
  if (current_yr %in% c(seq(2027, 2077, by = 5))) {
    temp[, current_age_group_id:=as.numeric(lapply(current_age_group_id, age_cohort))]
  }
  
  ## Increase former new prevalence by specified amount
  temp <- copy(temp) %>% 
    .[current_prev>annual_reduction, former_new_prev:=former_new_prev+annual_reduction] %>%
    .[current_prev<=annual_reduction, former_new_prev:=former_new_prev+current_prev]
  
  ## Reduce current prevalence by specified amount
  temp <- copy(temp) %>% 
    .[, current_prev:=current_prev-annual_reduction] %>%
    .[current_prev<0, current_prev:=0]
  
  if (diff_mortality){
    ## Account for differential mortality
    
    ## read in inputs
    current_mort <- fread(paste0(mort_path, "/result/current_",l,"_",s,".csv")) %>%
      .[,current_mort := exp(avg_rr)] %>%
      .[,avg_rr := NULL]
    
    former_mort <- fread(paste0(mort_path, "/result/former_",l,"_",s,".csv")) %>%
      .[,former_mort := exp(avg_rr)] %>%
      .[,avg_rr := NULL] 
    
    new_former_mort <- fread(paste0(mort_path, "/result/new_former_",l,"_",s,".csv")) %>%
      .[,new_former_mort := exp(avg_rr)] %>%
      .[,avg_rr := NULL]
    
    ## Fill in gaps
    for (a in c(7:10)){
      
      current_mort <- add_mort_rows(df = current_mort, a = a, type = "current")
      
      for (y in seq(2022, 2057, 5)){
        former_mort <- add_mort_rows(df = former_mort, a = a, y = y, type = "former")
        new_former_mort <- add_mort_rows(df = new_former_mort, a = a, y = y, type = "new_former")
        
      }
    }
    
    if (!current_yr %in% unique(former_mort$year_id)){
      former_mort <- former_mort %>%
        .[, former_mort := approx(year_id, former_mort, xout = current_yr)$y, by = c("variable", "age_group_id")] %>%
        .[,year_id := current_yr] %>%
        unique()
    } else {
      former_mort <- former_mort[year_id == current_yr]
    } 
    
    if (!current_yr %in% unique(new_former_mort$year_id)){
      new_former_mort <- new_former_mort %>%
        .[, new_former_mort := approx(year_id, new_former_mort, xout = current_yr)$y, by = c("variable", "age_group_id")] %>%
        .[,year_id := current_yr] %>%
        unique()
    } else {
      new_former_mort <- new_former_mort[year_id == current_yr]
    } 
    
    temp <- merge(temp, mort, by.x = c("current_age_group_id", "variable"), by.y = c("age_group_id", "variable"), all.x=T) %>%
      .[is.na(val), val:=0] %>%
      .[, never_smoker:=1-current_prev-former_prev-former_new_prev]
    
    temp <- temp %>% 
      merge(current_mort, by.x = c("location_id", "sex_id", "current_age_group_id", "variable"), by.y = c("location_id", "sex_id", "age_group_id", "variable"), all.x = T) %>%
      merge(former_mort, by.x = c("location_id", "year_id", "sex_id", "current_age_group_id", "variable"), by.y = c("location_id", "year_id","sex_id", "age_group_id", "variable"), all.x = T) %>%
      merge(new_former_mort, by.x = c("location_id", "year_id", "sex_id", "current_age_group_id", "variable"), by.y = c("location_id", "year_id", "sex_id", "age_group_id", "variable"), all.x = T)
    
    temp[, rhs := never_smoker+(former_mort*former_prev)+(new_former_mort*former_new_prev)+(current_mort*current_prev)] #never smoker base rate
    
    ## New prevalence
    temp <- temp %>%
      .[, ns_mort:=val/rhs] %>%
      .[, never_smoker:=never_smoker-(never_smoker*ns_mort)] %>%
      .[, former_prev:=former_prev-(former_prev*former_mort*ns_mort)] %>%
      .[, former_new_prev:=former_new_prev-(former_new_prev*new_former_mort*ns_mort)] %>%
      .[, current_prev:=current_prev-(current_prev*current_mort*ns_mort)]
    
    temp <- temp %>%
      .[never_smoker < 0, never_smoker:=0] %>%
      .[former_prev < 0, former_prev := 0] %>%
      .[former_new_prev < 0, former_new_prev:=0] %>%
      .[current_prev < 0, current_prev:=0]
    
    temp[, denom:=never_smoker+former_prev+former_new_prev+current_prev]
    temp[denom == 0, denom := 1e-6]
    
    temp <- temp %>%
      .[, never_smoker:=never_smoker/denom] %>%
      .[, former_prev:=former_prev/denom] %>%
      .[, former_new_prev:=former_new_prev/denom] %>%
      .[, current_prev:=current_prev/denom] 
    
    keep_track <- rbind(keep_track, temp)
    temp <- temp %>%
      .[, c("val", "never_smoker", "rhs", "ns_mort", "denom", "current_mort", "former_mort", "new_former_mort") := NULL] 
  }
  
    future_prev <- rbind(future_prev, temp, fill = T)
}

dir.create(prev_path, recursive = T)
fwrite(future_prev, paste0(prev_path, l, "_", s, ".csv"))
