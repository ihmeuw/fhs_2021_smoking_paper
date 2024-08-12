
library(tidyr)
library(data.table)
library(gtools)
library(magrittr)
library(parallel)
library(zoo)

source(paste0("FILEPATH", "get_location_metadata.R"))
source(paste0("FILEPATH", "utils.R"))

set.seed(124535)

#-- SET UP ARGS ----------------------------------------------------------------
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

args    <- commandArgs(trailingOnly = TRUE)
version        <- args[2]
n_draws        <- args[3] %>% as.numeric()
gbd_round_id   <- args[4] %>% as.numeric()
last_yr        <- args[5] %>% as.numeric()
reference      <- args[6] %>% as.logical()
ref_version    <- args[7]

cap <- !reference

if (gbd_round_id == 6){
  r_id <- 6
} else if (gbd_round_id == 7){
  r_id <- 9
}

paf_dir <- paste0("FILEPATH")
out_dir <- paste0("FILEPATH")

dir.create(out_dir)


### Other variables
locations   <- get_location_metadata(location_set_id = 39, release_id = r_id)[level >= 3,location_id]
sexes       <- c(1,2)

param_map   <- expand.grid(location = locations, sex = sexes) %>% as.data.table()
loc_id <- param_map[task_id, location]
sex_id <- param_map[task_id, sex]

## Causes that need to be run
input <- list.files(paf_dir) %>%
  as.data.table() %>%
  .[,`.` := gsub(".csv", "", `.`)] %>%
  separate(col = ".", into = c("location", "sex", "cause")) %>%
  .[,location := NULL] %>%
  .[sex == sex_id] %>%
  unique()

causes <- input$cause %>% unique

message("starting PAF processing for location_id ", loc_id)

for (cause_id in causes){
  #--PULL PAFs  ---------------------------------------------------------
  message("reading PAF draws")

  vars <- c("age_group_id", "sex_id", "location_id", "year_id", "cause_id", "cohort_start_age")
  
  paf <- rbindlist(lapply(paste0(paf_dir, loc_id, "_",sex_id,"_",cause_id,".csv"), fread))

  paf <- paf %>%
    .[measure_id == 4] %>%
    .[,c("rei_id", "measure_id") := NULL] %>%
    unique() %>% 
    melt(id.vars = vars,
         measure.vars = paste0("draw_", 0:(n_draws - 1)),
         variable.name = "draw", value.name = "paf") %>%
    .[, draw := as.numeric(gsub("draw_", "", draw))]
  
  if (!reference){
    ## Add young age groups
    
    start_age <- min(paf$cohort_start_age)
    cohort_ages_append <- c(6, 5, 4)
    
    for (c in cohort_ages_append){
      
      temp <- copy(paf[cohort_start_age == start_age])
      diff <- 7 - c
      
      temp <- temp %>%
        .[,cohort_start_age := c] %>%
        .[,age_group_id := age_group_id - diff ] %>%
        .[,paf := 0]
      
      paf <- rbind(paf, temp)
    }
    paf <- paf[age_group_id >= 11 & year_id <= 2054]
    rm(temp)
    
    ## Drop older age groups
    paf <- paf %>%
      .[age_group_id != 235, keep := cohort_start_age] %>%
      .[age_group_id == 235, keep := min(cohort_start_age), by = c("cause_id", "year_id")]
    paf <- paf[keep == cohort_start_age]
    paf[,c("keep")] <- NULL
    
  }

  ## Interpolate PAF
  have_yrs <- paf$year_id %>% unique
  if (reference){
    need_yrs <- seq(2019, last_yr, 1)
  } else {
    need_yrs <- seq(2022, last_yr, 1)
  }
  need_yrs <- setdiff(need_yrs, have_yrs)  
  
  paf_addyrs <- rbindlist(lapply(need_yrs, interp_yrs)) 
  paf <- rbind(paf, paf_addyrs)
  
  paf <- paf[year_id <= last_yr]
  
  
  if (!reference){
    
    ref_result <- fread("FILEPATH") %>%
        .[year_id %in% c(2019:2021)]
      paf <- rbind(ref_result, paf)
      table(paf$year_id, paf$age_group_id)
  
  }
  

## Cap alternative scenario to reference
if (cap){
  message("Capping!")
  ref_result <- fread(paste0("FILEPATH")) %>%
    .[,c("cohort_start_age") := NULL] %>%
    setnames(old = c("paf"), new = c("ref"))
  
  paf <- merge(paf, ref_result, by = c("location_id", "age_group_id", "year_id","sex_id", "cause_id","draw"), all.x = T) %>%
    .[,flag := 0] %>%
    .[paf > ref, flag := 1] 
  
  paf[flag == 1, end_yr := (year_id)]
  paf[, end_yr := max(end_yr, na.rm = T), by = c("location_id", "age_group_id","sex_id", "cause_id","draw")]
  
  paf_issue <- copy(paf) %>%
    .[flag == 1]
  
  paf[year_id <= end_yr, paf := ref]
  
  paf[,c("flag", "end_yr", "ref")] <- NULL
  
}

dir.create(paste0(out_dir, "/paf/"))
fwrite(paf, paste0(out_dir, "/paf/",loc_id, "_", sex_id,"_", cause_id,".csv"))

}
