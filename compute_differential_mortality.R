###############################################################################################################
## Purpose: Calculate Differential Mortality
###############################################################################################################

args <- commandArgs(trailingOnly = TRUE)

version      <- args[1]
gbd_round_id <- as.numeric(args[2])

input_path <- "FILEPATH"


##### 1. Call libraries -------------------------------------------------------------------------------
library(data.table)
library(magrittr)
library(tidyr)
library(parallel)

##### 2. Call Functions -------------------------------------------------------------------------------
source(paste0("FILEPATH", "get_draws.R"))
source(paste0("FILEPATH", "get_population.R"))
source(paste0("FILEPATH", "get_location_metadata.R"))
source(paste0("FILEPATH", "utils.R"))

##### 3. Define variables -------------------------------------------------------------------------------
if (gbd_round_id == 6){
  r_id = 6
} else if (gbd_round_id == 7){
  r_id = 9
}

locs <- get_location_metadata(39, release_id = r_id)[level >= 3]$location_id
start_yr <- 2022
y        <- 2019 
ages     <- c(9:20, 30, 31, 32, 235)
sexes    <- c(1,2)

demo_cols  <- c("location_id", "year_id", "age_group_id", "sex_id")
other_cols <- c(demo_cols, "cause_id")

param_map <- expand_grid(location_id = locs, 
                        sex_id = sexes) %>% as.data.table()

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
l <- param_map[task_id, location_id]
s <- param_map[task_id, sex_id]

print(l)
print(s)

##### 4. Read in mortality data -------------------------------------------------------------------------------

fhs_path <- "FILEPATH"
ref <- fread(paste0(fhs_path, "/paf_calc_ref.csv"))[,c("cause_id", "cause")]
if (s == 1) { ref<-ref[!(cause_id%in%c(429, 432))]}
if (s == 2) { ref<-ref[!(cause_id%in%c(438))]}

ref <- expand_causes_parent(ref)

causes <- unique(ref$cause_id)

# Pull number of deaths
pops <- get_population(location_id = l, sex_id = s, age_group_id = ages, year_id = 2019,
                       gbd_round_id = 6, decomp_step = "step5", run_id = 192) %>%
  .[,run_id := NULL]

mort <- get_draws("cause_id", gbd_id=c(294, causes), location_id= l, sex_id=s, age_group_id=ages, year_id = y,
                  num_workers=5, gbd_round_id=6, decomp_step = "step5", measure_id=1, metric_id=1, source="codcorrect", version_id = 135) %>%
  .[,c("version_id", "metric_id", "measure_id") := NULL] %>%
  melt(id.vars = other_cols, variable.name = "draw")

max_deaths <- copy(mort[cause_id == 294]) %>%
  .[,max_value := max(value), by = c("location_id", "year_id", "sex_id", "age_group_id")] %>%
  .[,c("cause_id","draw", "value") := NULL] %>%
  unique()

pops <- merge(pops, max_deaths, by = c("location_id", "year_id", "sex_id", "age_group_id"), all.x = T)
pops[max_value > population, population := max_value]
pops$max_value <- NULL

done <- mort$cause_id %>% unique

# compute number of non-smoking deaths
mort_smoking <- copy(mort[cause_id != 294]) %>%
  .[,smoking := sum(value), by = c("draw", demo_cols)] %>%
  .[,cause_id := NULL] %>%
  .[,value := NULL] %>%
  unique()

mort_nonsmoking <- copy(mort[cause_id == 294]) %>%
  merge(mort_smoking, by = c("draw", demo_cols)) %>%
  .[,non_smoking := value - smoking] %>%
  .[,cause_id := 999] %>%
  .[,value := NULL] %>%
  .[,smoking := NULL] %>%
  setnames(old = "non_smoking", new = "value")

rates <- rbind(mort[cause_id != 294], mort_nonsmoking)
rm(mort_smoking, mort)

# transform into rate space
rates <- merge(rates, pops, by = demo_cols, all.x = T) %>%
  .[,death_rate := value/population] %>%
  .[,draw := gsub("draw_", "", draw)] %>%
  .[,draw := as.numeric(draw)]
rates[,c("value", "population", "year_id")] <- NULL


##### 5. Read in weighted RR data -------------------------------------------------------------------------------
current_path    <- paste0(input_path, "/current/")
former_path     <- paste0(input_path, "/former/")
new_former_path <- paste0(input_path, "/new_former/")

current <- compute_diff_mortality(path = current_path, type = "current", l = l, s = s)
former <- compute_diff_mortality(path = former_path, type = "former", l = l, s = s)
new_former <- compute_diff_mortality(path = new_former_path, type = "former_new", l = l, s = s)

save_path <- paste0(input_path, "/result/")

dir.create(save_path)
fwrite(current, paste0(save_path,"/current_",l,"_",s,".csv"))
fwrite(former, paste0(save_path,"/former_",l,"_",s,".csv"))
fwrite(new_former, paste0(save_path,"/new_former_",l,"_",s,".csv"))

