
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
version       <- args[1]
n_draws       <- args[2] %>% as.numeric()
gbd_round_id  <- args[3] %>% as.numeric()

if (gbd_round_id == 6){
  r_id <- 6
} else if (gbd_round_id == 7){
  r_id <- 9
}

rr_max_file <- "FILEPATH"
out_dir <- paste0("FILEPATH")

### Other variables
locations   <- get_location_metadata(location_set_id = 39, release_id = r_id)[level >= 3,location_id]
loc_id      <- locations[task_id]

message("starting SEV calc for location_id ", loc_id)


#--PULL RRMAX ------------------------------------------------------------------

message("reading RRmax draws from ", rr_max_file)
rr_max <- fread(rr_max_file) %>%
  .[, draw := as.numeric(gsub("draw_", "", draw))]


#--PULL PAFs AND MERGE ---------------------------------------------------------
message("reading PAF draws")

paf_files <- list.files(paste0(out_dir, "/paf/"))
paf_files <-  paf_files[grep(paste0("^",loc_id,"_"), paf_files)] 

paf <- rbindlist(lapply(paste0(out_dir, "/paf/", paf_files), fread))

dt <- merge(paf, rr_max, by=c("age_group_id", "sex_id", "cause_id", "draw"))
dt <- unique(dt)


#--CALC SEV --------------------------------------------------------------------
message("calculating SEV")
dt[paf < 0, paf := 0]
dt[, sev := (paf / (1 - paf)) / (rr - 1)]
dt[rr <= 1, sev := 0]

if (nrow(dt[is.infinite(sev) | round(paf, 12) == 1, ]) > 0) stop("SEVs are infinite, PAFs of 100%?")
if (nrow(dt[is.na(sev), ]) > 0) warning("SEVs are NA")

dt[sev >= 1, sev := 1]
dt[, mean_sev := mean(sev, na.rm=T), by=c("location_id", "year_id", "age_group_id", "sex_id", "cause_id")]
dt[is.na(sev), sev := mean_sev][, mean_sev := NULL]

if (nrow(dt[is.na(sev), ]) > 0) stop("SEVs are still NA")


#--AVERAGE ACROSS CAUSE AND SAVE ---------------------------------------------------------------

# save risk_cause-specific sev if requested
save_draws(dt, n_draws, paste0(out_dir, "/risk_cause/"))

vars <- c("location_id", "year_id", "age_group_id", "sex_id", "draw", "cohort_start_age")
dt <- dt[, .(sev = mean(sev)), by=vars] 

save_draws(dt, n_draws, paste0(out_dir, "/risk/"))

