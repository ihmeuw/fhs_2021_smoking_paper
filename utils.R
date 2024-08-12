###############################################################################################################
## Purpose: Functions for smoking PAF & SEV pipeline
###############################################################################################################

##### Note: This file has not yet undergone code review.


### To create custom prevalence

age_cohort <- function(a) {
  if (a %in% c(7:19, 30, 31)) {
    return(a+1)
  }
  if (a %in% c(20)) {
    return(30)
  }
  if (a %in% c(32)) {
    return(235)
  }
  if (a %in% c(235)) {
    return(235)
  }
}

add_mort_rows <- function(df, a, type, y = 9999){
  
  if (type == "current"){
    need_to_add <- nrow(df[age_group_id == a]) == 0
  } else {
    need_to_add <- nrow(df[age_group_id == a & year_id == y]) == 0
  }
  
  if (need_to_add){
    
    df_sub <- expand_grid(location_id = l, age_group_id = a, year_id = y, sex_id = s, variable = c(0:999)) %>% 
      as.data.table %>%
      .[,variable := paste0("draw_", variable)] 
    
    if (type == "current"){
      df_sub[,current_mort := 1]
      df_sub$year_id <- NULL
    } else if (type == "former"){
      df_sub[,former_mort := 1]
    } else if (type == "new_former"){
      df_sub[,new_former_mort := 1]
    }
    
    df <- rbind(df, df_sub)
  }
  
  return(df)
}


### Custom PAFs

new_cess <- function(cess_orig, t){
  
  new_cess_dist <- function(v) cess_orig(v - t)
  
  return(new_cess_dist)
  
}

calc_integral_new_former <- function(unif_start = 0, unif_max, rr, minval, maxval) {
  integral<-integrate(function(x) (dunif(x, min = unif_start, max = unif_max) * rr(x)), minval, maxval, subdivisions = 300, stop.on.error = FALSE, rel.tol = 0.0000001)$value
  return(integral)
} 

format_weighted_rr <- function(integral){
  df <- as.data.frame(integral) %>% as.data.table() %>%
    .[, draw := 0:(length(integral)-1)] %>%
    .[, year_id := future_year] %>%
    .[, cohort := cohort_a] %>%
    .[, age_group_id := current_age] %>%
    setnames(old = "integral", new = "value")
  
  return(df)
  
}

compute_exposure_weighted_rr <- function(distribution_age, rr_age){
  print(paste0("Creating current smoker integral for those that quit at ",distribution_age, " and are now ", rr_age))
  
  rr_current_list<-split(copy(rr_current[age_group_id==rr_age]), by = "draw")
  rr_current_list<-lapply(rr_current_list, function(x) approxfun(x$exposure, x$rr, method = "linear", rule =2))
  exp_weighted_rr <- mapply(calc_integral, exp = get(paste0(exp_def, "_", distribution_age, "_", s, "_", years)), 
                            rr = rr_current_list, 0, cap_exposure(exp_def)) 
  
  return(exp_weighted_rr)
}

compute_scaled_former_rr <- function(scaling_factor_name, rr_age){

  new_max_temp<-data.table(new_max=get(scaling_factor_name), draw=0:999) 
  
  ## Pull former smoker RR curve + merge
  rr_former_temp <- copy(rr_former) %>%
    .[age_group_id==rr_age] %>% 
    .[, max:=max(rr), by = c("draw", "age_group_id", "sex_id")] %>%
    merge(new_max_temp, by = "draw")
  
  if (c==544) {
    rr_former_temp$rr_new<-mapply(transform_former, orig_val = rr_former_temp$rr, new_max = 1, new_min = rr_former_temp$new_max)
  } else {
    rr_former_temp$rr_new<-mapply(transform_former, orig_val = rr_former_temp$rr, new_max = rr_former_temp$new_max, orig_max = rr_former_temp$max)
  }
  
  rr_former_list<-split(rr_former_temp, by = "draw")
  rr_former_list<-lapply(rr_former_list, function(x) approxfun(x$exposure, x$rr_new, method = "linear", rule =2))
  
  return(rr_former_list)
  
}

expand_causes_parent <- function(out){
  # Expand causes
  if (297 %in% unique(out$cause_id)) { out<-expand_causes(data_in = out, cause_in = 297, causes_out = c(934, 946, 947, 954), drop_orig = T) }
  if (417 %in% unique(out$cause_id)) {out<-expand_causes(data_in = out, cause_in = 417, causes_out = c(418, 419, 420, 421, 996), drop_orig = T)}
  if (487 %in% unique(out$cause_id)) {out<-expand_causes(data_in = out, cause_in = 487, causes_out = c(845, 846, 847, 848, 943), drop_orig = T)}
  if (494 %in% unique(out$cause_id)) {out<-expand_causes(data_in = out, cause_in = 494, causes_out = c(495, 496, 497), drop_orig = T)}
  if (587 %in% unique(out$cause_id)) {out<-expand_causes(data_in = out, cause_in = 587, causes_out = c(976), drop_orig = T)}
  
  return(out)
}


### Differential Mortality

read_relative_risks <- function(file, path){
  
  print(file)
  df <- fread(paste0(path, file)) %>%
    .[, filename := gsub(".csv", "", file)] %>%
    separate(filename, into = c("location_id", "sex_id", "cause_id"))
  
  df <- copy(df) %>%
    .[age_group_id != 235, keep := cohort] %>%
    .[age_group_id == 235, keep := min(cohort), by = c("draw", "year_id")]
  table(df$age_group_id, df$year_id)
  df <- df[keep == cohort]
  df[,c("keep")] <- NULL
  
  ## Append other draw rows if only one per l/a/s/y/c
  if (nrow(df) < 140000){
    if (length(unique(df$draw)) == 1){
      if (unique(df$draw) == 0){
        temp <- copy(df)
        for (d in c(1:999)){
          temp[,draw := d]
          df <- rbind(df, temp)
        }
      }
    }
  }

  if (nrow(df) < 140000){
    have <- copy(df[draw == 0, c("age_group_id", "year_id")])
    need <- copy(df[draw == 1, c("age_group_id", "year_id")]) %>%
      .[,need := 1]
    
    have <- merge(have, need, by = c("age_group_id", "year_id"), all.x = T) %>%
      .[is.na(need)]
    if (nrow(have) > 0){
      for (i in nrow(have)){
        a <- have[i]$age_group_id
        y <- have[i]$year_id
        
        temp <- copy(df[year_id == y & age_group_id == a])
        if (nrow(temp) == 1){
          for (d in c(1:999)){
            temp[,draw := d]
            df <- rbind(temp, df)
          }
        }
      }
    }
  }
  
  
  df <- expand_causes_parent(df)
  
  return(df)
}

compute_diff_mortality <- function(path, type, l, s){
  files <- list.files(path)
  files <- files[files %like% paste0("^",l, "_",s, "_")]
  
  if (type %in% c("former", "former_new")){
    files <- files[!(files %like% paste0("^",l, "_",s, "_878"))]
    files <- files[!(files %like% paste0("^",l, "_",s, "_923"))]
  }
  
  start <- Sys.time()
  df <- rbindlist(lapply(files, read_relative_risks, path = path)) %>%
    .[location_id == l & sex_id == s] %>%
    .[,location_id := as.numeric(location_id)] %>%
    .[,sex_id := as.numeric(sex_id)]%>%
    .[,cause_id := as.numeric(cause_id)]
  end <- Sys.time()
  end - start
  
  if (type %in% c("former")){
    df <- df[!(cause_id %in% c(878, 923))]
    df <- df[cohort >= 9] 
  }
  if (type %in% c("former_new")){
    df <- df[!(cause_id %in% c(878, 923))]
    df <- df[cohort >= 9] 
    df <- df[year_id != 2022] 
  }
  
  ## Append non-smoker
  df_nonsmoke <- copy(df[cause_id == 322]) %>% 
    .[, cause_id := 999] %>%
    .[, value := 1]
  
  df <- rbind(df, df_nonsmoke) %>%
    merge(rates, by = c("location_id", "age_group_id", "sex_id", "cause_id", "draw"), all.x = T) %>% #not year
    .[is.na(death_rate), death_rate := 0] %>% #NAs are only measured for older ages
    .[,log_rr := log(value)] %>%
    .[,total_death_rate := sum(death_rate), by = c(demo_cols, "draw")] %>%
    .[,weighted_death_rate := (death_rate * log_rr)/total_death_rate] %>%
    .[,avg_rr := sum(weighted_death_rate), by = c(demo_cols, "draw")]
  
  df[,c("cause_id", "value", "death_rate", "log_rr", "total_death_rate", "weighted_death_rate", "cohort")] <- NULL
  df <- unique(df)
  
  if (type == "current"){
    df[, year_id := NULL]
    df <- unique(df)
  }
  
  for (a in c(9:20, 30:32, 235)){
    
    if (type == "current"){
      
      df_sub <- copy(df[age_group_id == a])
      if (nrow(df_sub) == 0){
        df_sub <- expand_grid(location_id = l,
                              age_group_id = a,
                              sex_id = s, 
                              draw = c(0:999),
                              avg_rr = 0)
        df <- rbind(df, df_sub)
      }
      
    } else {
      for (y in seq(start_yr, (start_yr + 5*7), 5)){
        df_sub <- copy(df[age_group_id == a & year_id == y])
        if (nrow(df_sub) == 0){
          df_sub <- expand_grid(location_id = l,
                                age_group_id = a,
                                sex_id = s, 
                                year_id = y, 
                                draw = c(0:999),
                                avg_rr = 0)
          df <- rbind(df, df_sub)
        }
      }
    }
    
  }
  
  df <- df %>%
    .[,draw := paste0("draw_", draw)] %>%
    setnames(old = "draw", new = "variable")
  
  return(df) 
}


### PAF post processing

interp_yrs <- function(year){
  print(year)
  
  if (year == 2021){
    prior_have_yr <- 2020
  } else if (year %in% c(2023:2025)){
    prior_have_yr <- 2022
  } else{
    optional_yrs <- year - c(1:5)
    prior_have_yr <- have_yrs[have_yrs %in% optional_yrs] %>% last()
  }
  
  sev_y <- copy(paf) %>%
    .[, paf := approx(year_id, paf, xout = year)$y, by = c("location_id", "sex_id", "draw", "cause_id", "age_group_id")] %>%
    .[year_id == prior_have_yr] %>%
    .[,year_id := year] 
  
  return(sev_y)
}


### SEV

save_draws <- function(dt, n_draws, out_dir) {
  
  id_cols <- c("location_id", "year_id", "age_group_id", "sex_id", "measure_id", "metric_id", "cohort_start_age")
  
  if ("cause_id" %in% names(dt)) {
    id_cols <- c("rei_id", "cause_id", id_cols)
  } else {
    id_cols <- c("rei_id", id_cols)
  }
  
  dt <- dt %>%
    .[, `:=` (rei_id=99, measure_id=29, metric_id=3)]
  dcast(as.formula(paste(paste(id_cols, collapse = "+"), "~ draw")), value.var = "sev") %>%
    setnames(paste0(0:(n_draws - 1)), paste0("draw_", 0:(n_draws - 1)))
  
  setorderv(dt, cols = id_cols)
  
  message("saving to ", out_dir, file)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  file <- paste0(loc_id, ".csv")
  write.csv(dt, paste0(out_dir, file), row.names = F)
}


