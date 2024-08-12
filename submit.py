import click
import sys
from pathlib import Path
from fhs_lib_orchestration_interface.lib.submit_cluster_job import submit_cluster_job_request

def submit_job(version, reference, reduction, reduction_type, target_yr, last_year, ref_version, draws, gbd_round_id):
  
  code_dir = str(Path(__file__).parents[1].absolute() / "lib" / "")
  
  prev_script = "create_custom_prevalence"
  paf_alt_script = "compute_pafs_custom"

  print(gbd_round_id)
  
  if gbd_round_id == '6':
      num_tasks_l = 346
  elif gbd_round_id == '7':
      num_tasks_l = 482
      
  num_tasks_l_s = num_tasks_l*2
  num_tasks_l_s_c = int(num_tasks_l_s * 35.5)
    
  if reference:
    
    step1_args = [version, gbd_round_id]
    last_jid = launch_job(
        args_list = step1_args,
        script = "compute_pafs_reference",
        memory = 50,
        cores = 10,
        j_hold = 1,
        num_tasks = num_tasks_l_s_c
      )
    
  else:

    ### Step 1: Prevalence without differential mortality
    step1_args = [version, 'FALSE', target_yr, reduction, gbd_round_id, reduction_type] #diff mortality = False
    step1_jid = launch_job(
        args_list = step1_args,
        script = prev_script,
        memory = 10,
        cores = 4,
        j_hold = 1,
        num_tasks = num_tasks_l_s
      )
    
    ### Step 2: Compute exposure-weighted RRs
    step2_args = [version, 'FALSE', target_yr, reduction, gbd_round_id] #diff mortality = False
    step2_jid = launch_job(
        args_list = step2_args,
        script = paf_alt_script,
        memory = 50,
        cores = 10,
        j_hold = step1_jid,
        num_tasks = num_tasks_l_s_c
      )

    ### Step 3: Compute mortality rates by smoking status
    step3_args = [version, gbd_round_id]
    step3_jid = launch_job(
        args_list = step3_args,
        script = "compute_differential_mortality",
        memory = 40,
        cores = 20,
        j_hold = step2_jid,
        num_tasks = num_tasks_l_s
      )

    ### Step 4: Prevalence with differential mortality
    step4_args = [version, 'TRUE', target_yr, reduction, gbd_round_id, reduction_type] #diff mortality = True
    step4_jid = launch_job(
        args_list = step4_args,
        script = prev_script,
        memory = 10,
        cores = 4,
        j_hold = step3_jid,
        num_tasks = num_tasks_l_s
      )

    ### Step 5: Compute PAFs
    step5_args = [version, 'TRUE', target_yr, reduction, gbd_round_id] #diff mortality = True
    last_jid = launch_job(
        args_list = step5_args,
        script = paf_alt_script,
        memory = 50,
        cores = 10,
        j_hold = step4_jid,
        num_tasks = num_tasks_l_s_c
      )
  
  
  ### Step 6: Clean up PAFs 
  step6_args = [code_dir, version, draws, gbd_round_id, last_year, str(reference).upper(), ref_version]
  step6_jid = launch_job(
      args_list = step6_args, 
      script = "paf_postprocessing", 
      memory = 20, 
      cores = 4, 
      j_hold = last_jid, 
      num_tasks = num_tasks_l_s
    )
    
  ### Step 7: Compute SEVs
  step7_args = [version, draws, gbd_round_id]
  step7_jid = launch_job(
      args_list = step7_args, 
      script = "compute_sevs", 
      memory = 20, 
      cores = 4, 
      j_hold = step6_jid, 
      num_tasks = num_tasks_l
    )
  
  
def launch_job(args_list, script, memory, cores, j_hold, num_tasks):

  hours = "02:00:00"
  que = "all.q"  
  project = "proj_forecasting"

  shell = "FILEPATH"
  image = "FILEPATH"
  
  code_dir = str(Path(__file__).parents[1].absolute() / "lib" / "")
  full_script = code_dir + "/" + script + ".R"
  
  arguments = " ".join(args_list)
  
  command = (
      f"{shell} "
      f"-i {image} "
      f"-s {full_script} "
      f"{arguments}"
    )

  jid = submit_cluster_job_request(
      job_name = script,
      memory = memory,
      cores = cores,
      command_to_run = command,
      runtime = hours,
      cluster_queue = que,
      cluster_project = project,
      run_on_archive_node = True,
      job_holds = [j_hold],
      number_of_array_tasks = num_tasks
    )
  
  return jid

