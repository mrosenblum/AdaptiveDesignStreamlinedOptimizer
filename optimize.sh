# Created by Michael Rosenblum 9/26/17 based on template from Josh Betz
# Script for running adaptive design optimizer pipeline 
# computing 1 stage and 2 stage optimized designs
# username <- from web (current username)
# rand <- from web (random folder)
#
# project_root: root directory of project
#	e.g. /users/jbetz/jhbc/loe/rosenblum/optimizer/
# data_dir: data subdirectory off project root
#	e.g. /(project.root)/data/
# code_dir: data containing optimizer R code
#	e.g. /(project.root)/code/2016-12/
# optimizer_file: .R file containing optimizer code
#
# 1. Create a directory to hold the results of a simulation:
# 2. Run optimizer
#

username=$1
random=$2
project_root=/users/agherman/rosenblum/optimizer/${username}/${random}
code_dir=${project_root}/code/
data_dir=${project_root}/data/
log_dir=${project_root}/log

run_date=$(date +%Y-%m-%d)
#run_name=${JOB_NAME}_${JOB_ID}
#run_dir=${data_dir}${run_date}/${run_name}
run_dir=${data_dir}

#setup_rda=optimizer_parameters_${run_date}.rda

pipeline_r=pipeline_streamlined.R
optimizer_file=design_optimizer.R
performance_file=ComputePerformanceMetrics.R
binary_search_file=Utility_BinarySearch.R
backend_1tvc_file=Backend2Populations1Arm.R
backend_2tvc_file=Backend2Populations2Arms.R
report_generator_file=optimizer_report_final.Rnw

function simulation_setup {
	if [ -d "$project_root" ]; then
		################################################################
		# 1. Create Directory Structure
		################################################################
		# If data directory does not exist, create it
		if [ ! -d "$run_dir/results" ]; then
			echo "Creating results directory ${run_dir}"
			#mkdir -p ${log_dir}
			mkdir -p ${run_dir}/results
		fi
		
		################################################################
		# 2. Run Pipeline
		################################################################
		Rscript ${code_dir}${pipeline_r} \
			"code.dir=\"$code_dir\"" \
			"data.dir=\"$run_dir\"" \
			"optimizer.file=\"$optimizer_file\"" \
			"performance.file=\"$performance_file\"" \
			"binsearch.file=\"$binary_search_file\"" \
			"backend.1tvc.file=\"$backend_1tvc_file\"" \
			"backend.2tvc.file=\"$backend_2tvc_file\"" \
			"report.generator.file=\"$report_generator_file\"" \
			"job.name=\"$JOB_NAME\"" \
			"job.id=$JOB_ID" \
			"run.date=\"$run_date\""\
			"initial.seed=12345"
		wait
		cd $run_dir/results
		latex trial_design_performance_report.tex
		pdflatex trial_design_performance_report
	else
		echo 'ERROR: project_root ($project_root) does not exist.'
	fi
}

simulation_setup
