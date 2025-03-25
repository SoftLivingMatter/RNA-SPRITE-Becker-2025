#!/bin/bash

base_dir='/path/to/slurm/logs'

start_dir=$PWD
cd $base_dir

# sub_dir='100_seq_retries'
# cd $sub_dir

python $start_dir/get_usage.py \
	<(reportseff -s CD --format jobid,elapsed,maxrss,reqcpus) \
	>report.txt

cd ..

cd $start_dir
