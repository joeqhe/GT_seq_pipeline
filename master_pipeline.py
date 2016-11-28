#!/usr/bin/python

"""

Master script to generate "optimal" primers for GT-seq in octoploid strawberry

1. identify densest SNP windows in vesca genome
2. filter to remove overlapping windows
3. reference SNPs to ensure they exist in >= 1 octoploid mapping population
4. identify 

Run with: 

working_directory = "/home/hejoe/GT_seq_primer_design"
os.chdir (working_directory)
os.chdir ("full_pipeline")
execfile ("master_pipeline.py")

"""

import csv
import copy
import os


##########################################################################################
#Input files
##########################################################################################


working_directory = "/home/hejoe/GT_seq_primer_design"
os.chdir (working_directory) # start in correct directory
vesca_SNP_file = "vesca_snp_positions.csv"

octoploid_mapping_populations = [
"EMxFE_map_2016-04-14",
"FLxCH_map_2016-04-14",
"RGxHA_map_2016-04-04"]


##########################################################################################
#Input parameters
##########################################################################################


window_size = 1000 # Size of window to search through in vesca
window_number = 250 # Number of windows to output


##########################################################################################
#Empty datasets
##########################################################################################


SNP_locations = {} # function_read_SNP_location: SNP location output dictionary 
LG_count = {} # function_read_SNP_location: counts number of LG in input file
bins = {} # function_generate_bins (SNP_locations, window_size): generates bins

densest_SNP_window = {}
# function_full_scan (SNP_locations, bins, window_size): output of results


all_map_position = {} # function_read_map_file: SNP location output dictionary

for map_file in octoploid_mapping_populations:
	os.chdir (working_directory)
	os.chdir (map_file)
	os.chdir ("uniqfinal")
	for LG in os.listdir(os.getcwd()):
		if LG.endswith (".map"):
			dictionary_name = map_file [0:5] + "_" + LG [:-4]
			all_map_position [dictionary_name] = {}
			
os.chdir (working_directory)
# function_read_map_file: populate all_map_position with dictionary for each LG


##########################################################################################
# Load functions
##########################################################################################


os.chdir (working_directory)
os.chdir ("full_pipeline")
execfile ("GT_seq_functions.py")


##########################################################################################
# Execute pipeline
##########################################################################################


os.chdir (working_directory)
os.chdir ("full_pipeline")
execfile ("SNP_cluster_identification.py")

#execfile ("verify_marker_existence.py")






