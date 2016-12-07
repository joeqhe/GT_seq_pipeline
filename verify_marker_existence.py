#!/usr/bin/python

"""

"""


##########################################################################################
# Read map files
##########################################################################################


for map_file in octoploid_mapping_populations:
	os.chdir (working_directory)
	os.chdir (map_file)
	os.chdir ("uniqfinal")
	for LG in os.listdir(os.getcwd()):
		if LG.endswith (".map"):
			dictionary_name = map_file [0:5] + "_" + LG [:-4]
			function_read_map_file (LG)


##########################################################################################
# Search for markers in octoploid maps
##########################################################################################


for linkage_group in all_map_position:
	for probeset_ID in all_map_position [linkage_group]:
		if probeset_ID not in all_probeset_ID:
			all_probeset_ID.append (probeset_ID [-8:])


all_probeset_ID.append ("88858788") # testing


for window in densest_SNP_windows:
	for SNPs_present in densest_SNP_windows [window] [0]:
		if str (SNPs_present [-8:]) in all_probeset_ID:
			print SNPs_present
			
			
			
			
