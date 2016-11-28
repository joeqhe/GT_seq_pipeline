#!/usr/bin/python

"""

"""


##########################################################################################
# 2. Read map files
##########################################################################################

for map_file in octoploid_mapping_populations:
	os.chdir (working_directory)
	os.chdir (map_file)
	os.chdir ("uniqfinal")
	for LG in os.listdir(os.getcwd()):
		if LG.endswith (".map"):
			dictionary_name = map_file [0:5] + "_" + LG [:-4]
			function_read_map_file (LG)