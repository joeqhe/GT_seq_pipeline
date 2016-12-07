#!/usr/bin/python

"""

# Functions for use in GT-seq primer design pipeline.

# Includes:
	1. Read SNP locations; function_read_SNP_location (vesca_SNP_file):
	2. Read map file; function_read_map_file (map_file):
	3. Generate bins; function_generate_bins (SNP_locations, window_size):
	4. Full scan for dense windows; function_full_scan (SNP_locations, bins, window_size):
	5. Remove ovelapping windows; def function_remove_overlapping_windows (densest_SNP_windows, window_size)):

"""

##########################################################################################
# 1. Read SNP locations
##########################################################################################


def function_read_SNP_location (vesca_SNP_file):
	"""
	# Reads CSV (vesca_SNP_file) containing [SNP name, linkage group, physical location]
	
	# Only works when max(physical location) <= 10**8 - 1
	
	# A dictionary, SNP_locations, is outputted {SNP_id:[(LG)12345678],...} where LG is 
	the linkage group and the numbers represent the SNP location. Also outputs the total
	number of SNPs in file, the number of SNPs per linkage group and max SNP position
	
	# Creates a copy of SNP_locations as SNP_locations_original so SNP_locations can be
	manipulated later in pipeling without losing data
	"""
	
	os.chdir (working_directory)
	with open(vesca_SNP_file) as csv_file:
		csv_data = csv.reader(csv_file, delimiter = ",")
		
		max_SNP_position = 0
		for row_count,row in enumerate (csv_data):
			SNP_id, LG, position = row
			LG = int (LG)
			SNP_locations [SNP_id] = LG * 100000000 + int(position)
#			SNP_locations key in form (LG)12345678 where the numbers is the SNP position
			if LG not in LG_count:
				LG_count [LG] = 1
			else:
				LG_count [LG] +=  1
			if (SNP_locations [SNP_id] - LG * 100000000) > max_SNP_position:
				max_SNP_position = (SNP_locations [SNP_id] - LG * 100000000)
	
		print "this file contains %d SNPs" % (row_count)
		
		for LG in LG_count:
			print "with %d SNPs in LG %d" % (LG_count [LG], LG)
			
		print "the highest positioned SNP is at %d" % (max_SNP_position)
		if max_SNP_position > (10 ** 9 - 1):
			print "ERROR: linkage group too large!"


##########################################################################################
# 2. Read map file
##########################################################################################


def function_read_map_file (map_file):
	"""
	# Reads a .map file and outputs a dictionary named (filepath wrt working_directory)
	in format {SNP_name: [maternal_pos, paternal_pos, combined_pos], ...} with NA where
	data not available. Also returns which file is read and number of SNPs
	
	# Outputs data to all_map_position {map_name_1: {SNP_name_1: [maternal_pos, paternal_pos, combined_pos], ...}, map_name_2: {...}}
	
	# Ignores first line of file (assumed to be header)
	"""
	
	with open (map_file) as file:
		map_file_data = {}
		
		for row_count, row in enumerate (file):
		
			if row_count > 0:
				
				while "  " in row:
					row = row.replace ("  "," ")
				row = row.split()
				
				for element in row:
					if element.isdigit():
						row [element] = int(element)
				map_file_data [row [0]] = row [1:]		
		print "%s has %d SNPs" % (dictionary_name, len (map_file_data))	
		all_map_position [dictionary_name] = map_file_data


##########################################################################################
# 3. Generate bins
##########################################################################################


def function_generate_bins (SNP_locations, window_size):

	"""
	Inputs:
		1. SNP_locations {SNP_id: "(LG)12345678} where LG is linkage group and 12345678 is SNP
		location
		2. window_size
	
	Generates bins {bin_id: {SNP_id:}} where bin_id is integer division of SNP_id by
	window_size. Populates bin_id with SNP_ids that are present in that bin

	"""
	
	
	for SNP_id in SNP_locations:
		bin_id  = SNP_locations [SNP_id] / window_size
		if bin_id not in bins:
			bins [bin_id] = set([SNP_id])
		else:
			bins [bin_id].add(SNP_id)


##########################################################################################
# 4. Full scan for dense windows
##########################################################################################


def function_full_scan (SNP_locations, bins, window_size):
	"""
	# Takes inputs:
		1. SNP_locations: {SNP_id:[(LG)12345678],...} where LG is the linkage group and
		the numbers represent the SNP location
		2. bins: {bin_id: [SNP_id,...]} where bin_id is (LG)12345678 / window_size
		3. window_size: integer of window size
		
	# Outputs densest_SNP_windows {window_start: [[SNP_id,...], SNP_count]}
	# Also deletes SNP_locations [SNP_id] of SNP_id in best_window
	"""
	
	
	best_window = [None, None, 0]
# best_window = [window_start, [SNP_id,...] , SNP_count]
	
	for bin_id in bins:
		SNP_set = bins [bin_id]
# SNP_set = [SNP_id,...] where SNP_id is in bin_id
		
		if bin_id + 1 in bins:
			SNP_set = SNP_set | bins [bin_id +1]
# SNP_set includes bin_id + 1 if it exists, to catch dense SNP regions overlapping bins
		
		for SNP_id in bins [bin_id]:
			window_start = SNP_locations [SNP_id]
			current_window = [window_start, [], 0]
# Define window parameters for each SNP_id in SNP_set
		
			for SNP_id2 in SNP_set:
				if window_start <= SNP_locations [SNP_id2] <= window_start + window_size:
					current_window [1].append (SNP_id2)
					current_window [2] += 1
# Search every SNP in SNP_set for closeness to window_start; add to current_window if close
		
				if best_window [2] < current_window [2]:
					best_window = current_window
# Update best_window with current_window if current_window has more SNPs
	
	densest_SNP_windows [best_window [0]] = [best_window [1], best_window [2]]
# Append densest_SNP_windows with best_window
	
	for SNP_id in best_window [1]:
		while SNP_id in bins [SNP_locations [SNP_id] / window_size]:
			bins [SNP_locations [SNP_id] / window_size] .remove (SNP_id)
		del SNP_locations [SNP_id]
# Deletes SNP_id from bins and SNP_locations to allow another iteration of this function


##########################################################################################
# 5. Remove overlapping windows
##########################################################################################


def function_remove_overlapping_windows (densest_SNP_windows, window_size):
	"""
	# Inputs:
		1. densest_SNP_windows: {window_start: [SNP_id,...], SNP_count}
		2. window_size: integer window size
	
	# Outputs:
		1. densest_SNP_windows: {window_start: [SNP_id,...], SNP_count} with overlapping
		windows set to [None, None, 0], keeping densest window
	
	"""
	
	
	for SNP_window in densest_SNP_windows:
		overlap_windows = {}
# overlap_windows = {SNP_window: SNP_count}
		
		for SNP_window2 in densest_SNP_windows:
			if SNP_window2 - window_size <= SNP_window <= SNP_window2 + window_size:
				overlap_windows [SNP_window2] = densest_SNP_windows [SNP_window2] [-1]
# Check every SNP_window2 for overlap with SNP_window
		
		if len (overlap_windows) > 1:
			densest_SNP_count = max (overlap_windows.values ())
			
			for overlap_window in overlap_windows:
				if  overlap_windows [overlap_window] == densest_SNP_count:
					keep_window = overlap_window
				
				if overlap_window != keep_window:
					densest_SNP_windows [overlap_window] = [None, None, 0]
					remove_counter += 1

	
##########################################################################################
# 6. Convert probeset ID to SNP ID
##########################################################################################


def function_probesetID_to_SNPID (probeset_ID):
	"""
	# Inputs:
		1. Probeset_ID: probeset_ID name
		
	# Outputs:
		1. Changes 
	
	
	
	"""

print "placeholder function"















