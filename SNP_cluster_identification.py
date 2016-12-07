#!/usr/bin/python

""""
Identifies the densest SNP locations
"""

##########################################################################################
# Generate bins containing SNPs
##########################################################################################


function_read_SNP_location (vesca_SNP_file)

SNP_locations_original = copy.deepcopy (SNP_locations)

function_generate_bins (SNP_locations, window_size)


##########################################################################################
# Scan vesca genome for dense SNP windows
##########################################################################################


for results in range (1, window_number + 1):
	print "calculating %d of %d densest %d nucleotides in %s" % (results, window_number, window_size, vesca_SNP_file )
	
	function_full_scan (SNP_locations, bins, window_size)

for location in densest_SNP_windows:
	densest_window = [None, 0]
# densest_window  = [window_start, SNP_count]
	
	if densest_SNP_windows [location] [-1] > densest_window [-1]:
		densest_window = densest_SNP_windows [location]
		
	
print "The densest SNP window contains %d SNPs with these markers %s" % (densest_window [-1],densest_window [0])


##########################################################################################
# Scan densest_SNP_windows for overlaping SNPs
##########################################################################################


function_remove_overlapping_windows (densest_SNP_windows, window_size)
	
print "removed %d of densest_SNP_windows result(s) due to overlap" % (remove_counter)
	
	
	
