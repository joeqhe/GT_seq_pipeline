#!/usr/bin/python

""""
Identifies the densest SNP locations
"""

##########################################################################################
# Generate bins containing SNPs
##########################################################################################


function_generate_bins (SNP_locations, window_size)
		
		
##########################################################################################
# Scan vesca genome for dense SNP windows
##########################################################################################

function_read_SNP_location (vesca_SNP_file)

for results in range (1, window_number + 1):
	print "calculating %d of %d densest %d nucleotides in %s" % (results, window_number, window_size, file)
	function_full_scan (SNP_locations, bins, window_size)

