# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 09:54:28 2017
This program takes NGS sequencing read data as input csv file with 
start and length columns.
The program takes loci file for specific locations on the sequence
and calculates the coverage for that location (number of times that location
is contained in the read fragments). The output is a csv file 
containing the coverage for the requested loci. 

@author: brnissen
"""

import pandas as pd


# Files used: Assume files located in same working directory as python program.
# Included option to save output file as different file from loci input.
input_file = 'reads.csv'
loci_file = 'loci.csv'
output_file = 'loci.csv'

# Input the csv data for the NGS file to a Pandas Dataframe. 
input_df = pd.read_csv(input_file)


# Preprocessing 
# Sort the reads by start and length
sorted_df = input_df.sort_values(['start','length']) 
# Reset the index so we can see what is going on easier
sorted_df = sorted_df.reset_index(drop = True)

# Add a column to our dataframe, calculating end position of each sequence
# Knowing the end position makes it easier to determine 
# if the loci are contained in the sequences. 
# The (-1) because the lenght includes the start and end location.
sorted_df['end'] = sorted_df.start + sorted_df.length-1


# Create the dataframe for the loci positions. 
# Coverage column shown, but contains no data.
loci_df = pd.read_csv(loci_file)

# Find the maximum and minimum values of the loci_df, 
# so we can make the reads to be searched smaller
maximum,minimum = max(loci_df.position),min(loci_df.position)

# Throw out reads that are outside of the range of loci
temp_df = sorted_df.mask(sorted_df.start>maximum)
temp_df2 = temp_df.mask(temp_df.end<minimum)
new_sorted = temp_df2.dropna().astype(int)


# For each position in the loci file, check if the location is contained 
# in each read in the read file.
# This is done by comparing each read's start and end position, 
# if it is between or equal to either, it is contained.
# This produces a dataframe query, and counts the reads that meet the 
# query criteria. The count of each loci position is appended to 
# a coverage column, to be exported to csv later.
# This is the slowest part of the program because it loops through each row.
def get_count(locus):  # Get the count for each locus
    less_than = new_sorted.start<=locus 
    greater_than = new_sorted.end>=locus
    
    count = new_sorted.start[less_than & greater_than].count()
    return count

# Build the coverage column of the dataframe
loci_df.position = loci_df.position.apply(get_count) 

# Export the loci_df to csv file
loci_df.to_csv(output_file,index=False) 
# index=False does not include index column in the csv output
