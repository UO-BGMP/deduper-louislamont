#!/usr/bin/env python3

import argparse
import re
import itertools

def get_arguments():
    #Give description of deduper, take UMI input, SAM input, and output files
    parser = argparse.ArgumentParser(description="This program takes a sorted \
.sam file and outputs a .sam file with PCR duplicates removed. A set of known \
UMIs can be supplied, otherwise all 8-mers wil be used.", add_help=True)
    parser.add_argument("-f", "--file", help="Takes a .sam file (must be sorted \
by SAMtools)", required=True, type=str)
    parser.add_argument("-u", "--umi", help="Takes a .txt file containing \
UMI sequences", required=False, type=argparse.FileType('r'))
    parser.add_argument("-p", "--paired", help="Designates whether the .sam file \
is paired end", required=False, action="store_true")
    return parser.parse_args()

# Correct position for soft clipping
def positionchange(pos, cigar):
    # If soft clipped bases are found at beginning of cigar string
    if re.match("^[0-9]+S", cigar):
        # Extract the number and subtract from position, otherwise return position
        soft=cigar[0:cigar.find("S")]
        return (int(pos)-int(soft))
    else:
        return int(pos)

# Get argparse arguments
args=get_arguments()

# Check if paired flag is active and quit if it is
if args.paired:
    print("Fatal error. No paired-end functionality, sucka. Exiting.")
    quit()

# Create empty set for UMIs
umilist = set()
# Check is UMI flag is active and read from file if it is
if args.umi:
    with args.umi as umi:
        for line in umi:
            umilist.add(line.strip())
# Otherwise, create a set of all possible 8-mers
else:
    print("No UMI file supplied. Using all possible 8-mer combinations.")
    bases=["A","C","T","G"]
    umilist=set(''.join(i) for i in itertools.product(bases, repeat=8))

# Define the reference UMI/Position as null
umiref=""
posref=0
chrref= ""

# Make dictionary to hold unique read info
seen_before={}

# Read file directory and create output name
infile=args.file
outname=infile+"_deduped"

### Main loop
# Open files for reading and writing
with open(outname, "w") as out:
    with open(infile, "r") as sam:
        for line in sam:
            stripped=line.strip()
            # If line is a header line (starts with @), output to file
            if stripped.startswith("@"):
                print(stripped, file=out)
            # else if it's not, check for PCR duplicates
            else:
                # Split line and read UMI (follows last : in first field)
                splitline=stripped.split("\t")
                lineumi=splitline[0][splitline[0].rfind(":")+1:]
                # Check if UMI exists in dict, if not, read next line
                if lineumi not in umilist:
                    continue
        
                # Extract position, cigar string from the line
                linepos=int(splitline[3])
                linechr=str(splitline[2])
                linecigar=str(splitline[5])
                
                # Adjust position with soft clipping from cigar string
                adjpos=positionchange(linepos, linecigar)

                # Check if line's position, chromosome, umi are in the dict
                # Increment value if it has been seen before
                line_id=linechr+"-"+str(linepos)+"-"+lineumi
                if line_id in seen_before:
                    seen_before[line_id]+=1
                    continue
                # Create entry in dictionary otherwise
                else:
                    seen_before[line_id]=seen_before.get(line_id, 1)
                    print(stripped, file=out)
