#!/usr/bin/env python3

def bitcheck(bit):
    '''Takes a bitwise flag and checks it for strandedness. \
Assumes read is mapped (otherwise returns None) and data are single-stranded. \
Returns "+" or "-" depending on strand.'''
    # Check if read is mapped, if not, return None
    if (bit &4) == 4:
        return None
    # First define strand as + strand
    strand = "+"
    # Check if bit flag 16 is true, if so, read is on - strand, so change
    if ((bit & 16) == 16):
        strand = "-"
    return strand
