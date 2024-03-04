#!/usr/bin/env python3
import argparse
import sys
import urllib.request

# Update gc.prt
url = 'ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt'
urllib.request.urlretrieve(url, '../config/gc.prt')

prs = argparse.ArgumentParser(description='Generate binary header output file gc.h for SANS by using the given genetic code file from NCBI')
prs.add_argument('-i', '--input', required=False, help='Input: genetic code file from NCBI. Default: ../config/gc.prt')
prs.add_argument('-o', '--output', required=False, help='Output: resulting hex dump gc.h for best performance. Default: ../src/gc.h')
args = prs.parse_args()

# Check for possible args
input_file = '../config/gc.prt'
output_file='../src/gc.h'

if args and args.input: input_file = args.input
if args and args.output: output_file = args.output

# read input file
with open(input_file, 'r') as input:
    data = input.read()

# write hexdump
out = []
out.append('unsigned char config_gc_prt[] = {')
list = [data[i:i+12] for i in range(0,len(data), 12)]
for i, x in enumerate(list):
    line = ', '.join(['0x{val:02x}'.format(val=ord(c)) for c in x])
    out.append('  {hexval}{comma}'.format(hexval=line, comma=',' if i < len(list)-1 else ''))

out.append('};')
out.append('unsigned int config_gc_prt_len = {len};'.format(len=len(data)))

# write output
with open(output_file, 'w') as output:
    output.write('\n'.join(out))