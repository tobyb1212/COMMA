from argparse import ArgumentParser
import gzip
import re
from collections import Counter

# Code is adapted from the following script https://github.com/yachielab/Interstellar/blob/main/utils/sciATAC_header_parse.py 

def get_option():
    argparser = ArgumentParser()
    argparser.add_argument('-i', '--input', type=str, required=True, help="Input FASTQ file")
    argparser.add_argument('-b', '--barcode_list', type=str, required=True, help="File containing unique barcodes (one per line)")
    argparser.add_argument('-q', '--baseQuality', type=int, default=38, help="Base quality of generated FASTQ")
    argparser.add_argument('-o', '--outname', type=str, default="barista", help="Output file name")
    argparser.add_argument('-d', '--outdir', type=str, default=".", help="Output directory")
    return argparser.parse_args()

def fastq_generator(input_file):
    """Generator to read FASTQ file line by line."""
    with gzip.open(input_file, mode="rt") as seqf:
        for c, line in enumerate(seqf):
            if c % 4 == 0:
                yield line.strip()

def load_unique_barcodes(barcode_file):
    """Load unique barcodes from a file into an iterator."""
    with open(barcode_file, "r") as f:
        barcodes = [line.strip() for line in f.readlines()]
    return iter(barcodes)

def count_top_barcodes(fastq_iter, top_n=737000):
    """Count the occurrences of unique barcodes from FASTQ headers and return the top N barcodes."""
    barcode_counter = Counter()
    
    for headline in fastq_iter:
        headline = headline.split(" ")
        original_barcode = re.sub(r":.+$", "", headline[1]) # Extract Barcode from .fastq Header
        barcode_counter[original_barcode] += 1 
    
    top_barcodes = dict(barcode_counter.most_common(top_n)) # Extract the top_n barcodes from the .fastq file to correspond with the number of unique 10X barcodes
    
    print(f"Total unique barcodes: {len(barcode_counter)}")
    return top_barcodes  # Return the list of barcodes and counts

def header_processor(fastq_iter, unique_barcodes, allow_barcodes, quality, outdir, outname):
    """Process FASTQ headers, assign new unique barcodes consistently, and write to output."""
    barcode_map = {}  # Used to store original barcode and map to new barcode
    I1 = []
    
    for c, headline in enumerate(fastq_iter, start=1):  # Start from 1 to prevent first-read skipping
        headline = headline.split(" ")
        header = headline[0]
        original_barcode = re.sub(r":.+$", "", headline[1])  # Extract original barcode from .fastq header
        
        if original_barcode in allow_barcodes:# Assign a unique barcode only if this is in the list of allowed barcodes (top_n barcodes)
            if original_barcode not in barcode_map:
                try:
                    barcode_map[original_barcode] = next(unique_barcodes) # Assign each original barcode a new 10X compatible barcode
                except StopIteration:
                    raise ValueError("Not enough unique barcodes in the list!")
        else:
            barcode_map[original_barcode] = 'A' * 16 # If the barcode is not in the top_n barcodes, assign the barcode as AAAAAAAAAAAAAAAA, not CellRanger compatible

        new_barcode = barcode_map[original_barcode]  # Ensure that the same original barcode receives the same new barcode
        
        # Construct new I2 .fastq file
        I1.append(header)
        I1.append(new_barcode)     # New barcode 
        I1.append(f"+{original_barcode}")             # Original barcode
        I1.append(chr(quality) * len(new_barcode))  # Quality string

        # Periodically write to file (every 4M reads)
        if c%4000000==0:
            print(c,"reads have been processed.")
            I1="\n".join(I1)+"\n"
            if c==0:
                with gzip.open(outdir+"/"+outname+"_S1_L001_I2_001.fastq.gz",mode="wt") as wI1:
                    wI1.write(I1)
            else:
                with gzip.open(outdir+"/"+outname+"_S1_L001_I2_001.fastq.gz",mode="at") as wI1:
                    wI1.write(I1)
            I1=[]

    print(c,"reads have been processed.")
    I1="\n".join(I1)+"\n"
    with gzip.open(outdir+"/"+outname+"_S1_L001_I2_001.fastq.gz",mode="at") as wI1:
        wI1.write(I1)


if __name__ == "__main__":
    opt = get_option()
    input_fastq = opt.input # Parse input fastq file
    quality = opt.baseQuality + 30

    unique_barcodes = load_unique_barcodes(opt.barcode_list)

    fastq_iter = fastq_generator(input_fastq)
    allow_barcodes = count_top_barcodes(fastq_iter)
    fastq_iter = fastq_generator(input_fastq)
    header_processor(fastq_iter, unique_barcodes, allow_barcodes, quality, opt.outdir, opt.outname) # Construct new I2 file.
