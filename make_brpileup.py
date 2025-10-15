#!/usr/bin/env python3

"""
Build an artificial SAM file of consensus mappings from split sequences.

This script takes split sequences, their original alignment regions, and coordinate information
to generate consensus alignments. It handles quality filtering and can process alignments
in parallel using multiple threads.

IMPORTANT: Script assumes PHRED + 33 quality encoding. Will break with other encodings.
"""

import argparse
import skbio.alignment as skaln
from multiprocessing import Pool
import re

def argparser():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--consensus', help="Set prefix of output 'sam' file for consensus alignments.")
    parser.add_argument('-t', '--threads', type=int, help='Set number of threads to use. Default 1', default=1)
    parser.add_argument('bed', help='Path to a bed file containing coordinate information matching the primary fasta file.')
    parser.add_argument('regions', help='Path to a fasta file containing reference alignment regions from the original round of alignment')
    parser.add_argument('split', help='Path to file containing the split sequences (.fa)')
    args = parser.parse_args()
    return args

def rand_aln_p(bl: int, sl: int) -> float:
    """
    Calculate probability of random perfect alignment between sequences.
    
    Args:
        bl: Length of base sequence
        sl: Length of sequence to align
        
    Returns:
        float: Probability of random perfect alignment
    """
    if sl < bl:
        return 1-(1-.25**sl)**(bl-sl)
    else:        
        return 1-(1-.25**bl)**(sl-bl)

def phred_code_qual(sym):
    """
    Convert between PHRED+33 quality scores and symbols.
    
    Args:
        sym: Quality score symbol (str) or value (int)
        
    Returns:
        Converted quality value or symbol
    """
    if isinstance(sym, str):
        return ord(sym) - 33
    elif isinstance(sym, int):
        return chr(sym+33)
    else:
        print("Error: quality score conversion not understood. Check input")
        raise ValueError

def bedread(bed_path: str) -> dict:
    """
    Read BED file and organize entries by read name.
    
    Args:
        bed_path: Path to BED file
        
    Returns:
        dict: Mapping of read names to lists of region tuples (chrom, start, end)
    """
    bd = {}
    with open(bed_path) as bedf:
        for entry in bedf:
            spent = entry.strip().split()
            read = spent[3]
            reg = tuple(spent[0:3])
            bd[read] = bd.get(read, [])
            if reg not in bd[read]:
                bd[read].append(reg)
    return bd

def regionread(region_path: str) -> dict:
    """
    Read region sequences from FASTA file.
    
    Args:
        region_path: Path to regions FASTA file
        
    Returns:
        dict: Mapping of sequence names to sequences
    """
    rd = {}
    cur = None
    with open(region_path) as regf:
        for entry in regf:
            if entry[0] == '>':
                cur = entry.strip()[1:].partition('::')[0]  # forwards compatibility
            else:
                rd[cur] = entry.strip()
    return rd

def splitread(split_path: str) -> dict:
    """
    Read split sequences and their quality scores.
    
    Args:
        split_path: Path to split sequences file
        
    Returns:
        dict: Mapping of read names to lists of [sequence, quality] pairs
    """
    sd = {}
    cur = None
    with open(split_path) as splf:
        for entry in splf:
            if entry[0] == '+':
                continue  # Skip quality score identifier lines
            if entry[0] == '@':
                cur = entry.partition('_')[0][1:]
                sd[cur] = sd.get(cur, [])  # Initialize entry
                tmp = []
            elif len(tmp) == 0:  # Add sequence
                tmp.append(entry.strip())
            elif len(tmp) == 1:  # Add quality scores
                tmp.append(entry.strip())
                sd[cur].append(tmp)
    return sd

def generate_consensus(region: str, seqs: list) -> tuple:
    """
    Generate consensus sequence from multiple aligned sequences.
    
    Args:
        region: Reference region sequence
        seqs: List of [sequence, quality] pairs
        
    Returns:
        tuple: ([consensus sequence, quality string], flag)
    """
    region_counts, flag = basecount(region, seqs)
    consensus = []
    quality = []  # Actually just count, makes stratifying easier
    for i, bcs in sorted(region_counts.items(), key=lambda x:x[0]):
        b, q = best_base(bcs)
        consensus.append(b)
        quality.append(str(q))
    return [''.join(consensus), ''.join(quality)], flag

def best_base(basecounts: dict) -> tuple:
    """
    Determine consensus base from base counts.
    
    Args:
        basecounts: Dictionary of base counts
        
    Returns:
        tuple: (consensus base, count)
    """
    if all([c==0 for c in basecounts.values()]):
        return ('N',0)
    nonz = sorted([(b,c) for b,c in basecounts.items() if c > 0], key=lambda x:x[1], reverse=True)
    if len(nonz) == 1:
        return nonz[0]
    elif len(nonz) > 1:
        return ('N', nonz[0][1])  # Use N if there's any controversy
    else:
        return ('N',0)

def basecount(region: str, seqs: list) -> tuple:
    """
    Count bases at each position across all sequences.
    
    Args:
        region: Reference region sequence
        seqs: List of [sequence, quality] pairs
        
    Returns:
        tuple: (position->base->count dict, error flag)
    """
    # Initialize counts for each position and base
    counts = {i:{b:0 for b in 'ACGT'} for i in range(len(region)+1)}
    flag = False  # Track alignment issues
    
    for seq, qual in seqs:
        # Filter very short sequences
        if rand_aln_p(len(region), len(seq)) <= .0001:
            # Align sequence to region
            seqaln = skaln.StripedSmithWaterman(
                seq, 
                gap_open_penalty=5, 
                gap_extend_penalty=2, 
                match_score=1, 
                mismatch_score=-3, 
                zero_index=False
            )
            aln = seqaln(region)
            
            # Skip sequences with indels
            if aln['cigar'].count('I') == 0 and aln['cigar'].count('D') == 0:
                matched = 0
                for sec in aln['cigar'].split("M"):
                    reg = re.split('[A-Z]',sec)
                    if len(reg[-1]) > 0:
                        matched += int(reg[-1])
            else:
                flag = True
                continue
                
            index = aln['target_begin']
            if index != -1:  # Valid mapping
                try:
                    rseq = aln.target_sequence[aln.target_begin:aln.target_end_optimal]
                except:
                    flag = True
                    continue
                    
                qseq = aln.query_sequence[aln.query_begin:aln.query_end]
                qqual = qual[aln.query_begin:aln.query_end]
                
                # Filter low quality alignments
                if aln.optimal_alignment_score / (len(rseq)+1) < .5:
                    flag = True
                    continue
                    
                # Count high quality bases
                for i, rb in enumerate(rseq):
                    assert len(qseq) == len(rseq)
                    try:
                        qb = qseq[i]
                        if qb in 'ACGT' and phred_code_qual(qqual[i]) > 12:  # Q13+ = <5% error rate
                            counts[index + i + 1][qb] += 1
                    except:
                        print("Error trying to count alignment to region", region, 
                              aln.aligned_target_sequence, aln.query_sequence, aln['cigar'])
                        
    return counts, flag

def sam_entry(seqqual: list, coord: tuple, name: str, flag: bool=False) -> str:
    """
    Generate SAM format entry.
    
    Args:
        seqqual: [sequence, quality] pair
        coord: Coordinate tuple (chrom, start, end)
        name: Read name
        flag: Error flag
        
    Returns:
        str: SAM format entry
    """
    coord = coord[0]
    fv = '512' if flag else '0'  # Flag low quality alignments
    return f"{name}\t{fv}\t{coord[0]}\t{coord[1]}\t60\t{len(seqqual[0])}M\t*\t0\t0\t{seqqual[0]}\t{seqqual[1]}"

def mapper(input_iter: tuple) -> str:
    """
    Process single read alignment for parallel execution.
    
    Args:
        input_iter: Tuple of (name, coordinates, region, split sequences)
        
    Returns:
        str: SAM format entry
    """
    k, coord, reg, splts = input_iter
    cons, flag = generate_consensus(reg, splts) #version of the reference with either Ns for ambiguous/no information or the 
    #cons = [reg, 'I' * len(reg)]  # Use reference sequence with placeholder quality
    #flag = False  # Assume no errors
    samstr = sam_entry(cons, coord, k, flag)
    return samstr

def main():
    """Main function to process alignments and generate consensus SAM file."""
    args = argparser()
    
    # Read input files into dictionaries
    bed_d = bedread(args.bed)  # Read coordinates
    reg_d = regionread(args.regions)  # Read reference sequences
    spl_d = splitread(args.split)  # Read split sequences
    
    # Process alignments in parallel
    with open(args.consensus, 'w+') as outf:
        with Pool(args.threads) as p:
            arguments = [(k, bed_d[k], reg_d[k], spl_d[k]) 
                        for k in spl_d.keys() 
                        if k in reg_d and k in bed_d]
            samstrs = p.imap_unordered(mapper, arguments)
            for s in samstrs:
                if s != '':
                    print(s, file=outf)

if __name__ == "__main__":
    main()
