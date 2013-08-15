import sys
import argparse

parser = argparse.ArgumentParser(description='Phase pairs of het-sites.')

parser.add_argument('samtools', nargs=1,
                    help='Output of samtools phase (use - for stdin).')

parser.add_argument('--min-read-count', 
                    default=10, action='store', type=int,
                    help='Minimum number of reads covering each phase of a '
                    'pair of het sites.')

args = parser.parse_args()

# Threshold constants
MIN_COUNT_THRESHOLD = args.min_read_count

if args.samtools[0] == '-':
    infile = sys.stdin
else:
    infile = open(args.samtools[0])


def split_input_to_records(infile):
    '''Split an input stream into records, yielding the lines within
    each record.'''
    cur_record = []
    for line in infile:
        if line.startswith('CC'):
            continue
        if line.startswith('//'):
            yield cur_record
            cur_record = []
        else:
            cur_record.append(line)

def parse_record(record):
    '''Parse the lines of a record and count the number of supporting
    reads for each pair of sites in the record.'''
    chrom = None # is assigned each M line, but they should match
    het_map = dict()
    read_counts = dict()
    for line in record:
        if line.startswith('M'):
            _,chrom,_,pos,_,_,hetidx,_,_,_,_ = line.split()
            het_map[int(hetidx)] = int(pos)
        
        if line.startswith('EV'):
            _,_,evchrom,hetidx,_,_,_,_,_,seq = line.split()[:10]
            assert evchrom == chrom
            hetidx = int(hetidx)
            pos1 = het_map[hetidx]
            all1 = cigar[0]
            for hetidx2 in xrange(1,len(seq)):
                all2 = cigar[hetidx2]
                if all2 == 'N': continue
                
                pos2 = het_map[hetidx + hetidx2]
                try:
                    read_counts[(pos1,pos2,all1,all2)] += 1
                except:
                    read_counts[(pos1,pos2,all1,all2)] = 1

    return chrom, read_counts

def split_counts_in_pairs(counts):
    '''Change the counts table, so we have counts for each pair in separate
    tables.  This makes it easier to output the phased pairs.'''
    pos_counts = dict()
    for (pos1,pos2,all1,all2),count in counts.items():
        tbl = pos_counts.setdefault( (pos1,pos2), {} )
        tbl[(all1,all2)] = count
    return pos_counts


for record in split_input_to_records(infile):
    chrom,counts = parse_record(record)
    pos_counts = split_counts_in_pairs(counts)
    pos_in_record = pos_counts.keys()
    pos_in_record.sort()
    
    for pos1,pos2 in pos_in_record:
        allele_counts = [(k,v) for (k,v) in pos_counts[(pos1,pos2)].items()
                                   if v >= MIN_COUNT_THRESHOLD ]
                                   
        if len(allele_counts) != 2:
            continue # skip if we can't count exactly two
        
        (a11,a12),count1 = allele_counts[0]
        (a21,a22),count2 = allele_counts[1]
        
        if a11 == a21 or a12 == a22:
            continue # one of the sites is not het (in the supported reads)
            
        print "%s\t%d\t%d\t%s%s\t%d\t%s%s\t%d" % (
                chrom, pos1, pos2,
                a11,a12, count1,
                a21,a22, count2
            )

