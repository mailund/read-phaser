import sys
import argparse

parser = argparse.ArgumentParser(description='''
Phase haplotypes from phased pairs.
''')

parser.add_argument('pairs', nargs=1,
                    help='List of phased pairs (use - for stdin).')

parser.add_argument('--buffer', 
                    default=1000, action='store', type=int,
                    help='''
    Number of pairs to read in before processing a batch of connected
    components. The default should be a good choice for must purposes.
    ''')


args = parser.parse_args()

# FIXME: handle missing files more gracefully
if args.pairs[0] == '-':
    infile = sys.stdin
else:
    infile = open(args.pairs[0])



class node(object):
    def __init__(self, chrom, pos):
        self.chrom = chrom
        self.pos = pos
        self.in_edges = []
        self.out_edges = []
        self.component = None

# FIXME: Doesn't handle switches from one chromosome to another!!!
def collect_buffer_of_nodes(buffer, max_src = 0):
    '''Collect a number of pairs into a buffer to be processed.
    
    Parameter *buffer* is the left-overs from the last buffer processed.
    This is the buffer that will be extended.
    
    Parameter *max_src* is the largest source node seen in the previous
    call.
    
    Returns the new buffer, the largest source node seen, and a status
    flag indicating if there are more pairs in the input file.
    '''
    number_inserted = 0
    for line in infile:
        chrom, pos1, pos2, phase1, _, phase2, _ = line.split()
        pos1, pos2 = int(pos1), int(pos2)

        try:
            src = buffer[(chrom,pos1)]
        except:
            src = node(chrom,pos1)
            buffer[(chrom,pos1)] = src
        
        try:
            dst = buffer[(chrom,pos2)]
        except:
            dst = node(chrom,pos2)
            buffer[(chrom,pos2)] = dst
        
        src.out_edges.append( (phase1,phase2,dst) )
        dst.in_edges.append(  (phase1,phase2,src) )

        if src.pos > max_src:
            max_src = src.pos

        number_inserted += 1
        if number_inserted >= args.buffer:
            return buffer, max_src, True
    return buffer, max_src, False

def split_in_components(nodes):
    '''Split a buffer of nodes into connected components.'''

    def assign_components(node, component):
        def dfs(n):
            if n.component is not None:
                assert n.component == component
            else:
                n.component = component
                for _,_,src in n.in_edges:
                    dfs(src)
                for _,_,dst in n.out_edges:
                    dfs(dst)
        dfs(node)

    loci = nodes.keys()
    loci.sort()

    component = 0
    for locus in loci:
        node = nodes[locus]
        if node.component is None:
            assign_components(node,component)
            component += 1

    components = [ list() for i in xrange(component) ]
    for locus in loci:
        node = nodes[locus]
        components[node.component].append(node)

    return components

def split_components(components, max_src):
    '''Split a list of components in those that are done
    and those that are potentially still incomplete (based on *max_src*).
    
    Returns the finished components as a list and the non-finished
    as nodes in a dictionary matching the buffer format.
    '''
    finished_components = []
    new_buffer = {}
    for component in components:
        max_pos = max(n.pos for n in component)
        if max_pos < max_src:
            finished_components.append(component)
        else:
            for n in component:
                n.component = None # don't save the assignment for next time
                new_buffer[(n.chrom,n.pos)] = n
    return finished_components, new_buffer

## Phase a connected component
class InconsistentComponent:
    pass # FIXME: give warning message here rather than to stderr
def phase_component(graph):
    '''Phase a finished component and write the phase to stdout.'''

    for idx,node in enumerate(graph):
        out_allele_1 = [phase1[0] for phase1,phase2,n in node.out_edges]
        out_allele_2 = [phase2[0] for phase1,phase2,n in node.out_edges]
        in_allele_1 = [phase1[1] for phase1,phase2,n in node.in_edges]
        in_allele_2 = [phase2[1] for phase1,phase2,n in node.in_edges]
        alleles = set(out_allele_1+out_allele_2+in_allele_1+in_allele_2)
        if len(alleles) != 2:
            print >> sys.stderr, "Non biallelic", alleles
            raise InconsistentComponent
        node.alleles = tuple(alleles)
        node.phased = False

    def dfs(node):
        assert node.phased
        for phase1,phase2,n in node.out_edges:
            if node.alleles[0] == phase1[0]:
                n_phase = phase1[1],phase2[1]
            else:
                n_phase = phase2[1],phase1[1]
            
            if n.phased:
                if n.alleles != n_phase:
                    print >> sys.stderr, "Inconsistent phasing:",
                    print >> sys.stderr, n.alleles, "!=",
                    print >> sys.stderr, n_phase
                    print >> sys.stderr, graph[0].chrom,
                    print >> sys.stderr, [x.pos for x in graph]
                    raise InconsistentComponent
            else:
                n.alleles = n_phase
                n.phased = True
                dfs(n)
            
        for phase1,phase2,n in node.in_edges:
            if node.alleles[0] == phase1[1]:
                n_phase = phase1[0],phase2[0]
            else:
                n_phase = phase2[0],phase1[0]
            
            if n.phased:
                if n.alleles != n_phase:
                    print >> sys.stderr, "Inconsistent phasing:",
                    print >> sys.stderr, n.alleles, "!=",
                    print >> sys.stderr, n_phase
                    print >> sys.stderr, graph[0].chrom,
                    print >> sys.stderr, [x.pos for x in graph]
                    raise InconsistentComponent
            else:
                n.alleles = n_phase
                n.phased = True
                dfs(n)        
    
    first = graph[0]
    first.phased = True # arbitrary phase
    dfs(first)

    last = graph[-1]
    indices = [n.pos for n in graph]
    hap1 = [n.alleles[0] for n in graph]
    hap2 = [n.alleles[1] for n in graph]
    
    print "%s\t%d\t%d\t%s\t%s\t%s" % \
        (first.chrom, first.pos, last.pos, 
         ','.join(map(str,indices)),
         ''.join(hap1), ''.join(hap2))


def phase_finished_components(buffer, max_src, flush = False):
    '''Build components, phase those that are completely read in,
    and return an updated buffer with the nodes that are not completed.'''

    components = split_in_components(buffer)
    if flush:
        finished, new_buffer = components, dict()
    else:
        finished, new_buffer = split_components(components, max_src)

    for component in finished:
        try:
            phase_component(component)
        except InconsistentComponent:
            pass
    
    return new_buffer



## MAIN LOOP, READING IN PAIRS AND PROCESSING BUFFERS
buffer, max_src, more_pairs = collect_buffer_of_nodes(dict())
new_buffer = phase_finished_components(buffer, max_src, not more_pairs)
while more_pairs:
    buffer, max_src, more_pairs = collect_buffer_of_nodes(new_buffer,max_src)
    new_buffer = phase_finished_components(buffer, max_src, not more_pairs)

#print 'digraph reads {'
#for locus in loci:
#    node = nodes[locus]
#    for phase1,phase2,dst in node.out_edges:
#        print '"%s:%d<%d>"' % (node.chrom,node.pos,node.component), '->',
#        print '"%s:%d<%d>"' % (dst.chrom,dst.pos,dst.component),
#        print '[label="%s|%s"]' % (phase1,phase2), ';'
#print '}'