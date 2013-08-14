import sys
import argparse

parser = argparse.ArgumentParser(description='''
Phase haplotypes from phased pairs.
''')

parser.add_argument('pairs', nargs=1,
                    help='List of phased pairs (use - for stdin).')

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
nodes = dict()

## COLLECT ALL PAIRS AND BUILD SMALL GRAPHS
for line in infile:
    chrom, pos1, pos2, phase1, _, phase2, _ = line.split()
    pos1, pos2 = int(pos1), int(pos2)

    try:
        src = nodes[(chrom,pos1)]
    except:    
        src = node(chrom,pos1)
        nodes[(chrom,pos1)] = src
        
    try:
        dst = nodes[(chrom,pos2)]
    except:
        dst = node(chrom,pos2)
        nodes[(chrom,pos2)] = dst
        
    src.out_edges.append( (phase1,phase2,dst) )
    dst.in_edges.append(  (phase1,phase2,src) )

loci = nodes.keys()
loci.sort()

## COLLECT GRAPHS IN CONNECTED COMPONENTS - FIXME: handle 
## component assignments on the fly so we don't have to contain all the
## data in memory before we can process it...
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


## Phase a connected component
class InconsistentComponent:
    pass # FIXME: give warning message here rather than to stderr
def phase_component(graph):

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


for graph in components:
    try:
        phase_component(graph)
    except InconsistentComponent:
        pass



#print 'digraph reads {'
#for locus in loci:
#    node = nodes[locus]
#    for phase1,phase2,dst in node.out_edges:
#        print '"%s:%d<%d>"' % (node.chrom,node.pos,node.component), '->',
#        print '"%s:%d<%d>"' % (dst.chrom,dst.pos,dst.component),
#        print '[label="%s|%s"]' % (phase1,phase2), ';'
#print '}'