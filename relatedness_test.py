# -*- coding: utf-8 -*-
import sys
import math
import lib23
from collections import namedtuple

NC = 80 # -- parameter controlling number of characters to use for
         #    graphical display of relatedness

UNRELATED = 0
HALF_IBD  = 1
FULL_IBD  = 2

pchar = {UNRELATED:'.',
         HALF_IBD:'░',
         FULL_IBD:'█',
}

def is_half_ibd(ss1, ss2, thresh):
    """Are two unphased diploid genomic regions, ss1 and ss2, compatible
       with sharing at least one haplotype?
    """
    num_opposite_homozygotes = 0
    is_homo = lambda gt: gt[0] == gt[1]
    for snp1, snp2 in zip(ss1, ss2):
        s1h = is_homo(snp1.genotype) and snp1.genotype[0] in 'ACTG'
        s2h = is_homo(snp2.genotype) and snp2.genotype[0] in 'ACTG'
        if not (s1h and s2h):
            continue
        num_opposite_homozygotes += 0 if snp1.genotype == snp2.genotype else 1

    return num_opposite_homozygotes <= thresh * len(ss1)

def is_full_ibd(ss1, ss2, thresh):
    """Do two unphased diploid genomic regions, ss1 and ss2, share both
       haplotypes?
    """
    num_not_identical = 0
    for snp1, snp2 in zip(ss1, ss2):
        s1g = snp1.genotype
        s2g = snp2.genotype
        if s1g[0] not in 'ACTG':
            continue
        if s2g[0] not in 'ACTG':
            continue
        num_not_identical += 0 if s1g == s2g else 1

    return num_not_identical <= thresh * len(ss1)


def main(args):

    genome1file, genome2file = args

    # -- load them up
    genome1 = lib23.HumanGenotype(genome1file)
    genome2 = lib23.HumanGenotype(genome2file)

    # -- pull out just chromosome 1
    subset1 = genome1.select(chrom=lambda c: c == "1")
    subset2 = genome2.select(chrom=lambda c: c == "1")

    # -- what markers do both of these genomes share?
    rsid1 = set(snp.rsid for snp in subset1)
    rsid2 = set(snp.rsid for snp in subset2)
    shared_rsids = rsid1 & rsid2

    # -- sort markers by position, just the intersection
    subset1 = sorted([snp for snp in subset1 if snp.rsid in shared_rsids],
                     key=lambda m: m.position)
    subset2 = sorted([snp for snp in subset2 if snp.rsid in shared_rsids],
                     key=lambda m: m.position)
    assert len(subset1) == len(subset2), "Subsets not of same length!"

    # -- figure out how many SNPs are in a window, then get boundaries
    window_size = len(subset1) / NC
    thresh = 1.0 / window_size
    window_boundaries = [(window_size * i, window_size * (i + 1)) for i in range(NC)]
    windows = [UNRELATED] * NC

    # -- check each (nonoverlapping) window for relatedness
    for (win_n, window) in enumerate(window_boundaries):
        ss1 = subset1[window[0]:window[1]]
        ss2 = subset2[window[0]:window[1]]
        windows[win_n] = HALF_IBD if is_half_ibd(ss1, ss2, thresh) else UNRELATED
        if windows[win_n] == HALF_IBD:
            windows[win_n] = FULL_IBD if is_full_ibd(ss1, ss2, thresh) else HALF_IBD

    num_unrelated = sum(c == UNRELATED for c in windows)
    num_half_ibd = sum(c == HALF_IBD for c in windows)
    num_full_ibd = sum(c == FULL_IBD for c in windows)

    print "Relatedness Summary on Chromosome 1"
    print "==================================="
    print
    print "Unrelated: %3d of %3d windows (%6.2f%%)" % (num_unrelated, NC,
                                                     100.0 * num_unrelated / NC)
    print "Half-IBD : %3d of %3d windows (%6.2f%%)" % (num_half_ibd, NC,
                                                     100.0 * num_half_ibd / NC)
    print "Full-IBD : %3d of %3d windows (%6.2f%%)" % (num_full_ibd, NC,
                                                     100.0 * num_full_ibd / NC)
    print
    print "Graphical Summary Along Chromosome 1"
    print "------------------------------------"
    print
    print "  (Legend: Unrelated: %s  Half-IBD: %s  Full-IBD: %s)" % (pchar[UNRELATED],
                                                                     pchar[HALF_IBD],
                                                                     pchar[FULL_IBD])
    print
    print "".join(pchar[c] for c in windows)



if __name__ == '__main__':
    main(sys.argv[1:])
