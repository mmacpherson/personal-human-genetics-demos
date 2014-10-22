import sys
import math
import lib23
from collections import namedtuple

AIM = namedtuple("AIM", "rsid ref var NW_ref_freq SE_ref_freq AJ_ref_freq")

def parse_aims(infile):

    aims = {}
    for line in open(infile):
        if line.startswith("SNP"):
            continue
        rsid, chrom, pos, ref, var, NW, SE, AJ = line.strip().split()
        aims[rsid] = AIM(rsid, ref, var, float(NW), float(SE), float(AJ))

    return aims

def prob_allele(allele, source, aim):
    if allele not in (aim.ref, aim.var):
        print >> sys.stderr, "[Warning] Found allele [%s] in marker [%s], and expected one of [%s, %s]." % (allele, aim.rsid, aim.ref, aim.var)
    if allele == aim.ref:
        return getattr(aim, "%s_ref_freq" % source)
    return 1.0 - getattr(aim, "%s_ref_freq" % source)


def main(args):

    genomefile, aimsfile = args
    genome = lib23.HumanGenotype(genomefile)
    aims = parse_aims(aimsfile)

    aims_markers = genome.select(rsid=lambda m: m in aims.keys())

    parentages = (('NW', 'NW'),
                  ('NW', 'SE'), ('SE', 'SE'),
                  ('NW', 'AJ'), ('SE', 'AJ'), ('AJ', 'AJ'),
    )
    lnls = [0.0 for p in parentages]

    for marker in aims_markers:
        a1, a2 = list(marker.genotype)
        aim = aims[marker.rsid]
        for i, (pa, pb) in enumerate(parentages):
            like = 0.5 * prob_allele(a1, pa, aim) * prob_allele(a2, pb, aim) + \
                   0.5 * prob_allele(a1, pb, aim) * prob_allele(a2, pa, aim)
            lnls[i] += math.log(like)

    likes = [math.exp(lnl) for lnl in lnls]
    likes_normed = [like / sum(likes) for like in likes]
    for p, lnl, lkn in zip(parentages, lnls, likes_normed):
        print p, "%9.3f" % lnl, "%7.2f%%" % (100 * lkn)

if __name__ == '__main__':
    main(sys.argv[1:])
