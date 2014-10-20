from collections import namedtuple
import gzip

SNP = namedtuple('SNP', 'rsid chrom position genotype')

yes = lambda x: True

class HumanGenotype(object):

    _snps = None

    def __init__(self, datafile):

        if datafile.endswith('.gz'):
            opener = gzip.open
        else:
            opener = open

        with opener(datafile) as fp:
            _snps_list = []
            for line in fp:
                # -- skip header
                if line.startswith('#'):
                    continue
                rsid, chrom, pos, gtype = line.strip().split()
                pos = int(pos)
                _snps_list.append(SNP(rsid, chrom, pos, gtype))
            self._snps = tuple(_snps_list) # think of genome as immutable

    def summarize(self):
        return "This dataset contains [%d] SNPs." % len(self._snps)

    def select(self, rsid=yes, chrom=yes, position=yes, genotype=yes):
        return tuple(snp for snp in self._snps
                     if all((rsid(snp.rsid),
                             chrom(snp.chrom),
                             position(snp.position),
                             genotype(snp.genotype))))

if __name__ == "__main__":
    import sys

    genome = HumanGenotype(sys.argv[1])
    print genome.summarize()
    print genome.select(rsid=(lambda s: s == 'rs4477212'))
    #print genome.select()
