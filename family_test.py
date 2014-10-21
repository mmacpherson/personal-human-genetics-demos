import sys
import math
import lib23
from collections import namedtuple


def main(args):

    genome1file, genome2file = args
    genome1 = lib23.HumanGenotype(genome1file)
    genome2 = lib23.HumanGenotype(genome2file)

    subset1 = genome1.select(position=lambda p: p % 19999 == 0)
    subset2 = genome1.select(position=lambda p: p % 19999 == 0)





if __name__ == '__main__':
    main(sys.argv[1:])
