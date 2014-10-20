import sys
import lib23
from collections import namedtuple

RiskModel = namedtuple("RiskModel", "rsid effects")
risk_models = (
    # -- values pulled from 23andMe's technical report
    RiskModel("rs2187668", dict(CC=0.48, CT=3.35)),
    RiskModel("rs6822844", dict(GT=0.79, GG=1.14)),
    RiskModel("rs6441961", dict(TT=1.43, CC=0.82, CT=1.09)),
    RiskModel("rs9851967", dict(CT=0.99, TT=0.81, CC=1.21)),
)
celiac_baseline_risk = 0.0012

def odds(p):
    return p / (1.0 - p)

def inv_odds(p):
    return p / (1.0 + p)

def main(args):

    infile = args[0]
    genome = lib23.HumanGenotype(infile)

    risk_markers = genome.select(rsid=lambda m: m in set(e.rsid for e in risk_models))

    print "=========================="
    print "Celiac Disease Risk Report"
    print "=========================="
    print

    print "Individual Marker Reports"
    print "--------------------------"
    print
    marker_odds = []
    for i in range(len(risk_markers)):
        this_risk_marker = risk_markers[i]
        this_risk_model = [e for e in risk_models if e.rsid == this_risk_marker.rsid][0]
        print "[Risk Marker %2d] rsid: %10s  genotype: %2s  odds ratio: %0.2f" % \
            (i + 1,
             this_risk_marker.rsid,
             this_risk_marker.genotype,
             this_risk_model.effects[this_risk_marker.genotype])
        marker_odds.append(this_risk_model.effects[this_risk_marker.genotype])
    print

    combined_odds = 1.0
    for modds in marker_odds:
        combined_odds *= modds

    print "Overall Risk"
    print "------------"
    print " Average Celiac Risk: [%0.2f%%]" % (100 * celiac_baseline_risk)
    print " Your Celiac Risk:    [%0.2f%%]" % \
        (100 * inv_odds(combined_odds * odds(celiac_baseline_risk)))


if __name__ == '__main__':
    main(sys.argv[1:])
