import argparse

versionStr = "You are using a Pandas based FitHiC written by Omega!"

def parse_args(args):
    parser = argparse.ArgumentParser(description="A pandas based FitHiC runner", add_help=False)
    req_args = parser.add_argument_group('Required Arguments')
    req_args.add_argument("-i", "--interactions", dest="intersfile",\
                      help="REQUIRED: interactions between fragment pairs are \
                      read from INTERSFILE", required=True)

    req_args.add_argument("-f", "--fragments", dest="fragsfile", \
                      help="REQUIRED: midpoints (or start indices) \
                      of the fragments are read from FRAGSFILE",\
                      required=True)

    req_args.add_argument("-o", "--outdir", dest="outdir", \
                      help="REQUIRED: where the output files\
                      will be written", required=True)

    req_args.add_argument("-r", "--resolution", dest="resolution", type=int,
                      help="REQUIRED: If the files are fixed size, please \
                      supply the resolution of the dataset here; otherwise, \
                      please use a value of 0 if the data is not fixed size." \
                      , required=True)
    
    misc_args = parser.add_argument_group('Miscellaneous Arguments')
    misc_args.add_argument("-L", "--lowerbound", dest="distLowThres", type=int,
                      help="RECOMMENDED: lower bound on the intra-chromosomal \
                      distance range (unit: base pairs). DEFAULT no limit. \
                      Suggested limit is 2x the resolution of the input files",
                      required=False, default=-1)
    misc_args.add_argument("-U", "--upperbound", dest="distUpThres", type=int,
                      help="RECOMMENDED: upper bound on the intra-chromosomal \
                      distance range (unit: base pairs). DEFAULT no limit. \
                      STRONGLY suggested to have a limit for large genomes,\
                      such as human/mouse. ex. '1000000, 5000000, etc.'",
                      required=False, default=-1)
    misc_args.add_argument("-b", "--noOfBins", dest="noOfBins", type=int, \
                      help="RECOMMENDED: number of equal-occupancy (count) \
                      bins. Default is 100", required=False, default=100)
    
    misc_args.add_argument("-p", "--passes", dest="noOfPasses",type=int,\
                        help="OPTIONAL: number of spline passes to run\
                        Default is 1", required=False, default=1)

    misc_args.add_argument("-m", "--mappabilityThres", dest="mappabilityThreshold", 
                        type=int, help="OPTIONAL: minimum number of hits \
                        per locus that has to exist to call it mappable. \
                        DEFAULT is 1.", required=False, default=1)

    misc_args.add_argument("-l", "--lib", dest="libname", help="OPTIONAL: Name \
                        of the library that is analyzed to be used for name of\
                        file prefixes. DEFAULT is OMEGA_FitHic",required=False, default='OMEGA_FitHiC')

    misc_args.add_argument("-log", action="store_true", dest="logger", \
                        help="OPTINAL: whether to writting to a log file", \
                        required=False, default=False)
    
    misc_args.add_argument("-x", "--contactType", dest="contactType",
                      help="OPTIONAL: use this flag to determine which \
                      chromosomal regions to study (All, interOnly, intraOnly). \
                      DEFAULT is intraOnly", choices=['intraOnly', 'interOnly', 'All'], 
                      required=False, default='intraOnly')
    
    misc_args.add_argument("-t", "--biases", dest="biasfile",\
                        help="OPTIONAL: biases calculated by\
                        ICE or KR norm for each locus are read from BIASFILE",\
                        required=False)
    misc_args.add_argument("-tL", "--biasLowerBound", dest="biasLowerBound",\
                      help="OPTIONAL: this flag is used to determine the lower bound\
                      of bias values to discard. DEFAULT is 0.5"\
                      , required=False, default=0.5)
    
    misc_args.add_argument("-tU", "--biasUpperBound", dest="biasUpperBound",\
                      help="OPTIONAL: this flag is used to determine the upper bound\
                      of bias values to discard. DEFAULT is 2"\
                      , required=False, default=2)

    
    
    parser.add_argument('--help', '-h', action='help', help="Print this help message and exit")
    parser.add_argument("--version", action="version", version=versionStr)

    return parser.parse_args(args)
