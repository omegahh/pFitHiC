import os
import sys
import time
import multiprocessing
from os.path import dirname, basename
try:
    from .pfithic import *
    from .parser import parse_args
except:
    from pfithic import *
    from parser import parse_args

def align_files(reads, frags):
    readsfiles = [f for f in os.listdir(reads) if f.find('count') > 0]
    fragsfiles = [f for f in os.listdir(frags) if f.find('frags') > 0]
    assert len(readsfiles) == len(fragsfiles), 'readsfile and fragsfile must be paired!'
    chr_reads = set(f.split('_')[0] for f in readsfiles)
    chr_frags = set(f.split('_')[0] for f in readsfiles)
    assert len(chr_frags-chr_reads) == 0, 'readsfile and fragsfile must be paired in chromosomes!'
    lib_reads = set(f.split('_')[1] for f in readsfiles)
    lib_frags = set(f.split('_')[1] for f in readsfiles)
    assert len(lib_frags-lib_reads) == 0, 'readsfile and fragsfile must be paired in libnames!'
    chr_list = sorted(list(chr_reads))
    lib_list = sorted(list(lib_reads))
    files_dict = {}.fromkeys(chr_list)
    for chrn in chr_list:
        temp_dict = {}.fromkeys(lib_list)
        for lib in lib_list:
            temp_dict[lib] = [os.path.join(reads, f'{chrn}_{lib}_count.gz'), os.path.join(reads, f'{chrn}_{lib}_frags.gz')]
        files_dict[chrn] = temp_dict
    return files_dict, chr_list, lib_list

def fithic_workflow(outdir, readsFile, fragsFile, resolution, distLowThres, distUpThres, noOfBins, mappThres, contactType, noOfPasses, libname, biasFile=None, biasLowerBound=0.5, biasUpperBound=2, logger=False):
    print(f'Will read {dirname(readsFile)} [{basename(readsFile)}] and [{basename(fragsFile)}]')
    if logger:
        f = open(os.path.join(outdir, f'zlog_{libname}.md'), 'w')
        sys.stdout = f
    print(f"GIVEN FIT-HI-C ARGUMENTS, LIB: {libname}")
    print("==============================")
    print(f"Output path will be: '{outdir}'")
    print(f"Reading single interactions file from: {readsFile}")
    print(f"Reading single fragments file from: {fragsFile}")
    if biasFile is not None:
        print(f"Reading bias file from {biasFile}")
        assert biasUpperBound > biasLowerBound, "Error: bias lower bound is greater than bias upper bound. Please fix."
        print(f"Bound of bias values is ({biasLowerBound}, {biasUpperBound})")
    if resolution == 0:
        print("Fixed size data not being used.")
    elif resolution > 0:
        print("Fixed size option detected, Fast version of FitHiC will be used")
        print(f"Resolution is {resolution//1000} kb")
    else:
        print("INVALID RESOLUTION ARGUMENT DETECTED")
        print("Please make sure the given resolution is a positive number greater than zero")
        print("User-given resolution: %s" % resolution)
        sys.exit(2)
    print(f"Disired Distance threshold is ({distLowThres}, {distUpThres}]")
    print(f"The maximal number of bins is {noOfBins}")
    print(f"The minimal number of reads required for valid interaction is {mappThres}")
    region_parser(contactType)
    print(f"Start FitHiC now...[{noOfPasses} Pass totally]")
    print("==============================\n")

    outliersline, outliersdist = None, None
    splineFDRx = {}.fromkeys(range(1, noOfPasses+1))
    splineFDRy = {}.fromkeys(range(1, noOfPasses+1))
    X = {}.fromkeys(range(1, noOfPasses+1))
    Y = {}.fromkeys(range(1, noOfPasses+1))
    Yerr = {}.fromkeys(range(1, noOfPasses+1))
    splineXs = {}.fromkeys(range(1, noOfPasses+1))
    splineYs = {}.fromkeys(range(1, noOfPasses+1))

    for passNo in range(1, noOfPasses+1):
        outfile = os.path.join(outdir, f'{libname}.pass{passNo}_fithic')
        splinefile = os.path.join(outdir, f'{libname}.pass{passNo}_spline')
        print(f"FIT-HI-C PROCESSING: {passNo}-TH PASS")
        print("==============================")
        allReads = read_countsFile(readsFile, distLowThres, distUpThres, outliersline)
        mainDic, observedInterAllSum, observedIntraAllSum, observedIntraInRangeSum = get_maindic(allReads, distLowThres, distUpThres)
        print("")
        binStats, bins = make_bins(mainDic, noOfBins, observedIntraInRangeSum, outliersdist)
        print("")
        allFrags = read_fragsFile(fragsFile, resolution=resolution)
        binStats, possibleIntraInRangeCount, possibleInterAllCount, interChrProb = make_fragpairs(binStats, allFrags, distLowThres, distUpThres, observedIntraInRangeSum, resolution)
        print("")
        x, y, yerr = calculateProbabilities(outfile, mainDic, binStats, observedIntraInRangeSum)
        print("")
        splineX, splineY, residual = fit_Spline(splinefile, mainDic, x, y, yerr, passNo)
        print("")
        outliesline, outliersdist, FDRx, FDRy = calculateSignificant(splinefile, readsFile, splineX, splineY, possibleIntraInRangeCount, observedIntraInRangeSum, possibleInterAllCount, observedInterAllSum, distLowThres, distUpThres, passNo, region=contactType)
        print("One pass finished!")
        print("==============================\n")
        splineFDRx[passNo] = FDRx
        splineFDRy[passNo] = FDRy
        X[passNo] = x
        Y[passNo] = y
        Yerr[passNo] = yerr
        splineXs[passNo] = splineX
        splineYs[passNo] = splineY

    if noOfPasses > 1:
        comparefdr = os.path.join(outdir, f'{libname}.{noOfPasses}passes_fdr')
        print(f'>>>> Ploting to {comparefdr}.svg')
        compare_Spline_FDR(comparefdr, splineFDRx, splineFDRy)
        comparespline = os.path.join(outdir, f'{libname}.{noOfPasses}passes_spline')
        print(f'>>>> Ploting to {comparespline}.svg')
        compareFits_Spline(comparespline, splineXs, splineYs, X, Y, Yerr)
    if logger:
        sys.stdout = sys.__stdout__
        f.close()

def main(arguments):
    pool_num = multiprocessing.cpu_count()
    args = parse_args(arguments)
    outdir = args.outdir
    reads = args.intersfile
    frags = args.fragsfile
    resolution = args.resolution
    distLowThres = args.distLowThres
    distUpThres = args.distUpThres
    noOfBins = args.noOfBins
    mappThres = args.mappabilityThreshold
    contactType = args.contactType
    noOfPasses = args.noOfPasses
    logger = args.logger

    if os.path.isfile(reads):
        os.makedirs(outdir, exist_ok=True)
        libname = args.libname
        fithic_workflow(outdir, reads, frags, resolution, distLowThres, distUpThres, noOfBins, mappThres, contactType, noOfPasses, libname, logger=logger)
    elif os.path.isdir(reads):
        files_dict, chr_list, lib_list = align_files(reads, frags)
        start = time.time()
        print(f'Start a multiprocess pool with processes = {pool_num} for pfithic analysis')
        pool = multiprocessing.Pool(processes=pool_num)
        print(f"Files: [{','.join(chr_list)}] x [{','.join(lib_list)}]")
        for chrn in chr_list:
            outdir_new = os.path.join(outdir, chrn)
            os.makedirs(outdir_new, exist_ok=True)
            for libname in lib_list:
                reads = files_dict[chrn][libname][0]
                frags = files_dict[chrn][libname][1]
                pool.apply_async(fithic_workflow, (outdir_new, reads, frags, resolution, distLowThres, distUpThres, noOfBins, mappThres, contactType, noOfPasses, libname,), {'logger': logger})
        pool.close()
        pool.join()
        print(f'All processing done. Running cost is {(time.time()-start)/60:.1f} min')

if __name__ == "__main__":
    main(sys.argv[1:])
