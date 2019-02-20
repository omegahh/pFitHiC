import math
import bisect
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.special import bdtrc
from scipy.interpolate import UnivariateSpline
from sklearn.isotonic import IsotonicRegression
try:
    from utils import *
    from parser import parse_args
except:
    from .utils import *
    from .parser import parse_args

toKb = 10**3
toMb = 10**6
toProb = 10**-5
# contactTypes = ['inter', 'intraInRange', 'intraShort', 'intraLong']

def read_countsFile(readsFile, lowThres, upThres, outlier=None, silence=False):
    if not silence:
        print("Reading the contact counts file...")
        print("--------------------------------------------------")
    colnames = ['chr1', 'locus1', 'chr2', 'locus2', 'contactCount']
    if readsFile.endswith('.gz'):
        allReads = pd.read_csv(readsFile, sep='\t', header=None, names=colnames,compression='gzip')
    else:
        allReads = pd.read_csv(readsFile, sep='\t', header=None, names=colnames)
    if outlier is not None:
        print(f"!!!! Discarding {len(outlier)} pairs outliers...")
        allReads = allReads.drop(outlier).reset_index(drop=True)

    allReads['distance'] = getDist(allReads.chr1, allReads.chr2, allReads.locus1, allReads.locus2)
    allReads['contactType'] = getType(allReads.distance, lowThres, upThres)
    if not silence:
        print(f"Distinguishing contactTypes according to desired genomic distance range: ({lowThres}, {upThres}]")
    return allReads

def read_fragsFile(fragsFile, resolution=0, mappThres=1):
    print("Reading the fragments mappability file...")
    print("--------------------------------------------------")
    colnames = ['chr', 'start', 'fragmid', 'mappable', 'unknown']
    if fragsFile.endswith('.gz'):
        rawFrags = pd.read_csv(fragsFile, sep='\t', header=None, names=colnames,compression='gzip')
    else:
        rawFrags = pd.read_csv(fragsFile, sep='\t', header=None, names=colnames)
    # clear fragments lower than mappThres
    rawFrags = rawFrags[rawFrags.mappable >= mappThres].reset_index()
    groupFrags = rawFrags.groupby('chr')
    allFrags = pd.DataFrame(columns=['fragmids','length','maxdist','mindist'])
    for chrn, data in groupFrags:
        fragments = sorted(data.fragmid.tolist())
        length = len(data)
        maxdist = fragments[-1] - fragments[0]
        mindist = min([fragments[i+1]-fragments[i] for i in range(length-1)])
        allFrags.loc[chrn] = [fragments, length, maxdist, mindist]
    return allFrags

def read_biasFile(biasFile, biasLowerBound, biasUpperBound):
    print("Reading the bias file...")
    if biasFile.endswith('.gz'):
        rawBias = pd.read_csv(biasFile, sep='\t', header=None, names=['chr', 'fragmid', 'bias'], compression='gzip')
    else:
        rawBias = pd.read_csv(biasFile, sep='\t', header=None, names=['chr', 'fragmid', 'bias'])

    quans = [0.05, 0.5, 0.95]
    botQ, med, topQ = rawBias[rawBias.bias != 1].bias.quantile(quans).tolist()
    print(f"5th quantile of biases: {botQ}")
    print(f"50th quantile of biases: {med}")
    print(f"95th quantile of biases: {topQ}")

    lowerIdx = rawBias[rawBias.bias<biasLowerBound].index
    upperIdx = rawBias[rawBias.bias>biasUpperBound].index
    rawBias.loc[lowerIdx, 'bias'] = -1 # botQ
    rawBias.loc[upperIdx, 'bias'] = -1 # topQ
    lowerDiscards = len(lowerIdx)
    upperDiscards = len(upperIdx)
    allBias = pd.pivot_table(rawBias, index=['chr', 'fragmid'])
    print(f"Out of {len(allBias)} total loci, {lowerDiscards+upperDiscards} were discarded with biases not in range [{biasLowerBound} {biasUpperBound}]")
    return allBias

def get_maindic(allReads, lowThres, upThres):
    dict_mapping = {'chr1': 'count', 'contactCount': 'sum'}
    pairsInfo = allReads.groupby('contactType').agg(dict_mapping)
    Infos = {'inter': 0, 'intraInRange': 0, 'intraShort': 0, 'intraLong': 0}
    Pairs = {'inter': 0, 'intraInRange': 0, 'intraShort': 0, 'intraLong': 0}
    for key in Infos.keys():
        if key in pairsInfo.index:
            Infos[key] = pairsInfo.loc[key, 'contactCount']
            Pairs[key] = pairsInfo.loc[key, 'chr1']
    observedInterAllSum = Infos['inter']
    observedIntraAllSum = Infos['intraInRange'] + Infos['intraShort'] + Infos['intraLong']
    observedIntraInRangeSum = Infos['intraInRange']
    observedInterAllCount = Pairs['inter']
    observedIntraAllCount = Pairs['intraInRange'] + Pairs['intraShort'] + Pairs['intraLong']
    observedIntraInRangeCount = Pairs['intraInRange']

    intraData = allReads[allReads.distance > 0]
    mainSeries = intraData.groupby('distance')['contactCount'].sum()
    maxObservedGenomicDist = intraData['distance'].max()
    minObservedGenomicDist = intraData['distance'].min()
    print(f"Range of desired genomic distances: ({lowThres}, {upThres}]")
    print(f"Range of observed genomic distances: [{minObservedGenomicDist}, {maxObservedGenomicDist}]")
    print(f"Observed, Intra-chr in range: pairs = {observedIntraInRangeCount}, totalCount = {observedIntraInRangeSum}")
    print(f"Observed, Intra-chr all: pairs = {observedIntraAllCount}, totalCount = {observedIntraAllSum}")
    print(f"Observed, Inter-chr all: pairs = {observedInterAllCount}, totalCount = {observedInterAllSum}")
    result = mainSeries.loc[lowThres+1:upThres] # exclude lowThres
    return result, observedInterAllSum, observedIntraAllSum, observedIntraInRangeSum

def make_bins(mainDic, nBins, observedIntraInRangeSum, outliersdist=None):
    print("Making equal occupancy bins...")
    print("--------------------------------------------------")
    noPerBin = observedIntraInRangeSum/nBins
    print(f"Observed intra-chr read counts in range {observedIntraInRangeSum}")

    n = 0 # bin counter so far
    totalContactCountSoFar = 0
    distsToGoInAbin = []
    bins = []
    desiredPerBin = (observedIntraInRangeSum)/nBins
    for dist, countsum in mainDic.iteritems():
        totalContactCountSoFar += countsum
        distsToGoInAbin.append(dist)
        if sum(mainDic.loc[distsToGoInAbin]) >= desiredPerBin:
            n += 1
            if n < nBins:
                desiredPerBin = (observedIntraInRangeSum - totalContactCountSoFar)/(nBins - n)
            bins.append(distsToGoInAbin)
            distsToGoInAbin = []
    print(f"Desired number of contacts per bin is {noPerBin} at first, bin number is {len(bins)}")
    # print(bins)
    columns = ['binLeft', 'binRight', 'possiblePairs', 'sumOverallCC', 'sumOverallDist', 'avgCC', 'avgDist', 'fragsCount', 'binContent']
    types = {'binLeft': int, 'binRight': int, 'possiblePairs': int, 'sumOverallCC': int, 'sumOverallDist': int, 'avgCC': int, 'avgDist': int, 'fragsCount': int}
    binStats = pd.DataFrame(columns=columns).astype(types)
    for i, b in enumerate(bins):
        lb = 0 if i==0 else bins[i-1][-1]+1
        # lb = lowThres+1 if i==0 else bins[i-1][-1]+1
        ub = b[-1]
        sumOverallCC = sum(mainDic.loc[b])
        binStats.loc[i] = [lb, ub, 0, sumOverallCC, 0, 0, 0, 0, b]

    if outliersdist is not None:
        print(f"!!!! Discarding {len(outliersdist)} pairs outliers...")
        binLefts = binStats.binLeft.tolist()
        binRights = binStats.binRight.tolist()
        for dist in outliersdist:
            lidx = bisect.bisect_left(binLefts, dist) - 1
            ridx = bisect.bisect_left(binRights, dist)
            assert lidx == ridx, f"lidx: {lidx} doesn't equal to ridx: {ridx} when dist = {dist}, please check bins"
            binStats.loc[lidx, 'fragsCount'] -= 1
            binStats.loc[lidx, 'possiblePairs'] -= 1

    print("Equal occupancy bins generated")
    # print(binStats)
    return binStats, bins

def make_fragpairs(binStats, allFrags, lowThres, upThres, observedIntraInRangeSum, resolution, mapThres=1):
    # scale just avoids overflow, but is necessary for large genomes
    scale = resolution if resolution>0 else 10_000

    possibleIntraAllCount = 0
    possibleInterAllCount = 0
    possibleIntraInRangeCount = 0

    N = allFrags.length.sum()
    maxPossibleGenomicDist = allFrags.maxdist.max()
    minPossibleGenomicDist = allFrags.mindist.min()

    for ch, frags in allFrags.iterrows():
        n = frags.length
        possibleIntraInRangeCountPerChr = 0
        if resolution:
            print("Looping through all possible fragment pairs in-range")
            possibleIntraAllCount += (n*(n-1))//2
            maxdist = frags.maxdist
            for i, dist in enumerate(range(0, maxdist+1, resolution)):
                if out_range_check(dist, lowThres, upThres):
                    continue
                npairs = n - i
                possibleIntraInRangeCountPerChr += npairs
                idx = binStats[(binStats.binLeft<dist)&(binStats.binRight>=dist)].index
                assert len(idx) == 1, f"{len(idx)} bin(s) founded, Expect to be exactly one bin"
                binStats.loc[idx, 'possiblePairs'] += npairs
                binStats.loc[idx, 'fragsCount'] += npairs
                binStats.loc[idx, 'sumOverallDist'] += (dist/scale)*npairs
        else:
            print("Enumerating all possible fragment pairs in-range")
            fragments = frags.fragmids
            binLefts = binStats.binLeft.tolist()
            binRights = binStats.binRight.tolist()
            possiblePairs = np.zeros(len(binStats), dtype=np.int)
            fragsCount = np.zeros(len(binStats), dtype=np.int)
            sumOverallDist = np.zeros(len(binStats))
            for x in range(n):
                for i, y in enumerate(range(x+1, n)):
                    possibleIntraAllCount += 1
                    dist = fragments[y] - fragments[x]
                    if out_range_check(dist, lowThres, upThres):
                        continue
                    possibleIntraInRangeCountPerChr += 1
                    lidx = bisect.bisect_left(binLefts, dist) - 1
                    ridx = bisect.bisect_left(binRights, dist)
                    assert lidx == ridx, f"lidx: {lidx} doesn't equal to ridx: {ridx} when dist = {dist}, please check bins"
                    possiblePairs[lidx] += 1
                    fragsCount[lidx] += 1
                    sumOverallDist[lidx] += dist/scale
            binStats.possiblePairs = possiblePairs
            binStats.fragsCount = fragsCount
            binStats.sumOverallDist = sumOverallDist
        print(f"Chromosome {ch}: \n\t{n} mappable fragments, {possibleIntraInRangeCountPerChr} possible intra-chr fragment pairs in range, {(N-n)*n//2} possible inter-chr fragment pairs")
        possibleInterAllCount += (n*(N-n))//2
        possibleIntraInRangeCount += possibleIntraInRangeCountPerChr
    # print(binStats)
    binStats['avgCC'] = binStats['sumOverallCC'] / binStats['possiblePairs'] / observedIntraInRangeSum
    binStats['avgDist'] = scale * (binStats['sumOverallDist'] / binStats['possiblePairs'])
    binStats['avgDist'] = binStats['avgDist'].astype(int)

    interChrProb = 1/possibleInterAllCount if possibleInterAllCount>0 else 0.
    intraChrProb = 1/possibleIntraAllCount

    print(f"Number of all fragments = {N}")
    print(f"Desired genomic distance range: ({lowThres} {upThres}]")
    print(f"Range of possible genomic distances: [{minPossibleGenomicDist} {maxPossibleGenomicDist}]")
    print(f"Possible, Intra-chr in range: pairs = {possibleIntraInRangeCount}")
    print(f"Possible, Intra-chr all: pairs = {possibleIntraAllCount}")
    print(f"Possible, Inter-chr all: pairs = {possibleInterAllCount}")
    print(f"Baseline intrachromosomal probability is {intraChrProb:.6e}")
    print(f"Baseline interchromosomal probability is {interChrProb:.6e}")

    return binStats, possibleIntraInRangeCount, possibleInterAllCount, interChrProb

def calculateProbabilities(outfilename, mainDic, binStats, observedIntraInRangeSum):
    print("Calculating probability means and std. errors of contact counts...")
    print("--------------------------------------------------")
    y = binStats.avgCC.values
    x = binStats.avgDist.values
    pairCounts = binStats.possiblePairs.values
    interactionTotals = binStats.sumOverallCC.values
    binContents = binStats.binContent.values

    yerr = []
    for possPairsInRange, binContent in zip(pairCounts, binContents):
        meanCountPerPair = 0
        M2 = 0
        for dist in binContent: # excluded the nonzero dist
            delta = mainDic.loc[dist] - meanCountPerPair
            meanCountPerPair += delta / possPairsInRange
            M2 += delta * (mainDic.loc[dist] - meanCountPerPair)
        var = M2 / (possPairsInRange-1)
        sd = math.sqrt(var)
        se = sd / math.sqrt(possPairsInRange)
        se_p = se / possPairsInRange / observedIntraInRangeSum
        yerr.append(se_p)
    yerr = np.array(yerr)

    print(f">>>> Writing {outfilename}.tsv")
    with open(outfilename+'.tsv', 'w') as outfile:
        outfile.write("avgGenomicDist\tcontactProbability\tstandardError\tnoOfLocusPairs\ttotalOfContactCounts\n")
        for i in range(len(x)):
            outfile.write(f"{x[i]}\t{y[i]:.4e}\t{yerr[i]:.4e}\t{pairCounts[i]}\t{interactionTotals[i]}\n")
    print(f">>>> Means and error written to {outfilename}.tsv")
    return x, y, yerr

def fit_Spline(outfilename, mainDic, x, y, yerr, passNo, visual=True):
    print("Fitting a univariate spline to the probability means...")
    print("--------------------------------------------------")
    y = np.array([f for _, f in sorted(zip(x,y), key=lambda pair: pair[0])])
    x.sort()
    # maximum residual allowed for spline is set to min(y)^2
    splineError = min(y) ** 2
    # fitpack2 method-fit on the real x and y from equal occupancy binning
    ius = UnivariateSpline(x, y, s=splineError)
    residual = sum([i*i for i in (y - ius(x))])

    splineX = np.sort(mainDic.index.values)
    if min(splineX)<min(x) or max(splineX)>max(x):
        print(f"Extrapolation happened.")
    temp_splineY = ius(splineX)

    ir = IsotonicRegression(increasing=False)
    splineY = ir.fit_transform(splineX, temp_splineY)

    if visual:
        print(f">>>> Plotting to {outfilename}.svg")
        plot_splines(outfilename, splineX, splineY, x, y, yerr, passNo)

    return splineX, splineY, residual

def calculateSignificant(outfilename, infilename, splineX, splineY, possibleIntraInRangeCount, observedIntraInRangeSum, possibleInterAllCount, observedInterAllSum, lowThres, upThres, passNo, region='intraOnly'):
    print("Calculating p-values and q-values for all pairs from input file...")
    print("--------------------------------------------------")
    newAllReads = read_countsFile(infilename, lowThres, upThres, silence=True)
    CCNT = newAllReads.contactCount.values
    DIST = newAllReads.distance.values
    ITYPE = newAllReads.contactType.values

    allReg, interOnly = region_parser(region, silence=True)

    p_vals=[]
    for cc, dist, itype in zip(CCNT, DIST, ITYPE):
        bias1 = 1.0; bias2 = 1.0 # not use bias file so far
        if (bias1 < 0 or bias2 < 0) and itype != 'inter':
            prior_p = 1.0
            p_val = 1.0
        elif itype == 'intraInRange' and not interOnly:
            i = bisect.bisect_left(splineX, dist)
            prior_p = splineY[i] * (bias1*bias2)
            p_val = bdtrc(cc-1, observedIntraInRangeSum, prior_p)
        elif itype == 'intraShort' and not interOnly:
            prior_p = 1.0
            p_val = 1.0
        elif itype == 'intraLong' and not interOnly:
            prior_p = 1.0
            p_val = 1.0
        else:
            if allReg or interOnly:
                prior_p = interChrProb * (bias1*bias2)
                p_val = bdtrc(cc-1, observedInterAllSum, prior_p)
            else:
                p_val = 1.0
        p_vals.append(p_val)

    # Do the BH FDR correction
    if allReg:
        totalValidCount = possibleIntraInRangeCount + possibleInterAllCount
    elif interOnly and not allReg:
        totalValidCount = possibleInterAllCount
    else:
        totalValidCount = possibleIntraInRangeCount
    outlierThres = 1.0 / totalValidCount
    q_vals = bh_correction(p_vals, totalValidCount)
    print(f'The calculation of p-values and q-values finished!')
    print(f'>>>> Writing to {outfilename}.significant.gz')
    newAllReads['p_vals'] = p_vals
    newAllReads['q_vals'] = q_vals
    newAllReads.to_csv(f'{outfilename}.significant.gz', sep='\t', index=False, compression='gzip')
    print(f'>>>> p-vals and q-vals written to {outfilename}.significant.gz')

    # Find all outliers
    outlierReads = newAllReads[newAllReads.p_vals < outlierThres]
    outliersline = sorted(outlierReads.index.tolist())
    outliersdist = sorted(outlierReads.distance.tolist())
    print(f'Outlier threshold is: {outlierThres:.6e}')
    print(f'Found outlier pairs: {len(outliersline)}')

    print(f'>>>> Plotting q-values to file {outfilename}.qplot.svg')
    FDRx, FDRy = plot_qvalues(outfilename, q_vals, minFDR=0, maxFDR=0.05, increment=1e-3)

    return outliersline, outliersdist, FDRx, FDRy

def plot_splines(figname, splineX, splineY, x, y, yerr, passNo):
    plt.clf()
    fig = plt.figure(figsize=[12,4])

    ax = fig.add_subplot(1,2,1)
    plt.plot(splineX/toKb, splineY/toProb, 'g-', label=f'spline-{passNo}', linewidth=2)
    plt.errorbar(x/toKb, y/toProb, yerr/toProb, fmt='r.', label="Mean with std. err", linewidth=2)
    plt.ylabel('Contact probability (x10$^{-5}$)')
    plt.xlabel('Genomic distance (Kb)')
    ax.legend(loc="upper right")

    ax = fig.add_subplot(1,2,2)
    plt.loglog(splineX, splineY, 'g-')
    plt.errorbar(x, y, yerr=yerr, fmt='r.') # Data
    plt.ylabel('Contact probability (log-scale)')
    plt.xlabel('Genomic distance (log-scale)')
    plt.savefig(figname+'.svg', format='svg')

def plot_qvalues(figname, q_values, minFDR, maxFDR, increment):
    qvalTicks = np.arange(minFDR, maxFDR+increment, increment)
    significantTicks = [0 for i in range(len(qvalTicks))]
    qvalBins = [-1 for i in range(len(q_values))]
    for i, q in enumerate(q_values):
        if math.isnan(q): q = 1 # make sure NaNs are set to 1
        qvalBins[i] = int(math.floor(q/increment))
    for b in qvalBins:
        if b >= len(qvalTicks):
            continue
        significantTicks[b] += 1

    # make it cumulative
    significantTicks = np.cumsum(significantTicks)
    # shift them by 1
    significantTicks = np.array(significantTicks[1:])
    qvalTicks = np.array(qvalTicks[1:])

    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.plot(qvalTicks, significantTicks, 'b*-')
    plt.xlabel('FDR threshold')
    plt.ylabel('Number of significant contacts')
    plt.savefig(f'{figname}.qplot.svg', format='svg')
    return qvalTicks, significantTicks

def compare_Spline_FDR(figname, splineFDRx, splineFDRy):
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    for key in sorted(splineFDRx.keys()):
        plt.plot(splineFDRx[key], splineFDRy[key]/toKb, '*-', label=f'spline-{key}')
    plt.xlabel('FDR Threshold')
    plt.ylabel('Significant contacts (x10$^{3}$)')
    ax.legend(loc='lower right')
    plt.savefig(f'{figname}.svg', format='svg')

def compareFits_Spline(figname, splineXs, splineYs, X, Y, Yerr):
    plt.clf()
    fig = plt.figure(figsize=[12, 4])

    ax = fig.add_subplot(1, 2, 1)
    for key in sorted(splineXs.keys()):
        plt.plot(splineXs[key]/toKb, splineYs[key]/toProb, label=f'spline-{key}')
        plt.errorbar(X[key]/toKb, Y[key]/toProb, Yerr[key]/toProb, fmt='r.')
    plt.ylabel('Contact probability (x10$^{-5}$)')
    plt.xlabel('Genomic distance (kb)')
    ax.legend(loc="upper right")

    ax = fig.add_subplot(1, 2, 2)
    for key in sorted(splineXs.keys()):
        plt.loglog(splineXs[key], splineYs[key], label=f'spline-{key}')
        plt.errorbar(X[key], Y[key], yerr=Yerr[key], fmt='r.')
    plt.ylabel('Contact probability (log-scale)')
    plt.xlabel('Genomic distance (log-scale)')
    ax.legend(loc='upper right')
    plt.savefig(f'{figname}.svg', format='svg')
