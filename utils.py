import numpy as np

def getType_op(distance, LowThres, UpThres):
    if distance < 0:
        contactType = 'inter'
    else:
        if (LowThres==-1 or (LowThres>-1 and distance>LowThres)) and\
           (UpThres==-1 or (UpThres>-1 and distance<=UpThres)):
            contactType = 'intraInRange'
        elif (LowThres>-1 and distance<=LowThres):
            contactType = 'intraShort'
        elif (UpThres>-1 and distance > UpThres):
            contactType = 'intraLong'
    return contactType

def getType(distances, LowThres, UpThres):
    contactTypes = []
    for distance in distances:
        contactTypes.append(getType_op(distance, LowThres, UpThres))
    return contactTypes

def getDist_op(chr1, chr2, locus1, locus2):
    if chr1 == chr2:
        return abs(locus2 - locus1)
    else:
        return -1

def getDist(*args):
    return list(map(getDist_op, *args))

def region_parser(chromosome_region='intraOnly', silence=False):
    interOnly=False
    allReg=False
    if chromosome_region == "All":
        allReg=True
        if not silence:
            print("All genomic regions will be analyzed")
    elif chromosome_region == "interOnly":
        interOnly=True
        if not silence:
            print("Only inter-chromosomal regions will be analyzed")
    elif chromosome_region == "intraOnly":
        interOnly=False
        allReg=False
        if not silence:
            print("Only intra-chromosomal regions will be analyzed")
    else:
        print("Invalid Option. Only options are 'All', 'interOnly', or 'intraOnly'")
    return allReg, interOnly

################### FUNC out_range_check  #####################################
###  Check whether the given interactionDistance is out of the range we are 
###  interested. Should only be used for intra-chromosomal interactions.
###############################################################################
def out_range_check(interactionDistance, distLowThres, distUpThres):
    if (distLowThres == -1 or (distLowThres>-1 and interactionDistance >distLowThres)) and (distUpThres ==-1 or (distUpThres>-1 and interactionDistance <= distUpThres)):
        return False
    return True


################### FUNC benjamini_hochberg_correction  #####################
### Given an array of p-values (not necessarily sorted) and the number of 
### total. Tests that were performed to gather these p-values, this function 
### performs the multiple hypothesis testing correction described by 
### Benjamini-Hochberg.
###
### If the number of tests are much more compared to p-value array and
### the omitted p-values are all supposed to be zero you should use it like:
### q_array = bh_correction([0.03,0.4,0.7,0.01],10)
### 
### If the number of tests are equal to the ones in p-values array then:
### p_array = [0.03,0.4,0.7,0.01]
### q_array = bh_correction(p_array,len(p_array))
##############################################################################
def bh_correction(p_values, num_total_tests):
    # assumes that p_values vector might not be sorted
    pvalsArray = np.array(p_values)
    order = pvalsArray.argsort()
    sorted_pvals = np.take(p_values, order)
    q_values = [1.0 for i in range(len(p_values))]
    prev_bh_value = 0
    for i, p_value in enumerate(sorted_pvals):
        if p_value == 1.0:
            bh_value = 1.0
        else:
            bh_value = p_value * num_total_tests / (i + 1)
            # Sometimes this correction can give values greater than 1,
            # so we set those values at 1
            bh_value = min(bh_value, 1)

        # To preserve monotonicity in the values, we take the
        # maximum of the previous value or this one, so that we
        # don't yield a value less than the previous.
        bh_value = max(bh_value, prev_bh_value)
        prev_bh_value = bh_value
        qIndex = order[i]
        q_values[qIndex] = bh_value
    #END for
    return q_values
