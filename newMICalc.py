import math
from numba import jit
import numpy as np

def mean(arr):
    mu = 0.0
    for i in arr:
        mu += i
    if len(arr) > 0:
        return mu/(len(arr))
    else:
        return 0

def std(arr):
    mean = mean(arr)
    std = 0.0
    for i in arr:
        std += math.pow(i - mean, 2)
    return math.sqrt(std)

def skewness(arr):
    n = len(arr)
    mean = mean(arr)
    std = std(arr) + 0.000001
    # E[X-mu/std)^3]
    sum = 0.0
    for i in arr:
        sum += math.pow(((i - mean)/std),3)
    return sum/n

def findSigma(arr):
    n = len(arr)
    if n < 3:
        n = 3
    return math.sqrt((6*(n-2))/((n+1)*(n+3)))

def doanes(arr):
    if len(arr) < 1:
        return 1
    sigma = findSigma(arr)
    skew = math.fabs(skewness(arr))
    n = len(arr)
    return math.ceil(1 +math.log(n,2) + math.log((1+(skew/sigma)), 2))

@jit(nopython=True)
def bins( k, arr, retArr):
    low = arr[0]
    high = arr[0]
    for i in arr:
        if i < low:
            low = i
        if i > high:
            high = i
    low = low - 0.001 #minus epsilon
    high = high + 0.001 #plus epsilon
    w = (high - low)/k
    for i in range(k[0]+1):
        retArr = np.append(retArr, w*i + low)
    return retArr

@jit(nopython=True)
def getw(arr,k, empty):
    low = arr[0]
    high = arr[0]
    for i in arr:
        if i < low:
            low = i
        if i > high:
            high = i
    return np.append(empty,(high-low)/(k[0]))

def binPosition1( A, bins):
    Apos = 0
    # print(A)
    # print(B)
    for i in range(1, len(bins)):
        if A < bins[i]:
            if A > bins[i - 1]:
                Apos = i - 1
    return Apos

def binPosition( A, B):
    Apos = 0
    Bpos = 0
    #print(A)
    #print(B)
    for i in range(1,len(Albins)):
        if A < Albins[i]:
            if A > Albins[i-1]:
                Apos = i - 1
    for i in range(1,len(Blbins)):
        if B < Blbins[i]:
            if B > Blbins[i-1]:
                Bpos = i - 1
    return [Apos,Bpos]

def binPosition4(Al, Bl, Ar, Br, Albins, Blbins, Arbins, Brbins):
    Alpos = 0
    Blpos = 0
    Arpos = 0
    Brpos = 0
    for i in range(1, len(Albins)):
        if Al < Albins[i]:
            if Al > Albins[i - 1]:
                Alpos = i - 1
    for i in range(1, len(Blbins)):
        if Bl < Blbins[i]:
            if Bl > Blbins[i - 1]:
                Blpos = i - 1
    for i in range(1, len(Arbins)):
        if Ar < Arbins[i]:
            if Ar > Arbins[i - 1]:
                Arpos = i - 1
    for i in range(1, len(Brbins)):
        if Br < Brbins[i]:
            if Br > Brbins[i - 1]:
                Brpos = i - 1
    return [Alpos, Blpos, Arpos, Brpos]

@jit(nopython=True)
def binarrSingle(data, datak, fullarr):
    singleBins = bins(datak, fullarr)
    retArr = []
    for i in range(datak):
        retArr.append(0.0)
    allPoints = float(len(data))
    for i in range(len(data)):
        binPos = binPosition1(data[i], singleBins)
        retArr[binPos] += 1.0/allPoints
    return retArr

def binArr( Als, Bls, Ars, Brs, Albins, Blbins, Arbins, Brbins):
    binnedABlrs = []
    for i in range(len(Albins)-1):
        currArrBl = []
        binnedABlrs.append(currArrBl)
        for j in range(len(Blbins)-1):
            currArrAr = []
            currArrBl.append(currArrAr)
            for k in range(len(Arbins)-1):
                currArrBr = [0.0]*(len(Brbins)-1)
                currArrAr.append(currArrBr)
    allPoints = float(len(Als))
    for i in range(len(Als)):
        binPos = binPosition4(Als[i], Bls[i], Ars[i], Brs[i], Albins, Blbins, Arbins, Brbins)
        binnedABlrs[binPos[0]][binPos[1]][binPos[2]][binPos[3]] += 1.0/allPoints
    return binnedABlrs

@jit(nopython=True)
def returnSingleEntropy( data, fullArr,k, empty):
    if len(data) < 1:
        return 0
    binned = binarrSingle(data, k[0], fullArr)
    #width = getw(data,datak)
    sum = 0.0
    for i in range(len(binned)):
        if binned[i] == 1.0:
            continue
        if binned[i] > 0:
            sum += -1*(binned[i] * math.log(binned[i], 2))# + math.log(width,2))
    return np.append(empty, sum)

def returnEntropy( pdf):
    sum = 0.0
    pdfSum = 0.0
    for i in range(len(pdf)):
        currArrBl = pdf[i]
        #print("new line")
        for j in range(len(currArrBl)):
            currArrAr = currArrBl[j]
            for k in range(len(currArrAr)):
                currArrBr = currArrAr[k]
                for l in range(len(currArrBr)):
                    currPDF = currArrBr[l]
                    pdfSum += currPDF
                    if currPDF > 0:
                        sum += -1 * (currPDF*math.log(currPDF, 2))
    ret = sum
    return ret
def binProbability( bins, arr, k):
    probabilities = [0.0] * k
    length = float(len(arr))
    for i in arr:
        for j in range(1, len(bins)):
            if i < bins[j]:
                if i > bins[j-1]:
                    probabilities[j-1] += 1.0/length
    return probabilities
def singleEntropyXgivenY( x, y):
    xk = doanes(x)
    yk = doanes(y)
    xk = 20
    yk = 20
    binnedX = bins(xk, x)
    binnedY = bins(yk, y)
    allArr = []
    for i in range(yk):
        allArr.append([])
    for i in range(len(x)):
        ypos = binPosition1(y[i], binnedY)
        allArr[ypos].append(x[i])
    Yprobs = binProbability(binnedY, y, yk)
    ret = 0.0
    for i in range(len(allArr)):
        ret += returnSingleEntropy(allArr[i], x)*Yprobs[i]
    return ret

#data is an independent array of variables
def entropyXgivenY( dataX, dataY):
    sum = 0.0
    for i in range(len(dataX)):
        sum += singleEntropyXgivenY(dataX[i], dataY[i])
    return sum

def entropyABlrgivenZ( xal, xbl, xar, xbr, yal, ybl, yar, ybr, XAlbins, XBlbins, XArbins, XBrbins, YAlbins, YBlbins, YArbins, YBrbins):
    allZBins = []
    for i in range(len(XAlbins)):
        iarr = []
        allZBins.append(iarr)
        for j in range(len(XBlbins)):
            jarr = []
            iarr.append(jarr)
            for k in range(len(XArbins)):
                currArr = []
                jarr.append(currArr)
                for l in range(len(XBrbins)):
                    currArr.append([[],[],[],[]])
    for i in range(len(yal)):
        Al, Bl, Ar, Br = binPosition4(xal[i], xbl[i], xar[i], xbr[i], XAlbins, XBlbins, XArbins, XBrbins)
        allZBins[Al][Bl][Ar][Br][0].append(yal[i])
        allZBins[Al][Bl][Ar][Br][1].append(ybl[i])
        allZBins[Al][Bl][Ar][Br][2].append(yar[i])
        allZBins[Al][Bl][Ar][Br][3].append(ybr[i])
    XPdfs = binArr(xal, xbl, xar, xbr, XAlbins, XBlbins, XArbins, XBrbins)
    sum = 0.0
    for i in range(len(allZBins)-1):
        for j in range(len(allZBins[i])-1):
            for k in range(len(allZBins[i][j])-1):
                for l in range(len(allZBins[i][j][k])-1):
                    yAls, yBls, yArs, yBrs = allZBins[i][j][k][l]
                    entr_Als = returnSingleEntropy(yAls, yal)
                    entr_Bls = returnSingleEntropy(yBls, ybl)
                    entr_Ars = returnSingleEntropy(yArs, yar)
                    entr_Brs = returnSingleEntropy(yBrs, ybr)
                    newEntropy = entr_Als + entr_Bls + entr_Ars + entr_Brs
                    sum += newEntropy*(XPdfs[i][j][k][l])
    return sum

def entropyABgivenZ(self):
    Zarr =Zbins
    Als = Als
    Bls = Bls
    Zs = Zs
    allPDFs = []
    for i in range(Zk):
        fullPDF = []
        for j in range(Alk):
            currArr = [0.0]*(blk)
            fullPDF.append(currArr)
        allPDFs.append(fullPDF)
    ZbinnedAValues = []
    ZbinnedBvalues = []
    for i in range(Zk):
        currArr1 = []
        currArr2 = []
        ZbinnedAValues.append(currArr1)
        ZbinnedBvalues.append(currArr2)
    for i in range(len(Als)):
        for j in range(1,len(Zarr)):
            if Zs[i] < Zarr[j]:
                if Zs[i] > Zarr[j - 1]:
                    ZbinnedAValues[j-1].append(Als[i])
                    ZbinnedBvalues[j-1].append(Bls[i])
    zBinProb = binProbability(Zarr,Zs,ZArs)
    sum = 0.0
    for i in range(len(ZbinnedAValues)):
        binnedABls = binArr(ZbinnedAValues[i], ZbinnedBvalues[i])
        sum += returnEntropy(binnedABls) * zBinProb[i]
    return sum

def MI( Als, Bls, Ars, Brs):
    fullPDF = binArr(Als, Bls, Ars, Brs)
    HAB = returnEntropy(fullPDF)
    HABgivenZ = entropyABlrgivenZ()
    return HAB - HABgivenZ

#assuming data is independent array of data
@jit(nopython=True)
def AltMI( dataX, dataY, k, empty):
    HXGY = 0.0
    for i in range(len(dataX)):
        entr = returnSingleEntropy(dataX[i], dataX[i], k, empty)
        entrY = singleEntropyXgivenY(dataX[i], dataY[i], k, empty)
        HXGY += entr - entrY
    return HXGY

def AltMI2(dataX, dataY):
    HXGY = 0.0
    x = dataX[0]
    for i in range(len(dataY)):
        entr = returnSingleEntropy(x)
        entrY = singleEntropyXgivenY(x, dataY[i])
        HXGY += entr - entrY
    return HXGY

@jit(nopython=True)
def dependentMI( xal, xbl, xar, xbr, yal, ybl, yar, ybr, empty, k):
    XAlbins = bins(k, xal, empty)
    XBlbins = bins(k, xbl, empty)
    XArbins = bins(k, xar, empty)
    XBrbins = bins(k, xbr, empty)
    YAlbins = bins(k, yal, empty)
    YBlbins = bins(k, ybl, empty)
    YArbins = bins(k, yar, empty)
    YBrbins = bins(k, ybr, empty)
    wal = getw(yal,k, empty)
    war = getw(ybl,k, empty)
    wbl = getw(yar,k, empty)
    wbr = getw(ybr,k, empty)
    Ybinned4dim = binArr(yal, ybl, yar, ybr, YAlbins, YBlbins, YArbins, YBrbins)
    HY = returnEntropy(Ybinned4dim)# + math.log(wal,2) + math.log(war,2) + math.log(wbl,2) + math.log(wbr,2)
    #print(HY)
    HYgivenX = entropyABlrgivenZ(xal, xbl, xar, xbr, yal, ybl, yar, ybr, XAlbins, XBlbins, XArbins, XBrbins,YAlbins, YBlbins, YArbins, YBrbins)
    #print(HYgivenX)
    return HY - HYgivenX
#n = number of receptors
#pX = distribution of concentraiton
#a = probability scaling factor
#probability of binding is equal to a*x/n
#def expected( n , px):