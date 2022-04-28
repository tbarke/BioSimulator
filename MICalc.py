import math
#from numba import jit
#import numpy as np

class MICalc(object):

    def mean(self, arr):
        mu = 0.0
        for i in arr:
            mu += i
        if len(arr) > 0:
            return mu/(len(arr))
        else:
            return 0

    def std(self, arr):
        mean = self.mean(arr)
        std = 0.0
        for i in arr:
            std += math.pow(i - mean, 2)
        return math.sqrt(std)

    def skewness(self, arr):
        n = len(arr)
        mean = self.mean(arr)
        std = self.std(arr) + 0.000001
        # E[X-mu/std)^3]
        sum = 0.0
        for i in arr:
            sum += math.pow(((i - mean)/std),3)
        return sum/n

    def findSigma(self, arr):
        n = len(arr)
        if n < 3:
            n = 3
        return math.sqrt((6*(n-2))/((n+1)*(n+3)))

    def doanes(self, arr):
        if len(arr) < 1:
            return 1
        sigma = self.findSigma(arr)
        skew = math.fabs(self.skewness(arr))
        n = len(arr)
        return math.ceil(1 +math.log(n,2) + math.log((1+(skew/sigma)), 2))

    def bins(self, k, arr):
        retArr = []
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
        for i in range(k+1):
            retArr.append(w*i + low)
        return retArr

    def getw(self,arr,k):
        low = arr[0]
        high = arr[0]
        for i in arr:
            if i < low:
                low = i
            if i > high:
                high = i
        return (high-low)/(k)

    def __init__(self):
        return

    def binPosition1(self, A, bins):
        Apos = 0
        # print(A)
        # print(B)
        for i in range(1, len(bins)):
            if A < bins[i]:
                if A > bins[i - 1]:
                    Apos = i - 1
        return Apos

    def binPosition(self, A, B):
        Apos = 0
        Bpos = 0
        #print(A)
        #print(B)
        for i in range(1,len(self.Albins)):
            if A < self.Albins[i]:
                if A > self.Albins[i-1]:
                    Apos = i - 1

        for i in range(1,len(self.Blbins)):
            if B < self.Blbins[i]:
                if B > self.Blbins[i-1]:
                    Bpos = i - 1
        return [Apos,Bpos]

    def binPosition4(self,Al, Bl, Ar, Br, Albins, Blbins, Arbins, Brbins):
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

    def binarrSingle(self, data, datak, fullarr):
        singleBins = self.bins(datak, fullarr)
        retArr = []
        for i in range(datak):
            retArr.append(0.0)
        allPoints = float(len(data))
        for i in range(len(data)):
            binPos = self.binPosition1(data[i], singleBins)
            retArr[binPos] += 1.0/allPoints
        return retArr

    def binArr(self, Als, Bls, Ars, Brs, Albins, Blbins, Arbins, Brbins):
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
            binPos = self.binPosition4(Als[i], Bls[i], Ars[i], Brs[i], Albins, Blbins, Arbins, Brbins)
            binnedABlrs[binPos[0]][binPos[1]][binPos[2]][binPos[3]] += 1.0/allPoints

        return binnedABlrs

    #currently used H(X) calculation
    def returnSingleEntropy(self, data, fullArr):
        if len(data) < 1:
            return 0
        datak = 20
        binned = self.binarrSingle(data, datak, fullArr)
        #width = self.getw(data,datak)
        sum = 0.0
        for i in range(len(binned)):
            if binned[i] == 1.0:
                continue
            if binned[i] > 0:
                sum += -1*(binned[i] * math.log(binned[i], 2))# + math.log(width,2))
        return sum

    def returnEntropy(self, pdf):
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

    def binProbability(self, bins, arr, k):
        probabilities = [0.0] * k
        length = float(len(arr))
        for i in arr:
            for j in range(1, len(bins)):
                if i < bins[j]:
                    if i > bins[j-1]:
                        probabilities[j-1] += 1.0/length
        return probabilities

    #currently used H(X|Y) calculation
    def singleEntropyXgivenY(self, x, y):
        xk = 30
        yk = 30
        binnedY = self.bins(yk, y)
        allArr = []
        for i in range(yk):
            allArr.append([])

        for i in range(len(x)):
            ypos = self.binPosition1(y[i], binnedY)
            allArr[ypos].append(x[i])

        Yprobs = self.binProbability(binnedY, y, yk)
        ret = 0.0
        for i in range(len(allArr)):
            ret += self.returnSingleEntropy(allArr[i], x)*Yprobs[i]
        return ret

    def entropyABlrgivenZ(self, xal, xbl, xar, xbr, yal, ybl, yar, ybr, XAlbins, XBlbins, XArbins, XBrbins, YAlbins, YBlbins, YArbins, YBrbins):
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
            Al, Bl, Ar, Br = self.binPosition4(xal[i], xbl[i], xar[i], xbr[i], XAlbins, XBlbins, XArbins, XBrbins)
            allZBins[Al][Bl][Ar][Br][0].append(yal[i])
            allZBins[Al][Bl][Ar][Br][1].append(ybl[i])
            allZBins[Al][Bl][Ar][Br][2].append(yar[i])
            allZBins[Al][Bl][Ar][Br][3].append(ybr[i])

        XPdfs = self.binArr(xal, xbl, xar, xbr, XAlbins, XBlbins, XArbins, XBrbins)
        sum = 0.0
        for i in range(len(allZBins)-1):
            for j in range(len(allZBins[i])-1):
                for k in range(len(allZBins[i][j])-1):
                    for l in range(len(allZBins[i][j][k])-1):
                        yAls, yBls, yArs, yBrs = allZBins[i][j][k][l]
                        entr_Als = self.returnSingleEntropy(yAls, yal)
                        entr_Bls = self.returnSingleEntropy(yBls, ybl)
                        entr_Ars = self.returnSingleEntropy(yArs, yar)
                        entr_Brs = self.returnSingleEntropy(yBrs, ybr)
                        newEntropy = entr_Als + entr_Bls + entr_Ars + entr_Brs
                        sum += newEntropy*(XPdfs[i][j][k][l])
        return sum

    def entropyABgivenZ(self):
        Zarr =self.Zbins
        Als = self.Als
        Bls = self.Bls
        Zs = self.Zs
        allPDFs = []
        for i in range(self.Zk):
            fullPDF = []
            for j in range(self.Alk):
                currArr = [0.0]*(self.blk)
                fullPDF.append(currArr)
            allPDFs.append(fullPDF)
        ZbinnedAValues = []
        ZbinnedBvalues = []
        for i in range(self.Zk):
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

        zBinProb = self.binProbability(Zarr,Zs,self.ZArs)
        sum = 0.0
        for i in range(len(ZbinnedAValues)):
            binnedABls = self.binArr(ZbinnedAValues[i], ZbinnedBvalues[i])
            sum += self.returnEntropy(binnedABls) * zBinProb[i]
        return sum

    def MI(self, Als, Bls, Ars, Brs):
        fullPDF = self.binArr(Als, Bls, Ars, Brs)
        HAB = self.returnEntropy(fullPDF)
        HABgivenZ = self.entropyABlrgivenZ()
        return HAB - HABgivenZ

    #assuming data is independent array of data
    def AltMI(self, dataX, dataY):
        HXGY = 0.0
        tot_entr = 0.0
        tot_entrXY = 0.0
        for i in range(len(dataX)):
            entr = self.returnSingleEntropy(dataX[i], dataX[i])
            tot_entr = tot_entr + entr
            entrY = self.singleEntropyXgivenY(dataX[i], dataY[i])
            tot_entrXY = tot_entrXY + entrY
            HXGY += entr - entrY
        return [HXGY, tot_entr, tot_entrXY]

    def entropy_allX(self, dataX):
        entr = 0.0
        for i in range(len(dataX)):
            entr = entr + self.returnSingleEntropy(dataX[i], dataX[i])
        return entr

    def AltMI2(self, dataX, dataY):
        HXGY = 0.0
        x = dataX[0]
        tot_entr = 0.0
        tot_entrXY = 0.0
        for i in range(len(dataY)):
            entr = self.returnSingleEntropy(x)
            tot_entr = tot_entr + entr
            entrY = self.singleEntropyXgivenY(x, dataY[i])
            tot_entrXY = tot_entrXY + entrY
            HXGY += entr - entrY
        return [HXGY, tot_entr, tot_entrXY]

    #MI of a dependent input and ouput set
    def dependentMI(self, xal, xbl, xar, xbr, yal, ybl, yar, ybr):
        k = 30
        XAlbins = self.bins(k, xal)
        XBlbins = self.bins(k, xbl)
        XArbins = self.bins(k, xar)
        XBrbins = self.bins(k, xbr)

        YAlbins = self.bins(k, yal)
        YBlbins = self.bins(k, ybl)
        YArbins = self.bins(k, yar)
        YBrbins = self.bins(k, ybr)
        Ybinned4dim = self.binArr(yal, ybl, yar, ybr, YAlbins, YBlbins, YArbins, YBrbins)
        HY = self.returnEntropy(Ybinned4dim)
        HYgivenX = self.entropyABlrgivenZ(xal, xbl, xar, xbr, yal, ybl, yar, ybr, XAlbins, XBlbins, XArbins, XBrbins,YAlbins, YBlbins, YArbins, YBrbins)
        return HY - HYgivenX