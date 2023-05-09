import colors
import utils
import math
import MICalc
import glob
import matplotlib.pyplot as plt
import numpy as np
import utils
import log
l = log.log()


def reduceSize(arr, k):
    new_arr = []
    for i in range(0, len(arr), k):
        new_arr.append(arr[i])
    return new_arr

def calcGrowth(file, config):
    growth_dat = utils.loadDataDate(file, False)
    growth = growth_dat[0]
    avg_growth = 0.0
    for i in range(0, len(growth) - 1):
        if growth[i] > 0 and growth[i + 1] > 0:
            avg_growth += math.log2(growth[i + 1] / growth[i])
    avg_growth = (avg_growth / len(growth)) * (1 / config.simParams.simTimeStep)
    return avg_growth

def MIMoveWeighted(concs, moves, inters, binsMI, binsintAB, oneD):
    # list of concs: AR, AL, BR, BL
    # list of moves: absolute velocity
    # list of internal states intA, intB
    # find max intA and intA
    print(len(moves[0]))

    def findMI(indexes):
        currConcsAR = []
        currConcsAL = []
        currConcsBR = []
        currConcsBL = []
        currMoves = []
        for ind in indexes:
            currConcsAR.append(concs[0][ind])
            currConcsAL.append(concs[1][ind])
            currConcsBR.append(concs[2][ind])
            currConcsBL.append(concs[3][ind])

            currMoves.append(moves[0][ind])
        MImove1 = utils.MI2D1D(currConcsAR, currConcsAL, currMoves, binsMI)[0]
        MImove2 = utils.MI2D1D(currConcsBR, currConcsBL, currMoves, binsMI)[0]
        return MImove1, MImove2


    def findBins(intA, intB):
        if not oneD:
            maxintA = max(intA)
            maxintB = max(intB)
            minintA = min(intA)
            minintB = min(intB)
            binWidthA = (maxintA-minintA)/binsintAB
            binWidthB = (maxintB-minintB)/binsintAB
            binedgesA = utils.list_empty([binsintAB+1], 0)
            binedgesB = utils.list_empty([binsintAB+1], 0)
            binAB = utils.list_empty([binsintAB, binsintAB], [])
            for i in range(binsintAB+1):
                binedgesA[i] = minintA + binWidthA * i
                binedgesB[i] = minintB + binWidthB * i
            binedgesA[0] = binedgesA[0] - 0.01
            binedgesB[0] = binedgesB[0] - 0.01

            binedgesA[binsintAB] = binedgesA[binsintAB] + 0.01
            binedgesB[binsintAB] = binedgesB[binsintAB] + 0.01
            return binedgesA, binedgesB, binAB
        else:
            intRatio = []
            for i in range(len(intA)):
                if intA[i] + intB[i] >0:
                    intRatio.append(intA[i]/(intA[i] + intB[i]))
                else:
                    intRatio.append(0.0)
            maxintABrat = max(intRatio)
            minintABrat = min(intRatio)
            binWidthArat = (maxintABrat - minintABrat) / binsintAB
            binedgesArat = utils.list_empty([binsintAB + 1], 0)
            binAB = utils.list_empty([binsintAB], [])
            for i in range(binsintAB+1):
                binedgesArat[i] = minintABrat + binWidthArat * i
            binedgesArat[0] = binedgesArat[0] - 0.01
            binedgesArat[binsintAB] = binedgesArat[binsintAB] + 0.01
            return binedgesArat, binedgesArat, binAB

    # divide into total bins used
    binedgeA, binedgeB, binAB = findBins(inters[0], inters[1])
    #print(binedgeA)
    #print(binedgeB)

    def findBin(intA, intB):
        if not oneD:
            for i in range(len(binedgeA)):
                if intA >= binedgeA[i] and intA < binedgeA[i + 1]:
                    for j in range(len(binedgeA)):
                        if intB >= binedgeB[j] and intB < binedgeB[j +1]:
                            return i, j
        else:
            if (intA+intB) > 0:
                intRat = intA/(intA+intB)
            else:
                intRat = 0.0
            for i in range(len(binedgeA)):
                if intRat >= binedgeA[i] and intRat < binedgeA[i+1]:
                    return i, i

    for i in range(len(inters[0])):
        curr_intA = inters[0][i]
        curr_intB = inters[1][i]
        indexA, indexB = findBin(curr_intA, curr_intB)
        if not oneD:
            binAB[indexA][indexB].append(i)
        else:
            binAB[indexA].append(i)

    if not oneD:
        binProbs = utils.list_empty([binsintAB, binsintAB], 0)
    else:
        binProbs = utils.list_empty([binsintAB], 0)
    for i in range(len(binAB)):
        if oneD:
            binProbs[i] = len(binAB[i])/len(moves[0])
            continue
        for j in range(len(binAB[0])):
            binProbs[i][j] = len(binAB[i][j])/len(moves[0])

    MIweight = 0
    if not oneD:
        binABind = utils.list_empty([binsintAB, binsintAB], [])
    else:
        binABind = utils.list_empty([binsintAB], [])
    for i in range(len(binAB)):
        if oneD:
            currMIA, currMIB = findMI(binAB[i])
            binABind[i] = [currMIA, currMIB]
            MIweight += (currMIA + currMIB) * binProbs[i]
            continue
        for j in range(len(binAB[0])):
            currMIA, currMIB = findMI(binAB[i][j])
            binABind[i][j] = [currMIA, currMIB]
            MIweight += (currMIA+currMIB)*binProbs[i][j]
    return MIweight, binABind

def MI_trad(ext, bound, k, bins):
    dataX = []
    dataY = []
    MI_obj = MICalc.MICalc()
    dataX.append(reduceSize(ext[0], k))
    dataX.append(reduceSize(ext[1], k))
    dataX.append(reduceSize(ext[2], k))
    dataX.append(reduceSize(ext[3], k))
    dataY.append(reduceSize(bound[0], k))
    dataY.append(reduceSize(bound[1], k))
    dataY.append(reduceSize(bound[2], k))
    dataY.append(reduceSize(bound[3], k))
    MI = MI_obj.AltMI(dataX, dataY, bins)[0]
    return MI

def CalcData(config, bins, MITradFlag, MI2d2dFlag, MIMoveFlag, growthFlag, intWeightFlag, ABintdynamicWeightFlag, intABFile = None, extABFile = None, moveFile = None, boundFile = None, totalCellsFile = None, recsFile = None):
    extAB = None
    move = None
    boundAB = None
    intAB = None
    if MITradFlag or MI2d2dFlag or MIMoveFlag or intWeightFlag or ABintdynamicWeightFlag:
        extAB = utils.loadDataDate(extABFile, False)
    if MIMoveFlag or intWeightFlag or ABintdynamicWeightFlag:
        move = utils.loadDataDate(moveFile, False)
    if MITradFlag or MI2d2dFlag or MIMoveFlag or intWeightFlag or ABintdynamicWeightFlag:
        boundAB = utils.loadDataDate(boundFile, False)

    if intWeightFlag or ABintdynamicWeightFlag:
        intAB = utils.loadDataDate(intABFile, False)
        recs = utils.loadDataDate(recsFile, False)

    ext_A = []
    ext_B = []
    move_arr = []

    extAL = []
    extAR = []
    extBL = []
    extBR = []
    boundAL = []
    boundAR = []
    boundBL = []
    boundBR = []

    boundA = []
    boundB = []

    count = 0
    #intA = []
    intAextAL = []
    intAextAR = []

    #intB = []
    intBextBL = []
    intBextBR = []

    intBMove = []
    intAMove = []

    intAall =[]
    intBall = []

    Arecs = []
    Brecs = []

    if MITradFlag or MI2d2dFlag or MIMoveFlag or intWeightFlag or ABintdynamicWeightFlag:
        for i in range(len(extAB[0])):
            for j in range(len(extAB[0][i])):
                count += 1
                if MIMoveFlag or intWeightFlag or ABintdynamicWeightFlag:
                    move_arr.append(move[0][i][j])
                if MI2d2dFlag or MIMoveFlag or intWeightFlag or ABintdynamicWeightFlag:
                    ext_A.append(extAB[0][i][j][1] - extAB[0][i][j][0])
                    ext_B.append(extAB[0][i][j][3] - extAB[0][i][j][2])
                if MI2d2dFlag or intWeightFlag or ABintdynamicWeightFlag:
                    boundAB_all = math.fabs(boundAB[0][i][j][1] - boundAB[0][i][j][0]) + math.fabs(boundAB[0][i][j][3] - boundAB[0][i][j][2])
                    if boundAB_all > 0:
                        boundA.append((boundAB[0][i][j][1] - boundAB[0][i][j][0])/boundAB_all)
                        boundB.append((boundAB[0][i][j][3] - boundAB[0][i][j][2])/boundAB_all)
                    else:
                        boundA.append(0)
                        boundB.append(0)

                if MITradFlag or MIMoveFlag or MI2d2dFlag or intWeightFlag or ABintdynamicWeightFlag:
                    extAL.append(extAB[0][i][j][0])
                    extAR.append(extAB[0][i][j][1])
                    #------------bound A------------
                    boundAL.append(boundAB[0][i][j][0])
                    boundAR.append(boundAB[0][i][j][1])
                    #------------------------------
                    extBL.append(extAB[0][i][j][2])
                    extBR.append(extAB[0][i][j][3])
                    # ------------bound A------------
                    boundBL.append(boundAB[0][i][j][2])
                    boundBR.append(boundAB[0][i][j][3])
                if intWeightFlag or ABintdynamicWeightFlag:
                    intA_curr = intAB[0][i][j][0]
                    intB_curr = intAB[0][i][j][1]
                    Arec_curr = recs[0][i][j][0]
                    Brec_curr = recs[0][i][j][1]
                    intAall.append(intA_curr)
                    intBall.append(intB_curr)
                    Arecs.append(Arec_curr)
                    Brecs.append(Brec_curr)
                    if intA_curr > intB_curr:
                        intBextBL.append(extAB[0][i][j][2])
                        intBextBR.append(extAB[0][i][j][3])
                        intBMove.append(move[0][i][j])
                    else:
                        intAextAL.append(extAB[0][i][j][0])
                        intAextAR.append(extAB[0][i][j][1])
                        intAMove.append(move[0][i][j])


    MI2D2D = None
    MITrad = None
    MImove = None
    intweightMIMove = None
    weightedMIind = None
    growth = None

    k = 1
    if count > 20000:
        print("this happened")
        k = math.ceil(count / 10000)

    def reduce(arr, k):
        new = []
        for i in range(0, len(arr), k):
            new.append(arr[i])
        return new

    """""

    Arecs = reduce(Arecs, k)
    Brecs = reduce(Brecs, k)

    extAR = reduce(extAR, k)
    extAL = reduce(extAL, k)
    extBR = reduce(extBR, k)
    extBL = reduce(extBL, k)

    boundAR = reduce(boundAR, k)
    boundAL = reduce(boundAL, k)
    boundBR = reduce(boundBR, k)
    boundBL = reduce(boundBL, k)

    move_arr = reduce(move_arr, k)

    intAall = reduce(intAall, k)
    intBall = reduce(intBall, k)

    intAextAL = reduce( intAextAL, k)
    intAextAR = reduce(intAextAR , k)
    intAMove = reduce(intAMove , k)

    intBextBL = reduce(intBextBL, k)
    intBextBR = reduce(intBextBR, k)
    intBMove = reduce(intBMove, k)
    """""

    if MITradFlag:
        MITrad = MI_trad([extAL, extAR, extBL, extBR], [boundAL, boundAR, boundBL, boundBR], k, bins)
    if MI2d2dFlag:

        MI2D2D1 = utils.MI2D2D(extAR, extAL, boundAR, boundAL, bins)[0]
        MI2D2D2 = utils.MI2D2D(extBR, extBL, boundBR, boundBL, bins)[0]
        MI2D2D = MI2D2D1 + MI2D2D2
    if MIMoveFlag:
        #TODO fix old MI move code
        #MImove1 = utils.MI2D1D(extAR, extAL, move_arr, bins)[0]
        #MImove2 = utils.MI2D1D(extBR, extBL, move_arr, bins)[0]
        #MImove = MImove1 + MImove2
        ARboundNorm = []
        ALboundNorm = []
        BRboundNorm = []
        BLboundNorm = []
        for i in range(len(boundAR)):
            if Arecs[i] > 0:
                ARboundNorm.append(boundAR[i]/(Arecs[i]/2))
                ALboundNorm.append(boundAL[i]/(Arecs[i]/2))
            else:
                ARboundNorm.append(0.0)
                ALboundNorm.append(0.0)

            if Brecs[i] > 0:
                BRboundNorm.append(boundBR[i] / (Brecs[i] / 2))
                BLboundNorm.append(boundBL[i] / (Brecs[i] / 2))
            else:
                BRboundNorm.append(0.0)
                BLboundNorm.append(0.0)
        MIsyntA = utils.MI2D2D(extAR, extAL, ARboundNorm, ALboundNorm, bins)[0]
        MIsyntB = utils.MI2D2D(extBR, extBL, BRboundNorm, BLboundNorm, bins)[0]
        MImove = MIsyntA + MIsyntB
    if growthFlag:
        growth = calcGrowth(totalCellsFile, config)
    if intWeightFlag:
        intWeightAMI = utils.MI2D1D(intAextAL, intAextAR, intAMove, bins)[0]
        intWeightBMI = utils.MI2D1D(intBextBL, intBextBR, intBMove, bins)[0]
        intweightMIMove = intWeightAMI + intWeightBMI
    if ABintdynamicWeightFlag:
        #concs, moves, inters, binsMI, binsintAB
        concs = [extAR, extAL, extBR, extBL]
        BoundRecs = [boundAR, boundAL, boundBR, boundBL]
        recs = [Arecs, Brecs]
        moves = [move_arr]
        inters = [intAall, intBall]
        weightedMI, binMI = MIMoveWeighted(concs, moves, inters, 30, 10, False)
        weightedMIind = weightedMI


    return MITrad, MI2D2D, MImove, growth, intweightMIMove, weightedMIind

def getCalcArray(config, run, date, stresses, strats, enviorns, bins):
    MI_tradAll = []
    MI2D2DAll = []
    MI2D1DAll = []
    growthAll = []
    MI2D1DinAll = []
    for i in range(len(strats)):
        MI_tradAll.append([])
        MI2D2DAll.append([])
        MI2D1DAll.append([])
        growthAll.append([])
        MI2D1DinAll.append([])
        for j in range(len(enviorns)):
            MI_tradAll[i].append([])
            MI2D2DAll[i].append([])
            MI2D1DAll[i].append([])
            growthAll[i].append([])
            MI2D1DinAll[i].append([])
            for k in range(len(stresses)):
                l.log(strats[i])
                l.log(enviorns[j])
                l.log(stresses[k])
                externalAB_raw = glob.glob(
                    'Data/' + date + '/' + run + '/' + str(stresses[k]) + '_' + strats[i] + '_' + enviorns[
                        j] + '/totalConcs*.gz')[0]
                movement_raw = glob.glob(
                    'Data/' + date + '/' + run + '/' + str(stresses[k]) + '_' + strats[i] + '_' + enviorns[
                        j] + '/cellMovement*.gz')[0]
                boundAB_raw = glob.glob(
                    'Data/' + date + '/' + run + '/' + str(stresses[k]) + '_' + strats[i] + '_' + enviorns[
                        j] + '/boundReceptors*.gz')[0]
                growth_raw = glob.glob(
                    'Data/' + date + '/' + run + '/' + str(stresses[k]) + '_' + strats[i] + '_' + enviorns[
                        j] + '/totalCells*.gz')[0]
                # glob.glob('Data/'+ date +'/'+ run+ '/' + str(stress) + '_' + strat + '_' +enviorn +'/')

                MITrad, MI2D2D, MImove, growth = CalcData(config, bins, True, True, True, True,
                                                                         externalAB_raw, movement_raw, boundAB_raw,
                                                                         growth_raw)
                MI_tradAll[i][j].append(MITrad)
                MI2D2DAll[i][j].append(MI2D2D[0])
                MI2D1DAll[i][j].append(MImove[0])
                MI2D1DinAll[i][j].append(MImove[1])
                growthAll[i][j].append(growth)

    return MI_tradAll, MI2D2DAll, MI2D1DAll, MI2D1DinAll, growthAll

def createFigures(strats, stresses, enviorns, finalDat, run, save = True):
    date = str(utils.getTodaysDate())
    if save:
        utils.createSingleDirectory('Results', date)
    plt.clf()
    for i in range(len(enviorns)):
        plt.clf()
        for j in range(len(strats)):
            plt.plot(stresses, finalDat[2][j][i], label=strats[j] + " Receptor", marker="o")
        plt.xlabel("Stress")
        plt.ylabel("$MI(conc.;movement)$ [bits]")
        plt.grid()
        plt.legend()
        plt.title(enviorns[i])
        if save:
            utils.createSingleDirectory('Results/' + date, run)
            plt.savefig('Results/' + date + '/' + run + '/MI_move_' + enviorns[i] + '.pdf')
        else:
            plt.show()

    plt.clf()
    for i in range(len(enviorns)):
        plt.clf()
        for j in range(len(strats)):
            plt.plot(stresses, finalDat[0][j][i], label=strats[j] + " Receptor", marker="o")
        plt.xlabel("Stress")
        plt.ylabel("$MI(conc.;bound)$ [bits]")
        plt.grid()
        plt.legend()
        plt.title(enviorns[i])
        if save:
            utils.createSingleDirectory('Results/' + date, run)
            plt.savefig('Results/' + date + '/' + run + '/MI_Trad_' + enviorns[i] + '.pdf')
        else:
            plt.show()

    plt.clf()
    for i in range(len(enviorns)):
        plt.clf()
        for j in range(len(strats)):
            plt.plot(stresses, finalDat[3][j][i], label=strats[j] + " Receptor", marker="o")
        plt.xlabel("Stress")
        plt.ylabel("Input Entropy right-left [bits]")
        plt.grid()
        plt.legend()
        plt.title(enviorns[i])
        if save:
            utils.createSingleDirectory('Results/' + date, run)
            plt.savefig('Results/' + date + '/' + run + '/MI_input_' + enviorns[i] + '.pdf')
        else:
            plt.show()

    plt.clf()
    for i in range(len(enviorns)):
        plt.clf()
        for j in range(len(strats)):
            plt.plot(stresses, finalDat[1][j][i], label=strats[j] + " Receptor", marker="o")
        plt.xlabel("Stress")
        plt.ylabel("$MI(conc.;bound)$ right - left [bits]")
        plt.grid()
        plt.legend()
        plt.title(enviorns[i])
        if save:
            utils.createSingleDirectory('Results/' + date, run)
            plt.savefig('Results/' + date + '/' + run + '/MI_2D2D_' + enviorns[i] + '.pdf')
        else:
            plt.show()

    plt.clf()
    for i in range(len(enviorns)):
        plt.clf()
        for j in range(len(strats)):
            plt.plot(stresses, finalDat[4][j][i], label=strats[j] + " Receptor", marker="o")
        plt.xlabel("Stress")
        plt.ylabel("Growth")
        plt.grid()
        plt.legend()
        plt.title(enviorns[i])
        if save:
            utils.createSingleDirectory('Results/' + date, run)
            plt.savefig('Results/' + date + '/' + run + '/Growth_' + enviorns[i] + '.pdf')
        else:
            plt.show()

def createFigureMIGrowth(growths, MI_moves, MI_trad, Strategies, stresses):
    color_Scheme = ["Blue", "Red"]

    max_stress = len(stresses)
    for i in range(len(growths)):
        for j in range(len(growths[i])-1):
            plt.plot([MI_moves[i][j], MI_moves[i][j+1]], [growths[i][j], growths[i][j+1]], color = colors.findcolor(max_stress,0,color_Scheme,j ), zorder=3)
        plt.scatter(MI_moves[i], growths[i], label = Strategies[i], zorder=3)


    from scipy import interpolate
    for i in range(len(stresses)):
        y = np.array(growths)[:, i]
        x = np.array(MI_moves)[:, i]

        model4 = np.poly1d(np.polyfit(x, y, 4))
        polyline = np.linspace(0, 10, 100)
        l.log(x)
        l.log(y)
        #tck, u = interpolate.splprep([x, y], s = 0)
        #plt.plot(polyline, model4(polyline), label = str(stresses[i]))
        #xnew, ynew = interpolate.splev(np.linspace(0, 1, 100), tck, der=0)
        #plt.plot(x, y, label = str(stresses[i]))
        #plt.plot(x, y, 'o', xnew, ynew)


    plt.xlabel("MI")
    plt.ylabel("Growth")
    plt.legend()
    plt.grid(zorder=0)
    #plt.xlim([0, 2])
    plt.ylim([-2, 5])
    plt.show()

def createFigureMIGrowthStress(growths, MI_moves, MI_trad, Strategies, stresses, stress):
    j = stresses.index(stress)
    for i in range(len(growths)):
        plt.scatter(MI_moves[i][j], growths[i][j], label = Strategies[i], zorder=3)

    plt.xlabel("MI")
    plt.ylabel("Growth")
    plt.legend()
    plt.grid(zorder=0)
    plt.xlim([0, .2])
    plt.ylim([-2, 2])
    plt.show()

def createFigureGrowthMISyntactic(config, outputObjects):
    colors_plot1 = []
    colors_plot2 = []
    MImoves = []
    MIs = []
    MITrad = []
    growths = []
    MIWeightind = []
    MIWeight = []
    count = -1
    ratioA = []
    ratioAint = []
    MImovesDiff = []
    MIWeightedindDiff = []
    MIsDiff = []
    zero_counts = []
    MI_tradDiff = []
    zero_count = -1

    red_color = [1, 0, 0]
    green_color = [0, 1, 0]
    blue_color = [0, 0, 1]
    yellow_color = [1, 1, 0]
    colors_4D = [red_color, blue_color, green_color, yellow_color]
    for i in range(len(config.runStats.cellRatioAEmphasis)):
        for j in range(len(config.runStats.cellRatioAIntEmphasis)):
            zero_count += 1
            if config.runStats.cellRatioAEmphasis[i] == 0:
                zero_counts.append(zero_count)
    for i in range(len(config.runStats.cellRatioAEmphasis)):
        for j in range(len(config.runStats.cellRatioAIntEmphasis)):
            print(config.runStats.cellRatioAEmphasis[i])
            print(config.runStats.cellRatioAIntEmphasis[j])
            print(outputObjects[count].PrimitiveOutput.totalcellsFile)
            print()
            # print("here" + str(output_objects[count].RunClacOut.MI2D2D))
            count += 1
            if outputObjects[count].RunClacOut.MI2D2D != '' and not isinstance(outputObjects[count].RunClacOut.growth, str):# and config.runStats.cellRatioAEmphasis[i] >= 0 and config.runStats.cellRatioAIntEmphasis[j] >= 0:
            #if outputObjects[count].RunClacOut.MI2D2D != '' and config.runStats.cellRatioAEmphasis[i] >= 0 and config.runStats.cellRatioAIntEmphasis[j] >= 0:
            #if outputObjects[count].RunClacOut.MI2D2D != '' and config.runStats.cellRatioAEmphasis[i] <= 0 and config.runStats.cellRatioAIntEmphasis[j] >= 0:
            #if outputObjects[count].RunClacOut.MI2D2D != '' and config.runStats.cellRatioAEmphasis[i] >= 0 and config.runStats.cellRatioAIntEmphasis[j] <= 0:
            #if outputObjects[count].RunClacOut.MI2D2D != '' and config.runStats.cellRatioAEmphasis[i] <= 0 and config.runStats.cellRatioAIntEmphasis[j] <= 0:
            #if outputObjects[count].RunClacOut.MI2D2D != '' and config.runStats.cellRatioAIntEmphasis[j] <= 0:
            #if outputObjects[count].RunClacOut.MI2D2D != '' and  config.runStats.cellRatioAEmphasis[i] >= 0:
                # if c.runStats.cellRatioAEmphasis[i] != 0:
                #    colors_plot.append([1,0,0])
                # else:
                #    colors_plot.append([0,0,1])
                colors_plot1.append(colors.findColor4D(colors_4D, (config.runStats.cellRatioAEmphasis[i] + 50) / 100,
                                                       (config.runStats.cellRatioAIntEmphasis[j] + 50) / 100))
                ratioA.append(config.runStats.cellRatioAEmphasis[i])
                ratioAint.append(config.runStats.cellRatioAIntEmphasis[j])
                color_Scheme1 = ['Black', 'Red']
                color_Scheme2 = ['Black', 'Blue']
                # colors_plot1.append(colors.findcolor(50, -50, color_Scheme1, c.runStats.cellRatioAEmphasis[i]))
                # colors_plot2.append(colors.findcolor(50, -50, color_Scheme2, c.runStats.cellRatioAIntEmphasis[j]))
                #mimove1 = float(outputObjects[count].RunClacOut.MI2D1D.split(',')[0][1:].strip())
                #mimove2 = float(outputObjects[count].RunClacOut.MI2D1D.split(',')[3][1:].strip())
                mimove_count = outputObjects[count].RunClacOut.MI2D1D

                #mimovezero1 = float(outputObjects[zero_counts[j]].RunClacOut.MI2D1D.split(',')[0][1:].strip())
                #mimovezero2 = float(outputObjects[zero_counts[j]].RunClacOut.MI2D1D.split(',')[3][1:].strip())
                mimove_zero = outputObjects[zero_counts[j]].RunClacOut.MI2D1D
                # MImoves.append(output_objects[count].RunClacOut.MI2D1D)
                MImoves.append(mimove_count)
                # MImovesDiff.append(output_objects[count].RunClacOut.MI2D1D - output_objects[zero_counts[j]].RunClacOut.MI2D1D)
                MImovesDiff.append(mimove_count - mimove_zero)

                #num_zero1 = float(outputObjects[zero_counts[j]].RunClacOut.MI2D2D.split(',')[0][1:].strip())
                #num_zero2 = float(outputObjects[zero_counts[j]].RunClacOut.MI2D2D.split(',')[3][1:].strip())
                #MI2d2d_zero = outputObjects[zero_counts[j]].RunClacOut.MI2D2D
                #num1 = float(outputObjects[count].RunClacOut.MI2D2D.split(',')[0][1:].strip())
                #num2 = float(outputObjects[count].RunClacOut.MI2D2D.split(',')[3][1:].strip())
                mi2d2d_count = outputObjects[count].RunClacOut.MI2D2D
                MI2d2d_zero = outputObjects[zero_counts[j]].RunClacOut.MI2D2D
                # print(num)
                MIs.append(mi2d2d_count)
                MIsDiff.append(mi2d2d_count - MI2d2d_zero)
                MITrad.append(outputObjects[count].RunClacOut.MItrad)
                #MI_tradDiff.append((outputObjects[zero_counts[j]].RunClacOut.MItrad) - (outputObjects[count].RunClacOut.MItrad))
                growths.append(outputObjects[count].RunClacOut.growth)
                MIWeightind.append(outputObjects[count].RunClacOut.intABWeightind)
                MIWeightedindDiff.append(outputObjects[count].RunClacOut.intABWeightind - outputObjects[zero_counts[j]].RunClacOut.intABWeightind)
                MIWeight.append(outputObjects[count].RunClacOut.intweightMImove)

    """""
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the data points as scatter plot
    print(len(growths))
    print(len(ratioA))
    print(growths)
    ax.scatter(ratioA, ratioAint, growths)

    # Create a surface from the data points
    surf = ax.plot_trisurf(ratioA, ratioAint, growths, cmap='viridis', edgecolor='none')

    # Add a color bar to the plot
    fig.colorbar(surf)

    # Set labels for the axes
    ax.set_xlabel('Receptor Allocation Sigmoid Coefficient')
    ax.set_ylabel('Strategy Sigmoid Coefficient')
    ax.set_zlabel('\'Useful\' Information (Adaptive) - \'Useful\' Information (Equal)')
    ax.set_zlabel('Growth')

    # Show the plot
    plt.show()
    l.exit()
    """""
    import scipy.stats as stats
    # MIs 0.28974572341238464
    # MImoves 0.289745723412386
    # MITrad
    # MIWeightedindDiff
    x = MIWeightedindDiff
    y = growths

    x_arr = np.array(x)
    y_arr = np.array(y)

    # Calculate the correlation coefficient and p-value
    corr_coef, p_value = stats.pearsonr(x_arr, y_arr)

    print("Correlation coefficient:", corr_coef)
    print("P-value:", p_value)

    rects = colors.drawRectangles(100, colors_4D)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    ax1.scatter(x, growths, color=colors_plot1, zorder=3)
    slope, intercept = np.polyfit(x, y, 1)
    # Plot the data and the line of best fit

    newpoints = []
    for i in range(len(x)):
        newpoints.append(slope * x[i] + intercept)
    ax1.plot(x, newpoints, color='red')
    # plt.scatter(2, 2, color = 'Blue', label = "Strategy Sigmoid Coefficient")
    # plt.scatter(2, 2, color = 'Red', label = "Receptor Allocation Sigmoid Coefficient")
    ax1.set_xlabel("MI Weighted Information")
    #ax1.set_xlabel("Syntactic Information")
    # plt.xlim([0, 0.6])
    ax1.set_ylabel("Growth")
    # plt.ylim([-.1, .25])
    ax1.grid(zorder=0)
    # plt.legend()

    for box in rects:
        ax2.add_patch(box)

    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    #ax2.set_xticks([0, 1], [-50, 50])
    ax2.set_xticks([0, 1], [-50, 50])
    #ax2.set_yticks([0, 1], [-50, 50])
    ax2.set_yticks([0, 1], [-50, 50])
    ax2.set_xlabel("Receptor Allocation Gain")
    ax2.set_ylabel("Strategy Gain")
    plt.show()

