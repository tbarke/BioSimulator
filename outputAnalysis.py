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

def MIMoveWeighted(concs, moves, inters, binsMI, binsintAB):
    # list of concs: AR, AL, BR, BL
    # list of moves: absolute velocity
    # list of internal states intA, intB
    # find max intA and intA

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

    # divide into total bins used
    binedgeA, binedgeB, binAB = findBins(inters[0], inters[1])
    #print(binedgeA)
    #print(binedgeB)

    def findBin(intA, intB):
        for i in range(len(binedgeA)):
            if intA >= binedgeA[i] and intA < binedgeA[i + 1]:
                for j in range(len(binedgeA)):
                    if intB >= binedgeB[j] and intB < binedgeB[j +1]:
                        return i, j

    for i in range(len(inters[0])):
        curr_intA = inters[0][i]
        curr_intB = inters[1][i]
        indexA, indexB = findBin(curr_intA, curr_intB)
        binAB[indexA][indexB].append(i)

    binProbs = utils.list_empty([binsintAB, binsintAB], 0)
    for i in range(len(binAB)):
        for j in range(len(binAB[0])):
            binProbs[i][j] = len(binAB[i][j])/len(moves[0])

    MIweight = 0
    binABind = utils.list_empty([binsintAB, binsintAB], [])
    for i in range(len(binAB)):
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

def CalcData(config, bins, MITradFlag, MI2d2dFlag, MIMoveFlag, growthFlag, intWeightFlag, ABintdynamicWeightFlag, intABFile = None, extABFile = None, moveFile = None, boundFile = None, totalCellsFile = None):
    extAB = None
    move = None
    boundAB = None
    intAB = None
    if MITradFlag or MI2d2dFlag or MIMoveFlag or intWeightFlag:
        extAB = utils.loadDataDate(extABFile, False)
    if MIMoveFlag or intWeightFlag:
        move = utils.loadDataDate(moveFile, False)
    if MITradFlag or MI2d2dFlag or MIMoveFlag:
        boundAB = utils.loadDataDate(boundFile, False)

    if intWeightFlag:
        intAB = utils.loadDataDate(intABFile, False)

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

    if MITradFlag or MI2d2dFlag or MIMoveFlag or intWeightFlag:
        for i in range(len(extAB[0])):
            for j in range(len(extAB[0][i])):
                count += 1
                if MIMoveFlag:
                    move_arr.append(move[0][i][j])
                if MI2d2dFlag or MIMoveFlag:
                    ext_A.append(extAB[0][i][j][1] - extAB[0][i][j][0])
                    ext_B.append(extAB[0][i][j][3] - extAB[0][i][j][2])
                if MI2d2dFlag:
                    boundAB_all = math.fabs(boundAB[0][i][j][1] - boundAB[0][i][j][0]) + math.fabs(boundAB[0][i][j][3] - boundAB[0][i][j][2])
                    if boundAB_all > 0:
                        boundA.append((boundAB[0][i][j][1] - boundAB[0][i][j][0])/boundAB_all)
                        boundB.append((boundAB[0][i][j][3] - boundAB[0][i][j][2])/boundAB_all)
                    else:
                        boundA.append(0)
                        boundB.append(0)

                if MITradFlag or MIMoveFlag or MI2d2dFlag:
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
                if intWeightFlag:
                    intA_curr = intAB[0][i][j][0]
                    intB_curr = intAB[0][i][j][1]
                    intAall.append(intA_curr)
                    intBall.append(intB_curr)
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
    growth = None

    k = 1
    if count > 20000:
        k = math.ceil(count / 10000)

    def reduce(arr, k):
        new = []
        for i in range(0, len(arr), k):
            new.append(arr[i])
        return new

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


    if MITradFlag:
        MITrad = MI_trad([extAL, extAR, extBL, extBR], [boundAL, boundAR, boundBL, boundBR], k, bins)
    if MI2d2dFlag:

        MI2D2D1 = utils.MI2D2D(extAR, extAL, boundAR, boundAL, bins)[0]
        MI2D2D2 = utils.MI2D2D(extBR, extBL, boundBR, boundBL, bins)[0]
        MI2D2D = MI2D2D1 + MI2D2D2
    if MIMoveFlag:
        MImove1 = utils.MI2D1D(extAR, extAL, move_arr, bins)[0]
        MImove2 = utils.MI2D1D(extBR, extBL, move_arr, bins)[0]
        MImove = MImove1 + MImove2
    if growthFlag:
        growth = calcGrowth(totalCellsFile, config)
    if intWeightFlag:
        intWeightAMI = utils.MI2D1D(intAextAL, intAextAR, intAMove, bins)[0]
        intWeightBMI = utils.MI2D1D(intBextBL, intBextBR, intBMove, bins)[0]
        intweightMIMove = intWeightAMI + intWeightBMI
    if ABintdynamicWeightFlag:
        #concs, moves, inters, binsMI, binsintAB
        concs = [extAR, extAL, extBR, extBL]
        moves = [move_arr]
        inters = [intAall, intBall]
        weightedMI, binMI = MIMoveWeighted(concs, moves, inters, 30, 20)


    return MITrad, MI2D2D, MImove, growth, intweightMIMove

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
            count += 1
            # print("here" + str(output_objects[count].RunClacOut.MI2D2D))

            if outputObjects[count].RunClacOut.MI2D2D != '':# and config.runStats.cellRatioAEmphasis[i] >= 0 and config.runStats.cellRatioAIntEmphasis[j] >= 0:
            #if outputObjects[count].RunClacOut.MI2D2D != '' and config.runStats.cellRatioAEmphasis[i] >= 0 and config.runStats.cellRatioAIntEmphasis[j] >= 0:
            #if outputObjects[count].RunClacOut.MI2D2D != '' and config.runStats.cellRatioAEmphasis[i] <= 0 and config.runStats.cellRatioAIntEmphasis[j] >= 0:
            #if outputObjects[count].RunClacOut.MI2D2D != '' and config.runStats.cellRatioAEmphasis[i] >= 0 and config.runStats.cellRatioAIntEmphasis[j] <= 0:
            #if outputObjects[count].RunClacOut.MI2D2D != '' and config.runStats.cellRatioAEmphasis[i] <= 0 and config.runStats.cellRatioAIntEmphasis[j] <= 0:
            #if outputObjects[count].RunClacOut.MI2D2D != '' and config.runStats.cellRatioAIntEmphasis[j] <= 0:
            #if outputObjects[count].RunClacOut.MI2D2D != '' and config.runStats.cellRatioAEmphasis[i] >= 0:
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
                MIWeightind.append(outputObjects[count].RunClacOut.intabweightind)
                MIWeightedindDiff.append(outputObjects[count].RunClacOut.intabweightind - outputObjects[zero_counts[j]].RunClacOut.intabweightind)
                MIWeight.append(outputObjects[count].RunClacOut.intweightMImove)

    """""
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the data points as scatter plot
    ax.scatter(ratioA, ratioAint, growths)

    # Create a surface from the data points
    surf = ax.plot_trisurf(ratioA, ratioAint, growths, cmap='viridis', edgecolor='none')

    # Add a color bar to the plot
    fig.colorbar(surf)

    # Set labels for the axes
    ax.set_xlabel('Receptor Allocation Sigmoid Coefficient')
    ax.set_ylabel('Strategy Sigmoid Coefficient')
    ax.set_zlabel('\'Useful\' Information (Adaptive) - \'Useful\' Information (Equal)')

    # Show the plot
    plt.show()
    l.exit()
    """""
    """""
    path = "/Users/tyler/Library/CloudStorage/OneDrive-UniversityofNebraska-Lincoln/Biosim/Data/" + outputObjects[0].PrimitiveOutput.enviornmentConcsFile
    concentrations = utils.loadDataDate(path, False)
    print(len(concentrations))
    print(len(concentrations[0]))
    print(len(concentrations[0][0]))
    A = concentrations[0][0][0]
    B = concentrations[0][0][1]

    utils.saveData("Results/NanoCom_Data", A, "concsA_smallDis.bin")
    utils.saveData("Results/NanoCom_Data", B, "concsB_smallDis.bin")
    exit()
    print("#growths")
    print(growths)
    utils.saveData("Results/NanoCom_Data", growths, "growths_noSub.bin")
    print("#MI 2D2D")
    print(MIs)
    utils.saveData("Results/NanoCom_Data", MIs, "MIs_noSub.bin")
    print("#weightedMI")
    print(MIWeightedindDiff)
    utils.saveData("Results/NanoCom_Data", MIWeightedindDiff, "MIweighted_noSub.bin")
    print("#colors")
    print(colors_plot1)
    utils.saveData("Results/NanoCom_Data", colors_plot1, "colors_noSub.bin")
    """""
    growths_largeDis = utils.loadData("Results/NanoCom_Data/growths_largeDis.bin")
    growths_smallDis = utils.loadData("Results/NanoCom_Data/growths_smallDis.bin")
    growths_noSub = utils.loadData("Results/NanoCom_Data/growths_noSub.bin")

    MIs_largeDis = utils.loadData("Results/NanoCom_Data/MIs_largeDis.bin")
    MIs_smallDis = utils.loadData("Results/NanoCom_Data/MIs_smallDis.bin")
    MIs_noSub = utils.loadData("Results/NanoCom_Data/MIs_noSub.bin")

    weightedMI_largeDis = utils.loadData("Results/NanoCom_Data/MIweighted_largeDis.bin")
    weightedMI_smallDis = utils.loadData("Results/NanoCom_Data/MIweighted_smallDis.bin")
    weightedMI_noSub = utils.loadData("Results/NanoCom_Data/MIweighted_noSub.bin")

    colors_largeDis = utils.loadData("Results/NanoCom_Data/colors_largeDis.bin")
    colors_smallDis = utils.loadData("Results/NanoCom_Data/colors_smallDis.bin")
    colors_noSub = utils.loadData("Results/NanoCom_Data/colors_noSub.bin")

    Aconc_largeDis = utils.loadData("Results/NanoCom_Data/concsA_largeDis.bin")
    Aconc_smallDis = utils.loadData("Results/NanoCom_Data/concsA_smallDis.bin")
    Aconc_noSub = utils.loadData("Results/NanoCom_Data/concsA_noSub.bin")

    Bconc_largeDis = utils.loadData("Results/NanoCom_Data/concsB_largeDis.bin")
    Bconc_smallDis = utils.loadData("Results/NanoCom_Data/concsB_smallDis.bin")
    Bconc_noSub = utils.loadData("Results/NanoCom_Data/concsB_noSub.bin")

    import matplotlib.patches as mpatches
    import matplotlib as mpl
    mpl.rcParams['font.family'] = 'Times New Roman'
    mpl.rcParams['pdf.fonttype'] = 42
    figure, axis = plt.subplots(4, 2, figsize=(13, 20))
    x = range(1000)
    fontsize1 = 18
    fontsize2 = 14


    axis[0, 0].bar(x, Aconc_smallDis, width = 1, label = "A concentration", alpha = 0.5, zorder=3)
    axis[0, 0].bar(x, Bconc_smallDis, width = 1, label = "B concentration", alpha = 0.5, zorder=3)
    xtick = []
    xticklabels = []
    for i in range(0, 1100, 100):
        xtick.append(i)
        xticklabels.append(round(i/10))
    axis[0,0].set_xticks(xtick, xticklabels, fontsize=fontsize2)
    axis[0, 0].tick_params(axis='y', labelsize=fontsize2)
    axis[0, 0].set_xlim([0, 1000])
    axis[0, 0].set_ylim([0, 110])
    axis[0, 0].set_xlabel("Location ($\\bar{x}$)", fontsize=fontsize1)
    axis[0, 0].set_ylabel("Concentration", fontsize=fontsize1)
    axis[0, 0].set_title('a.', loc='left', weight='bold', fontsize=fontsize1)
    axis[0, 0].legend(fontsize = 12)
    axis[0,0].grid(zorder=0)


    axis[0, 1].bar(x, Aconc_noSub, width = 1, label = "A concentration", alpha= 0.5, zorder=3)
    axis[0, 1].bar(x, Bconc_noSub, width = 1, label = "B concentration", alpha= 0.5, zorder=3)
    axis[0, 1].set_xticks(xtick, xticklabels, fontsize = fontsize2)
    axis[0, 1].tick_params(axis='y', labelsize=fontsize2)
    axis[0, 1].set_xlim([0, 1000])
    axis[0, 1].set_ylim([0, 110])
    axis[0, 1].set_xlabel("Location ($\\bar{x}$)", fontsize=fontsize1)
    axis[0, 1].set_ylabel("Concentration", fontsize=fontsize1)
    axis[0, 1].set_title('b.', loc='left', weight='bold', fontsize=fontsize1)
    axis[0,1].legend(fontsize = 12)
    axis[0, 1].grid(zorder=0)

    #-----------------------------------------------------------------------------

    # MI
    x = MIs_smallDis
    y = growths_smallDis

    x_arr = np.array(x)
    y_arr = np.array(y)
    import scipy.stats as stats

    corr_coef, p_value = stats.pearsonr(x_arr, y_arr)

    print("Correlation coefficient:", corr_coef)
    print("P-value:", p_value)

    axis[1, 0].scatter(x, y, color=colors_smallDis, zorder=3)
    slope, intercept = np.polyfit(x, y, 1)
    newpoints = []
    for i in range(len(x)):
        newpoints.append(slope * x[i] + intercept)

    axis[1, 0].plot(x, newpoints, color='red', zorder=3)
    axis[1, 0].set_xlabel("Syntactic Information ($I_{synt}$) [bits]", fontsize=fontsize1)
    axis[1, 0].set_ylabel("Growth", fontsize=fontsize1)
    axis[1, 0].tick_params(axis='x', labelsize=fontsize2)
    axis[1, 0].tick_params(axis='y', labelsize=fontsize2)
    axis[1, 0].grid(zorder=0)
    axis[1, 0].text(3.4, .275, "$R = $" + str(round(corr_coef * 100) / 100), fontsize=14,
                    bbox=dict(facecolor='red', alpha=0.5))
    axis[1, 0].text(8.4, .275, "$K = $" + str(2.0), fontsize=14,
                    bbox=dict(facecolor='red', alpha=0.5))
    axis[1, 0].set_title('c.', loc='left', weight='bold', fontsize=fontsize1)
    pos = axis[1, 0].get_position()
    new_pos = [pos.x0, pos.y0 + - 0.01, pos.width, pos.height]
    axis[1, 0].set_position(new_pos)

    # MI
    x = MIs_noSub
    y = growths_noSub

    x_arr = np.array(x)
    y_arr = np.array(y)
    import scipy.stats as stats

    corr_coef, p_value = stats.pearsonr(x_arr, y_arr)

    print("Correlation coefficient:", corr_coef)
    print("P-value:", p_value)

    axis[1, 1].scatter(x, y, color=colors_noSub, zorder=3)
    slope, intercept = np.polyfit(x, y, 1)
    newpoints = []
    for i in range(len(x)):
        newpoints.append(slope * x[i] + intercept)

    axis[1, 1].plot(x, newpoints, color='red', zorder=3)
    axis[1, 1].set_xlabel("Syntactic Information ($I_{synt}$) [bits]", fontsize=fontsize1)
    axis[1, 1].set_ylabel("Growth", fontsize=fontsize1)
    axis[1, 1].tick_params(axis='x', labelsize=fontsize2)
    axis[1, 1].tick_params(axis='y', labelsize=fontsize2)
    axis[1, 1].grid(zorder=0)
    axis[1, 1].text(1.407, .895, "$R = $" + str(round(corr_coef * 100) / 100), fontsize=14,
                    bbox=dict(facecolor='red', alpha=0.5))
    axis[1, 1].text(1.550, .895, "$K = $" + str(2.0), fontsize=14,
                    bbox=dict(facecolor='red', alpha=0.5))
    axis[1, 1].set_title('d.', loc='left', weight='bold', fontsize=fontsize1)
    pos = axis[1, 1].get_position()
    new_pos = [pos.x0, pos.y0 + - 0.01, pos.width, pos.height]
    axis[1, 1].set_position(new_pos)




    # MI
    x = weightedMI_smallDis
    y = growths_smallDis

    x_arr = np.array(x)
    y_arr = np.array(y)
    import scipy.stats as stats

    corr_coef, p_value = stats.pearsonr(x_arr, y_arr)

    print("Correlation coefficient:", corr_coef)
    print("P-value:", p_value)

    axis[2, 0].scatter(x, y, color=colors_smallDis, zorder=3)
    slope, intercept = np.polyfit(x, y, 1)
    newpoints = []
    for i in range(len(x)):
        newpoints.append(slope * x[i] + intercept)

    axis[2, 0].plot(x, newpoints, color='red', zorder=3)
    axis[2, 0].set_xlabel("Subjective Information ($I_{subj}$) [bits]", fontsize=fontsize1)
    axis[2, 0].set_ylabel("Growth", fontsize=fontsize1)
    axis[2, 0].tick_params(axis='x', labelsize=fontsize2)
    axis[2, 0].tick_params(axis='y', labelsize=fontsize2)
    axis[2, 0].grid(zorder=0)
    axis[2,0].text(-.05,.275,"$R = $" + str(round(corr_coef*100)/100), fontsize = 14 ,bbox = dict(facecolor = 'red', alpha = 0.5))
    axis[2, 0].text(0.8, .275, "$K = $" + str(2.0), fontsize=14,
                    bbox=dict(facecolor='red', alpha=0.5))
    axis[2, 0].set_title('e.', loc='left', weight='bold', fontsize=fontsize1)
    axis[2, 0].set_xlim([-0.15, 1.1])
    pos = axis[2, 0].get_position()
    new_pos = [pos.x0, pos.y0 + - (0.01*2), pos.width, pos.height]
    axis[2, 0].set_position(new_pos)
    #pos = axis[1, 0].get_position()
    #new_pos = [pos.x0, pos.y0 + .05, pos.width, pos.height]
    #axis[1, 0].set_position(new_pos)
    #axis[1, 0].set_title("$R = $" + str(round(corr_coef*100)/100))

    # MIWeightedindDiff
    x = weightedMI_noSub
    y = growths_noSub

    newx = []
    newy = []
    newcolors = []
    count = 0
    for i in range(len(x)):
        if x[i] > -0.5:
            newx.append(x[i])
            newy.append(y[i])
            newcolors.append(colors_noSub[i])
        else:
            count +=1

    print("count: " + str(count) + " ----------------------------------------------------------")
    x = newx
    y = newy
    colors_noSub = newcolors

    x_arr = np.array(x)
    y_arr = np.array(y)
    import scipy.stats as stats

    corr_coef, p_value = stats.pearsonr(x_arr, y_arr)

    print("Correlation coefficient:", corr_coef)
    print("P-value:", p_value)

    axis[2, 1].scatter(x, y, color=colors_noSub, zorder=3)
    slope, intercept = np.polyfit(x, y, 1)
    newpoints = []
    for i in range(len(x)):
        newpoints.append(slope * x[i] + intercept)
    axis[2, 1].plot(x, newpoints, color='red', zorder=3)
    axis[2, 1].set_xlabel("Subjective Information ($I_{subj}$) [bits]", fontsize=fontsize1)
    axis[2, 1].set_ylabel("Growth", fontsize=fontsize1)
    axis[2, 1].tick_params(axis='x', labelsize=fontsize2)
    axis[2, 1].tick_params(axis='y', labelsize=fontsize2)
    axis[2, 1].grid(zorder=0)
    axis[2, 1].text(-0.095, .895, "$R = $" + str(round(corr_coef * 100) / 100), fontsize=14,
                    bbox=dict(facecolor='red', alpha=0.5))
    axis[2, 1].text(0.050, .895, "$K = $" + str(2.0), fontsize=14,
                    bbox=dict(facecolor='red', alpha=0.5))
    #pos = axis[1, 1].get_position()
    #new_pos = [pos.x0, pos.y0 + .05 - 0.05, pos.width, pos.height]
    #axis[1, 1].set_position(new_pos)
    axis[2, 1].set_title('f.', loc='left',weight='bold', fontsize=fontsize1)
    pos = axis[2, 1].get_position()
    new_pos = [pos.x0, pos.y0 + - (0.01 * 2), pos.width, pos.height]
    axis[2, 1].set_position(new_pos)
    #----------------

    # MIWeightedindDiff
    x = weightedMI_largeDis
    y = growths_largeDis

    x_arr = np.array(x)
    y_arr = np.array(y)
    import scipy.stats as stats

    corr_coef, p_value = stats.pearsonr(x_arr, y_arr)

    print("Correlation coefficient:", corr_coef)
    print("P-value:", p_value)

    axis[3, 0].scatter(x, y, color=colors_largeDis, zorder=3)
    slope, intercept = np.polyfit(x, y, 1)
    newpoints = []
    for i in range(len(x)):
        newpoints.append(slope * x[i] + intercept)
    axis[3, 0].plot(x, newpoints, color='red', zorder=3)
    axis[3, 0].set_xlabel("Subjective Information ($I_{subj}$) [bits]", fontsize=fontsize1)
    axis[3, 0].set_ylabel("Growth", fontsize=fontsize1)
    axis[3, 0].tick_params(axis='x', labelsize=fontsize2)
    axis[3, 0].tick_params(axis='y', labelsize=fontsize2)
    axis[3, 0].grid(zorder=0)
    axis[3, 0].text(-.05, .245, "$R = $" + str(round(corr_coef * 100) / 100), fontsize=14,
                    bbox=dict(facecolor='red', alpha=0.5))
    axis[3, 0].text(0.8, .245, "$K = $" + str(10.0), fontsize=14,
                    bbox=dict(facecolor='red', alpha=0.5))
    # pos = axis[1, 1].get_position()
    # new_pos = [pos.x0, pos.y0 + .05 - 0.05, pos.width, pos.height]
    # axis[1, 1].set_position(new_pos)
    axis[3, 0].set_title('g.', loc='left', weight='bold', fontsize=fontsize1)
    axis[3, 0].set_xlim([-0.15, 1.1])
    pos = axis[3, 0].get_position()
    new_pos = [pos.x0, pos.y0 + - (0.01 * 3), pos.width, pos.height]
    axis[3, 0].set_position(new_pos)

    #fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    rects = colors.drawRectangles(100, colors_4D)
    for box in rects:
        axis[3, 1].add_patch(box)

    axis[3, 1].set_xlim(0, 1)
    axis[3, 1].set_ylim(0, 1)
    axis[3, 1].set_xticks([0, 0.2, 0.4, 0.6, 0.8,  1], [-50, -30, -10, 10, 30, 50], fontsize = fontsize2)
    # ax2.set_xticks([0, 1], [-50, 50])
    axis[3, 1].set_yticks([0, 0.2, 0.4, 0.6, 0.8,  1],  [-50, -30, -10, 10, 30, 50], fontsize = fontsize2)
    axis[3, 1].grid(color = 'black', alpha = 0.7, linestyle = '--', linewidth = 1)
    # ax2.set_yticks([0, 1], [-50, 50])


    axis[3, 1].set_xlabel("$Gain_{Acq}$", fontsize=fontsize1)
    axis[3, 1].set_ylabel("$Gain_{Proc}$", fontsize=fontsize1)
    axis[3, 1].set_title('h.', loc='left', weight='bold', fontsize=fontsize1)
    pos = axis[3, 1].get_position()
    new_pos = [pos.x0, pos.y0 + - (0.01 * 3), pos.width, pos.height]
    axis[3, 1].set_position(new_pos)

    #plt.show()

    plt.savefig("Results/full.pdf")

    #
    #fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Plot the data and the line of best fit


    # plt.scatter(2, 2, color = 'Blue', label = "Strategy Sigmoid Coefficient")
    # plt.scatter(2, 2, color = 'Red', label = "Receptor Allocation Sigmoid Coefficient")
    #ax1.set_xlabel("Subjective Information")

    # plt.legend()
    exit()
    plt.clf()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    rects = colors.drawRectangles(100, colors_4D)
    for box in rects:
        ax2.add_patch(box)

    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.set_xticks([0, 1], [-50, 50])
    #ax2.set_xticks([0, 1], [-50, 50])
    ax2.set_yticks([0, 1], [-50, 50])
    #ax2.set_yticks([0, 1], [-50, 50])
    ax2.set_xlabel("Receptor Allocation Gain")
    ax2.set_ylabel("Strategy Gain")

    plt.show()

