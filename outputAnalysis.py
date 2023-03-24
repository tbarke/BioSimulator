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

def CalcData(config, bins, MITradFlag, MI2d2dFlag, MIMoveFlag, growthFlag, intWeightFlag, intABFile = None, extABFile = None, moveFile = None, boundFile = None, totalCellsFile = None):
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


    #------------------------
    #k = 10
    #intAextAL = reduce( intAextAL, k)
    #intAextAR = reduce(intAextAR , k)
    #intAMove = reduce(intAMove , k)

    #intBextBL = reduce(intBextBL, k)
    #intBextBR = reduce(intBextBR, k)
    #intBMove = reduce(intBMove, k)
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
    count = -1
    ratioA = []
    ratioAint = []
    MImovesDiff = []
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
            # print("here" + str(output_objects[count].RunClacOut.MI2D2D))
            count += 1
            if outputObjects[
                count].RunClacOut.MI2D2D != '':  # and c.runStats.cellRatioAEmphasis[i] >= 0 and c.runStats.cellRatioAIntEmphasis[j] >= 0:
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
                #mimove_zero = outputObjects[zero_counts[j]].RunClacOut.MI2D1D
                # MImoves.append(output_objects[count].RunClacOut.MI2D1D)
                MImoves.append(mimove_count)
                # MImovesDiff.append(output_objects[count].RunClacOut.MI2D1D - output_objects[zero_counts[j]].RunClacOut.MI2D1D)
                #MImovesDiff.append(mimove_count - mimove_zero)

                #num_zero1 = float(outputObjects[zero_counts[j]].RunClacOut.MI2D2D.split(',')[0][1:].strip())
                #num_zero2 = float(outputObjects[zero_counts[j]].RunClacOut.MI2D2D.split(',')[3][1:].strip())
                #MI2d2d_zero = outputObjects[zero_counts[j]].RunClacOut.MI2D2D
                #num1 = float(outputObjects[count].RunClacOut.MI2D2D.split(',')[0][1:].strip())
                #num2 = float(outputObjects[count].RunClacOut.MI2D2D.split(',')[3][1:].strip())
                mi2d2d_count = outputObjects[count].RunClacOut.MI2D2D
                # print(num)
                MIs.append(mi2d2d_count)
                #MIsDiff.append(mi2d2d_count - MI2d2d_zero)
                MITrad.append(outputObjects[count].RunClacOut.MItrad)
                #MI_tradDiff.append((outputObjects[zero_counts[j]].RunClacOut.MItrad) - (outputObjects[count].RunClacOut.MItrad))
                growths.append(outputObjects[count].RunClacOut.growth)

    """""
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the data points as scatter plot
    ax.scatter(ratioA, ratioAint, MImovesDiff)

    # Create a surface from the data points
    surf = ax.plot_trisurf(ratioA, ratioAint, MImovesDiff, cmap='viridis', edgecolor='none')

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

    rects = colors.drawRectangles(100, colors_4D)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    ax1.scatter(MImoves, growths, color=colors_plot1, zorder=3)
    # plt.scatter(2, 2, color = 'Blue', label = "Strategy Sigmoid Coefficient")
    # plt.scatter(2, 2, color = 'Red', label = "Receptor Allocation Sigmoid Coefficient")
    # ax1.set_xlabel("\'Useful\' Information")
    ax1.set_xlabel("Syntactic Information")
    # plt.xlim([0, 0.6])
    ax1.set_ylabel("Growth")
    # plt.ylim([-.1, .25])
    ax1.grid(zorder=0)
    # plt.legend()

    for box in rects:
        ax2.add_patch(box)

    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.set_xticks([0, 1], [-50, 50])
    ax2.set_yticks([0, 1], [-50, 50])
    ax2.set_xlabel("Receptor Allocation Gain")
    ax2.set_ylabel("Strategy Gain")
    plt.show()

