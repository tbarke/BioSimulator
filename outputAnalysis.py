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

def CalcData(config, bins, MITradFlag, MI2d2dFlag, MIMoveFlag, growthFlag, extABFile = None, moveFile = None, boundFile = None, totalCellsFile = None):
    extAB = None
    move = None
    boundAB = None
    if MITradFlag or MI2d2dFlag or MIMoveFlag:
        extAB = utils.loadDataDate(extABFile, False)
    if MIMoveFlag:
        move = utils.loadDataDate(moveFile, False)
    if MITradFlag or MI2d2dFlag:
        boundAB = utils.loadDataDate(boundFile, False)

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
    if MITradFlag or MI2d2dFlag or MIMoveFlag:
        for i in range(len(extAB[0])):
            for j in range(len(extAB[0][i])):
                count += 1
                if MIMoveFlag:
                    move_arr.append(move[0][i][j])
                if MI2d2dFlag or MIMoveFlag:
                    ext_A.append(extAB[0][i][j][1] - extAB[0][i][j][0])
                    ext_B.append(extAB[0][i][j][3] - extAB[0][i][j][2])
                if MI2d2dFlag:
                    boundA.append(boundAB[0][i][j][1] - boundAB[0][i][j][0])
                    boundB.append(boundAB[0][i][j][3] - boundAB[0][i][j][2])

                if MITradFlag:
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

    MI2D2D = None
    MITrad = None
    MImove = None
    growth = None

    k = 1
    if count > 200000:
        k = math.ceil(count / 100000)

    if MITradFlag:
        MITrad = MI_trad([extAL, extAR, extBL, extBR], [boundAL, boundAR, boundBL, boundBR], k, bins)
    if MI2d2dFlag:
        MI2D2D = utils.MI2D2D(ext_A, ext_B, boundA, boundB, bins)
    if MIMoveFlag:
        MImove = utils.MI2D1D(ext_A, ext_B, move_arr, bins)
    if growthFlag:
        growth = calcGrowth(totalCellsFile, config)

    return MITrad, MI2D2D, MImove, growth

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