import math
import matplotlib.pyplot as plt
import random
import statistics
import colors
import numpy as np

def calculateBoundKinetic(receptors, concentration, dissocociationConstant, newRand, noise, mode):
    boundRec = 0

    if mode == "simulate":
        #print("simulate is running")
        if receptors == 0:
            return 0
        prob = concentration/(dissocociationConstant+concentration)
        return np.random.binomial(receptors, prob, 1)[0]

        #for i in range(receptors):
        #    out = random.choices([0,1], [prob, 1-prob])
        #    if out[0] == 0.0:
        #        boundRec += 1
        #return boundRec
    elif mode == "gaussian":
        prob = concentration / (dissocociationConstant + concentration)
        std = math.sqrt(receptors * prob * (1 - prob))
        std = std * noise
        mean = receptors * prob
        #print()
        #print(prob)
        #print(receptors)
        flag = True
        while(flag):
            if mean == 0:
                return 0
            sample = np.random.normal(mean, std, 1)[0]
            if sample >= 0 and sample < receptors:
                #print(round(sample))
                return round(sample)

    elif mode == "rand":
        if receptors == 0:
            return 0
        return np.random.randint(0, receptors, 1)[0]

    print("Incorrect binding type")

    return boundRec

def highest(arr):
    max = 0.0
    for i in arr:
        if i > max:
            max = i
    return max

#should be based on how much is already in the cell
def calculateAbsorbtion(concentration, rate):
    #print(rate)
    return math.floor(rate * concentration)

def createPositionHash(length):
    hash = {0: []}
    for i in range(1,length):
        hash[i] = []
    return hash

def ATPcomplex(A, B, Acost, Bcost):
    newA = math.floor(A/Acost)
    newB = math.floor(B/Bcost)
    ATPFormed = 0
    if newA < newB:
        ATPFormed = newA
    if newB < newA:
        ATPFormed = newB
    if newB == newA:
        ATPFormed = newA
    return [ATPFormed, ATPFormed*Acost, ATPFormed*Bcost]

def printPOSHASH(posHash):
    ret = ""
    for pos in posHash:
        cells = posHash[pos]
        for i in range(len(cells)):
            ret += ("c" + str(i) + ":")
        ret += ", "
    print(ret)

def findCellID(ID, list):
    for i in list:
        if i.ID == ID:
            return i
    return

def fileToArr(fileName):
    file = open(fileName, "r")
    arr = []
    for line in file:
        currentString = line.strip()
        try:
            num = float(currentString)
            arr.append(num)
        except:
            print("could not parse string")
    return arr

def plotAllCells(testEnviornment, run):
    allCells = testEnviornment.toArr()
    xAxis = []
    for i in range(len(allCells)):
        xAxis.append(i + 1)
    plt.bar(xAxis, allCells, color = ["blue"])
    plt.savefig("Cells/runCells"+str(run)+".png")
    plt.clf()

#do not run this function, it changes A and B concentration set
def plotAllConcentrat(testEnviornment, run):
    A = testEnviornment.Aconcentrations
    B = testEnviornment.Bconcentrations
    xAxis = []
    for i in range(len(B)):
        xAxis.append(i + 1)

    for i in range(len(A)):
        if i % 2 == 0:
            A[i] = 0.0
        if i %2 ==1:
            B[i] = 0.0

    plt.title("Concentration Profile")
    plt.xlabel("location")
    plt.ylabel("concentration (units)")
    plt.bar(xAxis, A, color = ["blue"])
    plt.bar(xAxis, B, color=["orange"])
    plt.savefig("Concentrations/runConcen"+str(run)+".png")
    plt.clf()

def plotReceptors(testEnviornment, run):
    AllRectp = testEnviornment.returnAverageRec()
    As = []
    Bs = []
    for i in AllRectp:
        As.append(i[0])
        Bs.append(i[1])
    xAxis = range(100)
    plt.bar(xAxis,As, color = ["red"])
    plt.savefig("Receptors/A/runReceptorsA" + str(run)+".png")
    plt.clf()
    plt.bar(xAxis, Bs, color=["orange"])
    plt.savefig("Receptors/B/runReceptorsB" + str(run) + ".png")
    plt.clf()

def plotGradients(testEnviornment):
    As = testEnviornment.allGradientsA
    Bs = testEnviornment.allGradientsB
    plt.hexbin(As, Bs)
    plt.show()

def plotVelocities(testEnviornment):
    vel = testEnviornment.cellVelocities
    xAxis = range(len(vel))
    plt.hist(vel)
    plt.show()

def countconcnetrations(arr):
    ret = 0.0
    for i in range(len(arr)):
        ret += arr[i]
    return ret

def checkzero(num):
    if num < 0.0:
        return 0.0
    return num

def outputStr(message, file):
    f = open(file, "a")
    f.write(message)
    f.close()

#outputs data in the arrays of datax, datay into file: file
def outputData(datax, datay, file):
    ret = ""
    for i in range(len(datay)):
        ret += str(datay) + ", "
    ret += "\n\n"
    for i in range(len(datax)):
        for j in range(len(datax[i])):
            ret += str(datax[i][j]) + ", "
        ret += "\n"
    outputStr(ret, file)


#for data organized in terms of runs in time (not samples in time)
def outputDataError(datax, datay, xlabel, ylabel, title, YMax, iteration, figureDIR, dataDIR, flag):
    #print(datax)
    timeLength = len(datax[0])
    finaldata = [0] * timeLength
    errtimeLength = 0
    if timeLength > 10:
        errtimeLength = 10
    else:
        errtimeLength = timeLength
    finalSTD = []
    errdatay = []
    errdata = []

    #finds final data (averaged)
    for i in range(timeLength):
        currSum = 0.0
        currlength = len(datax)
        currArr = []
        for j in range(currlength):
            currSum += datax[j][i]/currlength
            if i % math.floor(timeLength/errtimeLength) == 0:
                currArr.append(datax[j][i])

        if i % math.floor(timeLength / errtimeLength) == 0:
            finalSTD.append(statistics.stdev(currArr))
            errdatay.append(i*10)
            errdata.append(currSum)

        finaldata[i] = currSum
    #print(errdatay)
    #print(errdata)
    if flag:
        plt.errorbar(errdatay, errdata, yerr = finalSTD, fmt='.k', capsize=10, elinewidth=1)
    #print("datay")
    #print(datay)
    #print(finaldata)
    plt.plot(datay, finaldata)
    titleplotstring = title[0]
    fulltitleString = title[0]
    for i in range(1, len(title)):
        titleplotstring += "\n" + title[i]
        fulltitleString += " " + title[i]

    plt.title(titleplotstring)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid()
    #plt.ylim(0, YMax)
    plt.savefig(figureDIR+fulltitleString+iteration+".pdf")
    outputData(datax, datay, dataDIR+ fulltitleString +iteration + ".txt")
    plt.clf()
    #plt.show()

def heatMap(X,Y,Z,divisions, limit, brought_max, min, useMax, xlabel, ylabel, heatlabel,color_scheme):
    #plt.clf()
    def findCoord(x,y):
        if x>limit-1 or y > limit-1:
            return -1,-1
        mult = divisions/limit
        ret_x_init = x*mult
        ret_y_init = y*mult
        ret_x = math.floor(ret_x_init)
        ret_y = math.floor(ret_y_init)
        #if ret_x_init == ret_x and ret_x != divisions-1:
        #    ret_x += 1
        ##if ret_y_init == ret_y and ret_y != divisions-1:
        #    ret_y += 1
        return ret_x, ret_y

    def findcolor(max, min, color_scheme, num):
        portion = (num - min) / (max - min)
        index = portion * (len(color_scheme) - 1)
        index_int_bottom = math.floor(index)
        if index_int_bottom == len(color_scheme) - 1:
            col = colors.getColorCode(color_scheme[index_int_bottom])
            return col[0]/255, col[1]/255, col[2]/255
        index_int_top = index_int_bottom + 1
        index = index - index_int_bottom
        smallColor = colors.getColorCode(color_scheme[index_int_bottom])
        largetColor = colors.getColorCode(color_scheme[index_int_top])
        new_x = 0
        new_y = 0
        new_z = 0
        if smallColor[0] > largetColor[0]:
            new_x = (largetColor[0] + (smallColor[0] - largetColor[0]) * (1-index)) / 255
        else:
            new_x = (smallColor[0] + (math.fabs(largetColor[0] - smallColor[0]) * index)) / 255

        if smallColor[1] > largetColor[1]:
            new_y = (largetColor[1] + (smallColor[1] - largetColor[1]) * (1-index)) / 255
        else:
            new_y = (smallColor[1] + (math.fabs(largetColor[1] - smallColor[1]) * index)) / 255

        if smallColor[2] > largetColor[2]:
            new_z = (largetColor[2] + (smallColor[2] - largetColor[2]) * (1-index)) / 255
        else:
            new_z = (smallColor[2] + (math.fabs(largetColor[2] - smallColor[2]) * index)) / 255
        return [new_x, new_y, new_z]

    inter_XYZ = []
    final_XYZ = []
    for i in range(divisions):
        curr_arr = []
        curr_arr_final = []
        for k in range(divisions):
            curr_arr.append([])
            curr_arr_final.append(0.0)
        inter_XYZ.append(curr_arr)
        final_XYZ.append(curr_arr_final)

    index_offset = limit/divisions
    total = 0.0
    for i in range(len(X)):
        #print(X[i], Y[i])
        new_x, new_y = findCoord(X[i],Y[i])
        #print(new_x, new_y)
        #print()
        if new_x != -1:
            inter_XYZ[new_x][new_y].append(Z[i])
        total += Z[i]


    max = 0.0
    for i in range(len(inter_XYZ)):
        for j in range(len(inter_XYZ[i])):
            sum = 0.0
            for k in range(len(inter_XYZ[i][j])):
                sum += inter_XYZ[i][j][k]
            sum = sum/total
            if max < sum:
                max = sum
            #print(sum)
            final_XYZ[i][j] = sum

    if useMax:
        max = brought_max

    #print(max)
    fig, ax = plt.subplots(1, 2, gridspec_kw={'width_ratios': [30, 1]})

    fig.tight_layout()
    ax[1].yaxis.tick_right()
    ax[1].get_xaxis().set_visible(False)
    x = [0, 1]
    length = 200
    #color_scheme = ['Purple', 'Blue', 'Green', 'Yellow']
    #color_scheme = ['p1', 'g1', 'y1']
    for i in range(length):
        y = [1 + i, 1 + i]
        y2 = [i, i]
        col = findcolor(length, 0, color_scheme, i)
        #print(col)
        ax[1].fill_between(x, y, y2, facecolor=col)
    ax[1].set_xlim(0.5, 0.75)
    ax[1].set_ylim(1, length - 1)

    #print(len(final_XYZ))
    for i in range(len(final_XYZ)):
        for j in range(len(final_XYZ[i])):
            index_i = i*index_offset
            index_j = j*index_offset
            x = [index_i, index_i + index_offset]
            y = [index_j + index_offset, index_j + index_offset]
            y2 = [index_j, index_j]
            #print(final_XYZ[i][j])
            col = findcolor(max,min,color_scheme, final_XYZ[i][j])
            #print(col)
            ax[0].fill_between(x, y, y2, facecolor = col)

    ax[0].axis('square')
    plt.subplots_adjust(wspace=0)
    tick_divisions = 5
    tick_in = max/tick_divisions
    ticks = []
    tick_labels = []
    for i in range(tick_divisions-1):
        ticks.append(((tick_in*(i+1))/max)*200)
        lab = '{:.3}'.format(tick_in*(i+1))
        tick_labels.append(lab)
    #labels = ['first', 'second', 'third']
    ax[1].set_yticks(ticks)
    ax[1].set_yticklabels(tick_labels)
    ax[0].set_xlim(0, limit)
    ax[0].set_ylim(0, limit)
    ax[0].set_xlabel(xlabel)
    ax[0].set_ylabel(ylabel)
    ax[1].set_ylabel(heatlabel, rotation = 270)
    ax[1].yaxis.set_label_coords(6, 0.5)
    ax[1].yaxis.set_label_position("right")
    x1, y1 = [25, 100], [25, 25]
    x2, y2 = [25, 25], [25, 100]
    ax[0].plot(x1, y1, x2, y2, color="black", lw=2, alpha = 0.6)
    plt.tight_layout()
    plt.show()
    return max

