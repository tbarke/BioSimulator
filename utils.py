import math
import matplotlib.pyplot as plt
import random
import statistics
import colors
import os
import pickle
import datetime
from datetime import date
import gzip
import numpy as np
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
import io
from PIL import Image
import imageio
import log
l = log.log()

def calculateBoundKinetic(receptors, concentration, dissocociationConstant, newRand, noise, mode):
    boundRec = 0
    #print(mode)
    #return concentration

    if mode == "simulate":
        #print("simulate is running")
        if receptors == 0:
            return 0
        #print(concentration)
        #print(dissocociationConstant)
        prob = concentration/(dissocociationConstant+concentration)
        #print(prob)
        #print()
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
        #print("------------------------------------H--------------------------------")
        flag = True
        while(flag):
            if mean == 0:
                return 0
            sample = np.random.normal(mean, std, 1)[0]
            #print(receptors)
            # print(sample)
            if sample >= 0 and sample < receptors:
                #print(round(sample))
                return round(sample)

    elif mode == "rand":
        if receptors == 0:
            return 0
        return np.random.randint(0, receptors, 1)[0]

    l.log("Incorrect binding type")

    return boundRec

def highest(arr):
    max = 0.0
    for i in arr:
        if i > max:
            max = i
    return max

#should be based on how much is already in the cell
def calculateAbsorbtion(concentration, rate):
    return rate * concentration

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
    l.log(ret)

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
            l.log("could not parse string")
    return arr

def plotAllCells(testEnviornment, run):
    allCells = testEnviornment.toArr()
    xAxis = []
    for i in range(len(allCells)):
        xAxis.append(i + 1)
    plt.bar(xAxis, allCells, color = ["blue"])
    l.log("Cells/runCells"+str(run)+".png")
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

def threeDmap(inputX, inputY, Z):
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator
    import numpy as np
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    # Make data.
    Xs = np.array(inputX)
    Ys = np.array(inputY)
    Xs, Ys = np.meshgrid(Xs,Ys)
    Z1 = []
    for i in range(len(Xs)):
        newZ = []
        for j in range(len(Ys)):
            newZ.append(Z[i][j])
        Z1.append(newZ)
    Zs = np.array(Z1)

    surf = ax.plot_surface(Xs, Ys, Zs, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    ax.zaxis.set_major_locator(LinearLocator(10))
    # A StrMethodFormatter is used automatically
    ax.zaxis.set_major_formatter('{x:.02f}')
    #ax.set_zlim(3, 8)
    ax.set_ylabel("Noise (as Portion of Binomal Variance)")
    ax.set_xlabel("Cell Stress")
    ax.set_zlabel("Growth")


    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

def MI2D2D(X1, X2, Y1, Y2, bins):
    def matchBinPos(binPos, y):
        for i in range(len(binPos) - 1):
            if y >= binPos[i] and y < binPos[i + 1]:
                return i

    def list_empty(shapeArr, instaniate=None):
        arr = []
        for i in range(shapeArr[0]):
            obj = instaniate.copy()
            arr.append(obj)
        if len(shapeArr) == 1:
            return arr
        newShapeArr = shapeArr[1:len(shapeArr)]
        for i, a in enumerate(arr):
            arr[i] = list_empty(newShapeArr, instaniate)
        return arr

    def entropy1d(data, edges):
        # sum over bins -p(bin)*log2(p(bin)))
        prob_arr = np.zeros(len(edges) - 1)
        for dat in data:
            prob_arr[matchBinPos(edges, dat)] += 1
        prob_arr /= len(data)
        sum1 = 0.0
        for prob in prob_arr:
            if prob > 0:
                sum1 += prob * math.log(prob, 2)
        sum1 = sum1 * -1
        return sum1

    def createBinPos(data, bins):
        max_dat = max(data)
        min_dat = min(data)
        epsilon = math.fabs(0.01 * ((max_dat - min_dat) / bins))
        return np.linspace(min_dat - epsilon, max_dat + epsilon, bins + 1)

    def entropy2d(data1, data2, edges1, edges2):
        prob_arr1 = np.zeros(len(edges1) - 1)
        prob_arr2 = np.zeros(len(edges2) - 1)
        for dat in data1:
            prob_arr1[matchBinPos(edges1, dat)] += 1
        for dat in data2:
            prob_arr2[matchBinPos(edges2, dat)] += 1
        if len(data1) > 0:
            prob_arr1 /= len(data1)
        if len(data2) > 0:
            prob_arr2 /= len(data2)
        sum1 = 0.0
        for x in prob_arr1:
            for y in prob_arr2:
                if x * y > 0:
                    sum1 += x * y * math.log(x * y, 2)
        sum1 = sum1 * -1
        return sum1

    dat_arr = list_empty([bins, bins,2], [])
    binPosY1 = createBinPos(Y1, bins)
    binPosY2 = createBinPos(Y2, bins)

    for i, y in enumerate(Y1):
        i_index = matchBinPos(binPosY1, y)
        j_index = matchBinPos(binPosY2, Y2[i])
        dat_arr[i_index][j_index][0].append(X1[i])
        dat_arr[i_index][j_index][1].append(X2[i])

    X1Edges = createBinPos(X1, bins)
    X2Edges = createBinPos(X2, bins)
    HX = entropy2d(np.asarray(X1), np.asarray(X2), X1Edges, X2Edges)
    HXgivenY = 0.0
    Y_prob, Yedges1, Yedges2 = np.histogram2d(Y1, Y2, bins, density = True)
    binwidthx = Yedges1[1] - Yedges1[0]
    binwidthy = Yedges2[1] - Yedges2[0]
    Y_prob = Y_prob*binwidthx*binwidthy
    for i, dat_arr_x in enumerate(dat_arr):
        for j, dat_arrxy in enumerate(dat_arr_x):
            HXgivenY += Y_prob[i][j]*entropy2d(np.asarray(dat_arrxy[0]), np.asarray(dat_arrxy[1]), X1Edges, X2Edges)

    return HX-HXgivenY, HX, HXgivenY


def list_empty(shapeArr, instaniate=None):
    primitive = (int, str, bool, float)

    def is_primitive(thing):
        return isinstance(thing, primitive)
    arr = []
    for i in range(shapeArr[0]):
        obj = instaniate
        if not is_primitive(instaniate):
            obj = instaniate.copy()
        arr.append(obj)
    if len(shapeArr) == 1:
        return arr
    newShapeArr = shapeArr[1:len(shapeArr)]
    for i, a in enumerate(arr):
        arr[i] = list_empty(newShapeArr, instaniate)
    return arr


def MI2D1D(X1,X2,Y,bins):
    def matchBinPos(binPos, y):
        for i in range(len(binPos) - 1):
            if y >= binPos[i] and y < binPos[i + 1]:
                return i

    def list_empty(shapeArr, instaniate=None):
        arr = []
        for i in range(shapeArr[0]):
            obj = instaniate.copy()
            arr.append(obj)
        if len(shapeArr) == 1:
            return arr
        newShapeArr = shapeArr[1:len(shapeArr)]
        for i, a in enumerate(arr):
            arr[i] = list_empty(newShapeArr, instaniate)
        return arr

    def entropy1d(data, edges):
        # sum over bins -p(bin)*log2(p(bin)))
        prob_arr = np.zeros(len(edges) - 1)
        for dat in data:
            prob_arr[matchBinPos(edges, dat)] += 1
        prob_arr /= len(data)
        sum1 = 0.0
        for prob in prob_arr:
            if prob > 0:
                sum1 += prob * math.log(prob, 2)
        sum1 = sum1 * -1
        return sum1

    def createBinPos(data, bins):
        max_dat = max(data)
        min_dat = min(data)
        epsilon = math.fabs(0.01 * ((max_dat - min_dat) / bins))
        return np.linspace(min_dat - epsilon, max_dat + epsilon, bins + 1)

    def entropy2d(data1, data2, edges1, edges2):
        prob_arr1 = np.zeros(len(edges1) - 1)
        prob_arr2 = np.zeros(len(edges2) - 1)
        for dat in data1:
            prob_arr1[matchBinPos(edges1, dat)] += 1
        for dat in data2:
            prob_arr2[matchBinPos(edges2, dat)] += 1
        if len(data1) > 0:
            prob_arr1 /= len(data1)
        if len(data2) > 0:
            prob_arr2 /= len(data2)
        sum1 = 0.0
        for x in prob_arr1:
            for y in prob_arr2:
                if x * y > 0:
                    sum1 += x * y * math.log(x * y, 2)
        sum1 = sum1 * -1
        return sum1

    dat_arr = list_empty([bins,2], [])
    binPos = createBinPos(Y, bins)

    for i, y in enumerate(Y):

        i_index = matchBinPos(binPos, y)
        dat_arr[i_index][0].append(X1[i])
        dat_arr[i_index][1].append(X2[i])

    X1Edges = createBinPos(X1, bins)
    X2Edges = createBinPos(X2, bins)
    HX = entropy2d(np.asarray(X1), np.asarray(X2), X1Edges, X2Edges)
    HXgivenY = 0.0
    Y_prob, edges = np.histogram(Y, bins, density = True)
    binwidth = edges[1] - edges[0]
    Y_prob = Y_prob*binwidth
    for i, y in enumerate(dat_arr):
        HXgivenY += Y_prob[i]*entropy2d(np.asarray(y[0]), np.asarray(y[1]), X1Edges, X2Edges)

    return HX-HXgivenY, HX, HXgivenY

def getTodaysDate():
    return date.today()

def getTime():
    from datetime import datetime
    now = datetime.now()
    return now.strftime("%H_%M_%S")

def createSingleDirectory(prepath, name):
    path = prepath +'/'+ name
    isPath = os.path.isdir(path)
    if not isPath:
        os.mkdir(path)

def createDirectory(run, date, config):
    path = config.runStats.saveDir + "/" + str(date) + "/" + run
    isPath = os.path.isdir(config.runStats.saveDir + "/" + str(date))
    if not isPath:
        os.mkdir(config.runStats.saveDir + "/" + str(date))
    isPath = os.path.isdir(path)
    if not isPath:
        os.mkdir(path)
    return str(date)

def saveData(dir, data, fileName):
    full_path = dir + "/" + fileName
    with open(full_path, 'wb') as f:
        pickle.dump(data, f)
    return full_path

def loadData(filename):
    with open(filename, 'rb') as f:
        new_data = pickle.load(f)
    return new_data

def saveDataDate(run, date, dataName, data_array, compress, prePath):
    path = prePath + "/" + date + "/" + run
    isPath = os.path.isdir(prePath + "/" + date)
    if not isPath:
        os.mkdir(prePath + "/" + date)
    isPath = os.path.isdir(path)
    if not isPath:
        os.mkdir(path)
    file_name =dataName + "_"+ str(datetime.datetime.now().time()).replace(":", "_")
    full_uncompressed_path = path + "/" +file_name + ".bin"
    full_compressed_path = path + "/" +file_name + ".gz"
    with open(full_uncompressed_path, 'wb') as f:
        pickle.dump(data_array, f)
    if compress:
        f_in = open(full_uncompressed_path, mode = 'rb')
        open(full_compressed_path, "w").close()
        f_out = gzip.open(full_compressed_path,  mode='wb', compresslevel=9, encoding=None, errors=None, newline=None)
        f_out.writelines(f_in)
        f_out.close()
        f_in.close()
        os.remove(full_uncompressed_path)
        return full_compressed_path
    else:
        return full_uncompressed_path

def loadDataDate(file_name, compressed, remove = False, removeCompress = False):
    path = ""
    if compressed:
        with gzip.open(file_name, 'rb') as f:
            file_content = f.read()
        tmpFileName = os.path.splitext(file_name)[0] + "tmp.bin"
        tmp = open(tmpFileName, "wb")
        tmp.write(file_content)
        tmp.close()
        with open(tmpFileName, 'rb') as f:
            new_data = pickle.load(f)
        path = tmpFileName
    else:
        with open(file_name, 'rb') as f:
            new_data = pickle.load(f)
        path = file_name
    if remove:
        os.remove(path)
    if removeCompress and compressed:
        os.remove(file_name)
    return new_data, path

def createGifBar(config, arr, celllocs, name):
    l.log("creating gif for: " + name)
    x = range(len(arr[0][0]))
    ims = []
    index = math.ceil(len(arr)/200)

    for i in range(0, len(arr), index):
        currTime = (i*config.simParams.simTimeStep)
        plt.clf()
        plt.ylim([0,100])
        plt.bar(x, arr[i][0], width=1, color = 'blue')
        plt.bar(x, arr[i][1], width=1, color = 'orange')
        plt.bar(x, celllocs[i], width = 1, color = 'red')
        plt.title("Time: " + str(currTime))
        plt.xlabel("Location")
        plt.ylabel("Concentration / Cell Count")
        buf = io.BytesIO()
        plt.savefig(buf, format='png')
        ims.append(buf)

    ims2 = []
    for buf in ims:
        buf.seek(0)
        ims2.append(imageio.imread(buf))
        buf.close()

    imageio.mimsave('sim_gifs/' + name + '.gif', ims2)

def loadMetaData(filename = None):
    metaDict = {}

    def Map(line):
        class metaFileException(Exception):
            pass
        arr = line.split('=')
        if len(arr) != 2:
            l.log("error: meta format is not correct on line: " + line)
            raise metaFileException()

        metaDict[arr[0]] = arr[1]


    if not filename:
        filename = "SimMeta.txt"
    with open(filename, "r") as file:
        count = 0
        for line in file:
            stripped_line = line.replace(" ", "").replace("\t", "").replace("\n", "")
            Map(stripped_line)

    return metaDict
