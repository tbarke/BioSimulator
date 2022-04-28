import cell
import utils
import MICalc
import math
from scipy import special
import matplotlib.pyplot as plt
import numpy as np

#number of receptors either A or B on both sides of the cel


x = [1,2]
y = [2,2]
y2 = [1,1]

a = [1,0]
b = [0,1]
h = plt.hist2d(a,b,bins = 50)
#plt.clf()
#y = [index_j, index_j + (max/divisions)]
#plt.fill(x, y, color = (0.5*shade,0,0.5*(1-shade)), alpha = shade)
fig, ax = plt.subplots()
ax.fill_between(x, y, y2, facecolor = (0,0,1,0.5), edgecolor = (0,0,0,0))

fig.colorbar(h[3], ax=ax)
plt.show()
exit(-1)


def VonMises( std, maximum, offset, length):
    length = length + 1
    arr = []
    kappa = 1 / std
    mu = 2 * math.pi / length * offset
    bessel = (special.iv(0, kappa))
    sum = 0.0
    increment = 2 * math.pi / length
    for i in range(math.floor(-1 * (length - (length / 2))), math.floor(length - (length / 2))):
        x = i * increment
        numerator = math.exp(kappa * math.cos(x - mu))
        arr.append(numerator / (2 * math.pi * bessel))
        sum += numerator / (2 * math.pi * bessel)
    arr = (arr / sum) * maximum
    return arr

def output(array, receptors, file):
    f = open(file, "a")
    f.write(str(receptors))
    f.write(",")
    for i in array:
        f.write(str(i))
        f.write(",")
    f.write("\n")
    f.close()

def runSNR(conc):
    recNum = 100
    dissocociationConstant = 2
    expectedBound = []
    rep = 1000
    for k in range(200):
        SNR = []
        recNum = k
        print(str(recNum)+"-----------------------------------------------------------------------------------")
        for j in range(len(conc)):
            expectedBound = recNum * (conc[j] / (dissocociationConstant + conc[j]))
            expected = 0
            error = 0
            print(j)
            for i in range(rep):
                numBound = utils.calculateBoundKinetic(recNum, conc[j], dissocociationConstant, 0)
                error += math.pow((expectedBound - numBound), 2)
                expected += math.pow(expectedBound, 2)
            error = math.sqrt(error / rep)
            expected = math.sqrt(expected / rep)
            if error == 0.0:
                SNR.append(0.0)
            else:
                SNR.append(math.pow(expected / error, 2))

        for i in range(len(SNR)):
            if SNR[i] > 0:
                SNR[i] = 10*math.log(SNR[i],10)

        output(SNR, recNum, "allSNRtable.txt")

def plotSNR(conc):
    f = open("Output/SNR/SNRperReceptor.txt", "r")
    lines = f.readlines()
    x = []
    rec = []
    for i in range(99):
        x.append(i)
    arrays = []
    print(len(lines))
    for i in range(0,len(lines), 30):
        arrline = lines[i].split(",")
        arrline = arrline[0:len(arrline)-1]
        for j in range(len(arrline)):
            arrline[j] = float(arrline[j])

            array = arrline[1:len(arrline)]
        rec.append(arrline[0])
        arrays.append(array)
    f.close()
    handles1 = []
    for i in range(1,len(arrays)):
        #p, = plt.scatter(conc, arrays[i], label = rec[i])
        plt.scatter(conc, arrays[i], label = rec[i], marker = '.')
        curr_coefficients = np.polyfit(conc, arrays[i], 3)
        curr_poly = np.poly1d(curr_coefficients)
        new_x = np.linspace(conc[0], conc[-1])
        new_y = curr_poly(new_x)

        #new_new_x = new_x.tolist()
        #new_new_y = new_y.tolist()
        #print(type(new_new_x))

        p, = plt.plot(new_x, new_y, label = rec[i])
        handles1.append(p)
    plt.xlabel("concentration")
    plt.ylabel("dB")


    #new_x = np.linspace(x[0], x[-1])
    #Calculate
    #new
    #x and y
    #values


    #

    #plt.plot(x, y, "o", new_x, new_y)

    #rec[0] = str(rec[0]) + " receptors"
    print(rec)
    handles1 = handles1[::-1]
    #for i in handles1:
    #    handles2.append(handles1[len(handles1)-1-i])
    plt.legend(handles=handles1, title='receptors', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.show()

def signalDerivative(conc):
    total_points = 1000
    interpolated = total_points/len(conc);
    interpolated_int = math.floor(total_points/len(conc));
    new_conc = []
    for i in range(1, len(conc)-1, 2):
        delta = conc[i+1] - conc[i]
        for j in range(interpolated_int):
            #print(delta * 1/interpolated*j)
            new_conc.append(conc[i] + (delta * 1/interpolated*j))

    max = conc[0]
    min = conc[0]
    for i in conc:
        if i < min:
            min = i
        if i > max:
            max = i
    delta = max - min
    print(min)
    linear_conc = []
    for i in range(total_points):
        linear_conc.append(min + (delta * 1/total_points * i))


    for i in range(len(linear_conc)):
        if linear_conc[i] != 0:
            linear_conc[i] = math.log(linear_conc[i], 10)
    rec_num = 50
    dissocociationConstant = 2
    expected = []
    for i in linear_conc:
        expected.append(rec_num * (i / (dissocociationConstant + i)))

    derivative = []
    for i in range(0,len(expected)-1,2):
        derivative.append((expected[i+1] - expected[i])/(expected[i+1] - expected[i]))

    #print(derivative)
    plt.plot(linear_conc)
    plt.show()


conc = VonMises( 10 * ((2*math.pi)/100),  500, 0, 100)
conc.sort()

testMI = MICalc.MICalc()
print(testMI.returnSingleEntropy(conc,conc)*4)

#plotSNR(conc)
#signalDerivative(conc)
#runSNR(conc)
#plotSNR(conc)
#RMS = @(x) sqrt(mean(x.^2));
#x_snr = (RMS(x) / RMS(n)) ^ 2;

#poisson SNR = N/sqrt(N)