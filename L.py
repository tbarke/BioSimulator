from scipy import special
import cell
import enviornment
import utils
import simulation
import MICalc
import numpy
import random
import randomClass
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm
import math
import numpy as np
import gc
import run


SimParams = [100,-1]
# [diffCoeff, absorb/reflect/periodic, repeatFrequency, magnitude, locationA (negative if random), locationB]
ConcParams = [5, "periodic", 0, 100000.0, 35, 70, "static", 10 * ((2*math.pi)/SimParams[0]), 500, -25.5, 25.5, 0, 50]

#cell meta stats: [ Absorbtion rate, receptorConsumptionRate, survivalCost, velocities, mutate bool, decisiontype ("naive", "measured", "prediction", "non"-no change)]
    #Stats: [Arec, Brec, MaxRec, Amol, Bmol, ATP, biomass, generation, distTrav, StartID]
    #cell location Stats: [array of Locations (empty if everywhere)(zero indexed), #TODO add later: cell counts at the locations]

CellMetaStats = [5, 0.0, 5, 1.0, "non", "non"]
CellStats = [200,200,400,0,0,0,5,1,0,0]
Celllocations = []
CellLocationStats = [Celllocations]

# Enviornemnt paramter array: [simulationLengthTime, full divide, full Die]
EnviornmentParams = [10, False, True]

CellMetaStats[5] = "measured"
simLength = EnviornmentParams[0]
simRepeat = 1
fianlTimeDivx = range(simLength)
finalTimeDiv = [0]*simLength
survivaly = []

Amol = []
Bmol = []

def findCoord(Amol, Bmol):
    if(Amol >= max):
        return -1, -1
    if(Bmol >= max):
        return -1, -1
    Aeffective = float(Amol)/float(max/divisions)
    Beffective = float(Bmol)/float(max/divisions)
    return math.trunc(Aeffective), math.trunc(Beffective)

locs = []
divides = []
Amol = []
splitLoc = []
Bmol = []
concs = []
for i in range(simRepeat):
    print(str(i) + ": Repeat Number")
    init = run.moveTrackRun(SimParams, ConcParams, CellMetaStats, CellStats, Celllocations, EnviornmentParams, True, False)
    locs = init[0]
    divides = init[1]
    splitLoc = init[2]
    Amol = init[3]
    Bmol = init[4]
    Acons = init[5]
    Bcons = init[6]

print(Acons)
print(Bcons)

print(Amol)
print(Bmol)
print(divides)
count = 0;
new_Amol = [];
new_Bmol = [];
for i in range(len(Amol)):
    if Amol[i] < 100:
        if Bmol[i] < 100:
            new_Amol.append(Amol[i])
            new_Bmol.append(Bmol[i])
    if Amol[i] < 10:
        count = count + 1
    if Bmol[i] <10:
        count = count + 1
print(count)

#for i in range(len(new_Amol)):
#    new_Amol[i] = new_Amol[i]/len(new_Amol)
#    new_Bmol[i] = new_Bmol[i]/len(new_Bmol)


#fig, ax = plt.subplots()
#h = plt.hist2d(new_Amol, new_Bmol, bins = 25)
#cbar = fig.colorbar(h[3], ax=ax)
#lab = "Cumulative cell density\n(time-step $200$)"
#cbar.set_label(lab, rotation = 270,labelpad = 25)
#plt.axis('square')
#plt.xlim(0, 100)
#plt.ylim(0, 100)
#x1, y1 = [25, 100], [25, 25]
#x2, y2 = [25, 25], [25, 100]
#plt.plot(x1, y1, x2, y2, color="black", lw = 2, alpha = 0.8)
##plt.plot(Amol, Bmol, marker = 'o')
#plt.xlabel("Internal A Molecule Count")
#plt.ylabel("Internal B Molecule Count")
#plt.show()


#Amol = [12, 15, 12, 9, 7, 6, 9, 12, 31, 13, 16, 16, 14, 12, 10, 8, 5, 3]
#Bmol = [15, 32, 97, 161, 224, 266, 283, 300, 302, 148, 198, 225, 276, 335, 394, 453, 518, 580]
#divides = [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]

print("divides len " + str(len(divides)))
print("Amol len " + str(len(Amol)))

print(Amol)
print(Bmol)
print(divides)
plt.rc('axes', axisbelow=True)
currArrA = [Amol[0], Amol[1]]
currArrB = [Bmol[0], Bmol[1]]
#plt.grid()
#plt.plot(currArrA, currArrB, color='orange')
plt.axis('square')
plt.arrow(Amol[0], Bmol[0], Amol[1]-Amol[0], Bmol[1] - Bmol[0], color = 'orange', head_width = 2, head_length = 2, length_includes_head = True)
currArrA = []
currArrB = []
#index = 1
#plt.show()
#exit(-1)
for i in range(1,len(divides)-1,1):

    if divides[i] == 1:
        #print(index)
        #currArrA.append(Amol[index])
        #currArrB.append(Bmol[index])
        #plt.plot(currArrA, currArrB, color = 'orange')
        plt.arrow(Amol[i], Bmol[i], Amol[i+1] - Amol[i], Bmol[i+1] - Bmol[i], color='black', head_width=2, head_length=2,length_includes_head=True)
        #currArrA = [Amol[index], Amol[index+1]]
        #currArrB = [Bmol[index], Bmol[index+1]]
        #plt.plot(currArrA, currArrB, color='black', marker='o')
        #print(currArrA)
        #print(currArrB)
        #currArrA = [Amol[index+1]]
        #currArrB = [Bmol[index+1]]
    else:
        plt.arrow(Amol[i], Bmol[i], Amol[i+1] - Amol[i], Bmol[i+1] - Bmol[i], color='orange', head_width=2, head_length=2,length_includes_head=True)
    #index = index + 1
#plt.plot(currArrA, currArrB, color='orange')

plt.xlim(0, 100)
plt.ylim(0, 100)
x1, y1 = [25, 100], [25, 25]
x2, y2 = [25, 25], [25, 100]
plt.plot(x1, y1, x2, y2, color="black", lw = 2)
#ssplt.plot(Amol, Bmol, marker = 'o')
plt.xlabel("Internal A Molecule Count")
plt.ylabel("Internal B Molecule Count")
plt.grid()


plt.show()
exit(-1)

colors = ['#4a4bbf','#6f51c2','#7c47be','#a24fc1','#af4abf','#bc48bf','#bd43a2','#be45ad','#ba418e','#be477f','#c14f80','#bd435e', '#bf4a4e']
Ac_high = utils.highest(Acons)
Bc_high = utils.highest(Bcons)

#print(Amol)
#print(Bmol)

#plt.plot(Amol, Bmol)


#for i in range(len(Amol)):
    #plt.plot(Amol[i], Bmol[i],marker='o', color = (Acons[i] / Ac_high,0.0,Bcons[i] / Bc_high))
    #if i != len(Amol) -1:
        #plt.arrow(Amol[i], Bmol[i], Amol[i+1]-Amol[i], Bmol[i+1]-Bmol[i], head_width = 4 , ec ='black')

max = 100
divisions = 100

ext_A = [] * divisions
ext_A_count = [] * divisions
ext_B = [] * divisions
ext_B_count = [] * divisions
for i in range(divisions):
    ext_A.append([0.0]* divisions)
    ext_A_count.append([0] *divisions)
    ext_B.append([0.0]* divisions)
    ext_B_count.append([0] *divisions)

for i in range(len(Amol)):
    A_in, B_in = findCoord(Amol[i], Bmol[i])
    if A_in >= 0:
        ext_A[A_in][B_in] = ext_A[A_in][B_in] + Acons[i]
        ext_B[A_in][B_in] = ext_B[A_in][B_in] + Bcons[i]
        ext_B_count[A_in][B_in] = ext_B_count[A_in][B_in] + 1
        ext_A_count[A_in][B_in] = ext_A_count[A_in][B_in] + 1

#for i in range(len(ext_A)):
#    for j in range(len(ext_A[i])):
#        if ext_A_count[i][j] != 0:
#            ext_A[i][j] = ext_A[i][j]/ext_A_count[i][j]
#        if ext_B_count[i][j] != 0:
#            ext_B[i][j] = ext_B[i][j]/ext_B_count[i][j]

Aall = []
Ball = []
high = Ac_high
#a = [high,0]
#b = [0,high]
#h = plt.hist2d(a,b,bins = 50)
#fig, ax = plt.subplots()
#fig.colorbar(h[3], ax=ax)
for i in range(len(ext_A)):
    for j in range(len(ext_B[i])):
        #print(ext_A[i][j])
        for k in range(int(ext_A[i][j])):
            Aall.append(i)
            Ball.append(j)

        #shade = (ext_A[i][j]/high)
        #index_i = i*max/divisions
        #index_j = j * max / divisions
        #x = [index_i, index_i + (max/divisions)]
        #x = [index_i, index_i + (max/divisions)]
        #y = [index_j + (max/divisions), index_j + (max/divisions)]
        #y2 = [index_j, index_j]
        #y = [index_j, index_j + (max/divisions)]
        #plt.fill(x, y, color = (0.5*shade,0,0.5*(1-shade)), alpha = shade)
        #ax.fill_between(x, y, y2, facecolor = (0.188+((1-0.188)*shade),0.098+((1-0.098)*shade),0.2039-(0.2039*shade)))
        #plt.fill_between(x, y, y2, facecolor = (0,0,shade))

fig, ax = plt.subplots()
h = plt.hist2d(Aall,Ball,bins = 25)
cbar = fig.colorbar(h[3], ax=ax)
lab = "Cumulative external A concentration\n(time-step $200$)"
cbar.set_label(lab, rotation = 270,labelpad = 25)

#fig, ax = plt.subplots()
#h = plt.hist2d(new_Amol, new_Bmol, bins = 50)
#fig.colorbar(h[3], ax=ax)
plt.axis('square')
plt.xlim(0, max)
plt.ylim(0, max)
x1, y1 = [25, max], [25, 25]
x2, y2 = [25, 25], [25, max]
plt.plot(x1, y1, x2, y2, color="black", lw = 2, alpha = 0.8)
plt.xlabel("Internal A molecule count")
plt.ylabel("Internal B molecule count")

#legend_elements = [Line2D([0], [0], marker='o', color='w', label='High External A Concentration', markerfacecolor='red', markersize=15), Line2D([0], [0], marker='o', color='w', label='High External B Concentration', markerfacecolor='blue', markersize=15)]
#plt.legend(handles=legend_elements)
#plt.grid()
plt.show()

Aall = []
Ball = []
high = Ac_high
#a = [high,0]
#b = [0,high]
#h = plt.hist2d(a,b,bins = 50)
#fig, ax = plt.subplots()
#fig.colorbar(h[3], ax=ax)
for i in range(len(ext_A)):
    for j in range(len(ext_B[i])):
        #print(ext_A[i][j])
        for k in range(int(ext_B[i][j])):
            Aall.append(i)
            Ball.append(j)

fig, ax = plt.subplots()
h = plt.hist2d(Aall,Ball,bins = 25)
cbar = fig.colorbar(h[3], ax=ax)
lab = "Cumulative external B concentration\n(time-step $200$)"
cbar.set_label(lab, rotation = 270,labelpad = 25)

#fig, ax = plt.subplots()
#h = plt.hist2d(new_Amol, new_Bmol, bins = 50)
#fig.colorbar(h[3], ax=ax)
plt.axis('square')
plt.xlim(0, max)
plt.ylim(0, max)
x1, y1 = [25, max], [25, 25]
x2, y2 = [25, 25], [25, max]
plt.plot(x1, y1, x2, y2, color="black", lw = 2, alpha = 0.8)
plt.xlabel("Internal A molecule count")
plt.ylabel("Internal B molecule count")
plt.show()

exit(-1)

max = 200
divisions = 50
print(len(Amol))

gradients = [] * divisions
gradients_count = [] * divisions
for i in range(divisions):
    gradients.append([0.0]* divisions)
    gradients_count.append([0] *divisions)

for i in range(len(Amol)-1):
    A_in, B_in = findCoord(Amol[i], Bmol[i])
    if A_in >= 0:
        if float(Amol[i+1] - Amol[i]) == 0:
            slope = 1
            print("this happend")
        else:
            #slope =  float(Bmol[i + 1] - Bmol[i])/float(Amol[i + 1] - Amol[i])
            slope = math.atan(math.fabs(float(Bmol[i + 1] - Bmol[i]))/math.fabs(float(Amol[i + 1] - Amol[i])))
            if Bmol[i+1] - Bmol[i] < 0:
                if Amol[i+1] - Amol[i] < 0:
                    slope = slope + math.pi
                else:
                    slope = 2*math.pi - slope
            else:
                if Amol[i+1] - Amol[i] < 0:
                    slope = math.pi - slope
            if A_in > 1:
                if B_in > 1:
                    print(slope)
        gradients[A_in][B_in] = gradients[A_in][B_in] + slope
        gradients_count[A_in][B_in] = gradients_count[A_in][B_in] + 1

for i in range(len(gradients)):
    curr_gradient = gradients[i]
    curr_gradient_count = gradients_count[i]
    for j in range(len(curr_gradient)):
        if curr_gradient_count[j] > 0:
            curr_gradient[j] = float(curr_gradient[j])/float(curr_gradient_count[j])

largest_gradient = 0.0
for i in range(len(gradients)):
    for j in range(len(gradients[i])):
        if gradients[i][j] > largest_gradient:
            largest_gradient = gradients[i][j]


for i in range(divisions):
    for j in range(divisions):
        if gradients[i][j] != 0.0:
            #print(gradients[i][j])
            magnitude = (max / (divisions*2))+(max/(divisions*2))
            Arrow_magnitude = (max/(divisions*3.5*2))+(max/(divisions*3.5*2))

            A_new = (magnitude)*math.cos(gradients[i][j])
            B_new = (magnitude)*math.sin(gradients[i][j])
            #plt.arrow(i*max/divisions, j*max/divisions, A_new, B_new, head_width=Arrow_magnitude, color = 'black', length_includes_head = True, ec='black')
            index_i = i*(max/divisions) + max/(divisions*2) - A_new/2
            index_j = j*(max/divisions) + max/(divisions*2) - B_new/2
            plt.arrow(index_i, index_j, A_new, B_new, head_width=Arrow_magnitude, color = 'black', length_includes_head = True, ec='black')
        if gradients[i][j] == 0.0:
            #print(gradients[i][j])
            plt.plot(i*max/divisions,j*max/divisions, marker ='o', color ='black', markersize = 1)

plt.xlabel("Internal A Molecule Count")
plt.ylabel("Internal B Molecule Count")
plt.xlim(-1, max)
plt.ylim(-1, max)
x1, y1 = [25, max], [25, 25]
x2, y2 = [25, 25], [25, max]
plt.plot(x1, y1, x2, y2, color="red", lw = 2)
plt.grid()
plt.show()

#plt.xlim(0, 200)
#plt.ylim(0, 200)
#plt.xlabel("Internal A Molecule Count")
#plt.ylabel("Internal B Molecule Count")

#plt.grid()

#plt.show()