import math
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import MICalc



A = [ 2.87215757, 3.17401621, 3.50621655, 3.8701379,  4.26680696, 4.69679018
, 5.16007924, 5.65597384, 6.18296729, 6.73864149, 7.31957889, 7.92129975
, 8.53823294, 9.16372833, 9.79011719 ,10.40882549, 11.01054163, 11.58543739,
 12.12343668, 12.61452277, 13.04907136, 13.41819327, 13.71406885, 13.93025513,
 14.06194744, 14.10617924, 14.06194744, 13.93025513, 13.71406885, 13.41819327,
 13.04907136, 12.61452277, 12.12343668, 11.58543739, 11.01054163, 10.40882549
, 9.79011719, 9.16372833, 8.53823294, 7.92129975, 7.31957889, 6.73864149
, 6.18296729, 5.65597384, 5.16007924, 4.69679018, 4.26680696, 3.8701379
, 3.50621655, 3.17401621, 2.87215757, 2.59900661, 2.35276087, 2.1315233
, 1.93336356, 1.75636739, 1.59867489, 1.45850906, 1.33419582, 1.22417688
, 1.1270169,  1.041406, 0.96615883, 0.90021101, 0.84261393, 0.79252833
, 0.74921738, 0.71203951, 0.68044147, 0.65395174, 0.63217442, 0.61478389
, 0.60152018, 0.59218507, 0.58663917, 0.58479968, 0.58663917, 0.59218507
, 0.60152018, 0.61478389, 0.63217442, 0.65395174, 0.68044147, 0.71203951
, 0.74921738, 0.79252833, 0.84261393, 0.90021101, 0.96615883, 1.041406
, 1.1270169,  1.22417688, 1.33419582, 1.45850906, 1.59867489, 1.75636739
, 1.93336356, 2.1315233,  2.35276087, 2.59900661]

B = [ 2.87215757, 2.59900661, 2.35276087, 2.1315233,  1.93336356, 1.75636739
, 1.59867489, 1.45850906, 1.33419582, 1.22417688, 1.1270169,  1.041406
, 0.96615883, 0.90021101, 0.84261393, 0.79252833, 0.74921738, 0.71203951
, 0.68044147, 0.65395174, 0.63217442, 0.61478389, 0.60152018, 0.59218507
, 0.58663917, 0.58479968, 0.58663917, 0.59218507, 0.60152018, 0.61478389
, 0.63217442, 0.65395174, 0.68044147, 0.71203951, 0.74921738, 0.79252833
, 0.84261393, 0.90021101, 0.96615883, 1.041406, 1.1270169,  1.22417688
, 1.33419582, 1.45850906, 1.59867489, 1.75636739, 1.93336356, 2.1315233
, 2.35276087, 2.59900661, 2.87215757, 3.17401621, 3.50621655, 3.8701379
, 4.26680696, 4.69679018, 5.16007924, 5.65597384, 6.18296729, 6.73864149
, 7.31957889, 7.92129975, 8.53823294, 9.16372833, 9.79011719, 10.40882549,
 11.01054163 ,11.58543739, 12.12343668, 12.61452277 ,13.04907136 ,13.41819327,
 13.71406885, 13.93025513 ,14.06194744 ,14.10617924, 14.06194744, 13.93025513,
 13.71406885 ,13.41819327 ,13.04907136, 12.61452277, 12.12343668, 11.58543739,
 11.01054163 ,10.40882549, 9.79011719, 9.16372833, 8.53823294, 7.92129975
, 7.31957889, 6.73864149, 6.18296729, 5.65597384, 5.16007924, 4.69679018
, 4.26680696, 3.8701379,  3.50621655, 3.17401621]

stress = [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.8, 0.9, 0.95, 1]
stress = []
tot = 1000
max = 1
for i in range(tot):
 print(0.2 + i * ((max - 0.2) / (tot-1)))
 stress.append(0.2+i*((max-0.2)/(tot-1)))

#stress = [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.8, 0.9, 0.95, 1]
#stress = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
growth = []
for i in range(len(stress)):
 total_growth = 0.0
 absorb = 1/stress[i]
 sur = 1.0
 curr_growth = 0.0
 n = len(A)
 for j in range(len(A)):
  x1 = A[j]
  x2 = B[j]
  #print()
  #print(absorb)
  #print(sur)
  #print(x1)
  #print(x2)
  if x1 < x2:
   rate = (absorb * x1 - sur) / (2.5 * sur)
   #print("rate: " + str(rate))
   if(rate <= 0.0):
    curr_growth = 0
    n = n -1
   else:
    curr_growth = math.log(rate,2)
  else:
   rate = (absorb * x1 - sur) / (2.5 * sur)
   #print("rate: " + str(rate))
   if rate <= 0.0:
    curr_growth = 0.0
    n = n - 1
   else:
    curr_growth = math.log(rate, 2)
  total_growth = total_growth + curr_growth
 print(n)
 total_growth = total_growth/n
 #print(total_growth)
 growth.append(total_growth)
plt.plot(stress,growth)
plt.xlabel("cell Stress: (1/absorption rate)")
plt.ylabel("Cell Growth rate")
plt.grid()
plt.show()

exit(-1)

receptors = 100
prob = 0.0001
std = math.sqrt(receptors * prob * (1 - prob))

probs = []
for i in range(len(concA)):
 prob = concA[i] / (2.0 + concA[i])
 probs.append(prob)
plt.hist(probs)
plt.show()
plt.clf()

stds = []
ran = 100
final_std = 0.0
for i in range(len(probs)):
 prob = probs[i]
 stds.append(math.sqrt(receptors * prob * (1 - prob)))
 final_std = final_std + math.sqrt(receptors * prob * (1 - prob))
print(final_std/len(probs))
#plt.hist(stds)
#plt.show()
exit(-1)

print(std)
std = std * 4
mean = receptors * prob
flag = True
samples = []
for i in range(10000):
 sample = np.random.normal(mean, std, 1)[0]
 if sample >= 0 and sample < receptors:
  samples.append(sample)

plt.hist(samples)
plt.show()
exit(-1)


testMI = MICalc.MICalc()
ret = testMI.AltMI(dataX,dataY)



n = 100
p = 0.2

std = math.sqrt(n*p*(1-p))
noise_array = [0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
mean = n*p
ret = []
for i in range(len(noise_array)):
 ret = []
 std = math.sqrt(n * p * (1 - p))
 std = std * noise_array[i]
 for j in range(1000):
  curr = np.random.normal(mean, std, 1)
  ret.append(curr[0])
#ret = np.random.normal(mean, std, 1000)
 print(ret)

 plt.hist(ret)
 plt.show()
exit(-1)

# Create Map
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# Make data.
noise_array = [0.0, 0.3, 0.7, 1]
cell_stress_array = [0.2, 0.4, 0.6, 0.8]
full_measured_growth = [[6.932529919349607, 5.844643785650846, 4.555216445560742, 3.5004600041087497], [6.995793093552119, 5.811183859424624, 4.633669365812751, 3.5325763866386106], [6.9409292456610725, 5.915455510002327, 4.479599809860369, 3.447238940167634], [6.908560964361998, 5.878963678823775, 4.588215586444321, 3.532219778053206]]

Xs = np.array(noise_array)
Ys = np.array(cell_stress_array)
Xs, Ys = np.meshgrid(Xs,Ys)
Z = []
for i in range(len(Xs)):
 newZ = []
 for j in range(len(Ys)):
  newZ.append(full_measured_growth[i][j])
 Z.append(newZ)

print(Xs)
print(Z)

Zs = np.array(Z)
#print(Zs)
#X = np.array(Xs)
#Y = np.array(Ys)
#Z = np.array(Zs)
#print(Z)
#X = np.arange(-5, 5, 0.25)
#Y = np.arange(-5, 5, 0.25)
#X, Y = np.meshgrid(X, Y)
#R = np.sqrt(X**2 + Y**2)
#Z = np.sin(R)
#print(Z)

# Plot the surface.
surf = ax.plot_surface(Xs, Ys, Zs, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
#ax.set_zlim(-1.01, 1.01)
#ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
#ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

exit(-1)

n = 100
p = 0.2

std = math.sqrt(n*p*(1-p))
std = 0
mean = n*p
ret = []
for i in range(1):
 curr = np.random.normal(mean, std, 1)
 print(curr[0])
 ret.append(curr[0])
ret = np.random.normal(mean, std, 1000)
print(ret)

plt.hist(ret)
plt.show()
exit(-1)






