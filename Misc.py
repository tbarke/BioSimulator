# single array variance code-------------------------------------------
def get_growth_entropy(totalCells, entropy ):
    growth_arr = []
    all_entropy = []
    for i in range(math.ceil(1 / EnviornmentParams[3]), len(totalCells), math.ceil(1 / EnviornmentParams[3])):
        growth_arr.append(math.log2(totalCells[i] / totalCells[i - math.ceil(1 / EnviornmentParams[3])]))
        all_entropy.append(entropy[i])
    ave_growth = 0.0
    ave_entropy = 0.0
    for i in range(len(growth_arr)):
        ave_growth = ave_growth + growth_arr[i]
        ave_entropy = ave_entropy + all_entropy[i]

    return ave_growth, ave_entropy

def parse_arr(string_g):
    string_g = string_g.replace(" ", "")
    string_g = string_g.replace("[", "")
    string_g = string_g.replace("]", "")
    arr_g = string_g.split(",")
    num_arr_g = []
    for i in range(len(arr_g)):
        num_arr_g.append(float(arr_g[i]))
    return num_arr_g

f = open("Data_8_1_2022.txt", "r")
i = 0
growth_arr_num = []
ave_growth_arr = []
ave_entropy_arr = []
num_MI_arr = []
ratio_arr = []

globali = -1
for x in f:
    globali = globali + 1
    if i == 0:
        str_ratio = x.replace(" ", "")
        ratio_arr.append(float(str_ratio))
    if i == 1:
        growth_arr_num = parse_arr(x)
    if i == 2:
        ave_growth, ave_entropy = get_growth_entropy(growth_arr_num, parse_arr(x))
        ave_growth_arr.append(ave_growth)
        ave_entropy_arr.append(ave_entropy)
    if i == 3:
        str_MI = x.replace(" ", "")
        num_MI = float(str_MI)
        num_MI_arr.append(num_MI)
    i = (i + 1) % 9
# single array variance code ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#TestSuite stuff -------------------------------------------------------
presets =[]
presets.append([])
presets.append([])
presets.append([])
presets.append([])

import sys
for i in range(10):
    for j in range(5):
        for k in range(5):
            adaptive_ratio = j*0.2
            CellMetaStats[11] = adaptive_ratio
            if k == 0 and j == 0:
                preset_bool = False
            else:
                preset_bool = True
            totalCells, entropy, MI, presets = run.testRun(SimParams, ConcParams, CellMetaStats, CellStats,
                                                           Celllocations, EnviornmentParams, True, preset_bool, presets)
            original_stdout = sys.stdout
            with open('Data_8_1_2022.txt', 'a') as f:
                sys.stdout = f  # Change the standard output to the file we created.
                print(adaptive_ratio)
                print(totalCells)
                print(entropy)
                print(MI[0])
                print(presets[0])
                print(presets[1])
                print(presets[2])
                print(presets[3])
                print()
                sys.stdout = original_stdout

#change cell meta stats[11]
#totalCells, entropy, MI, presets = run.testRun(SimParams, ConcParams, CellMetaStats, CellStats, Celllocations, EnviornmentParams, True)

all_growth = []
all_entropy = []
time = []
count = -1
for i in range(math.ceil(1/EnviornmentParams[3]), len(totalCells), math.ceil(1/EnviornmentParams[3])):
    count = count+1
    time.append(count)
    all_growth.append(math.log2(totalCells[i]/totalCells[i-math.ceil(1/EnviornmentParams[3])]))
    all_entropy.append(entropy[i])

print(time)
print(all_growth)
print(MI[0])
ave_growth = 0.0
ave_entropy = 0.0
for i in range(len(all_growth)):
    ave_growth = ave_growth + all_growth[i]
    ave_entropy = ave_entropy + all_entropy[i]
ave_growth = ave_growth/len(all_growth)
ave_entropy = ave_entropy/len(all_entropy)

print(ave_growth)
print(ave_entropy)
exit(-1)
print(totalCells)
sum = 0
growth = []
gprime = []

for i in range(len(totalCells)-1):
    if totalCells[i] == 0:
        growth.append(0)
        sum += 0
    else:
        growth.append(math.log2(totalCells[i+1]/totalCells[i]))
        sum += math.log2(totalCells[i+1]/totalCells[i])

for i in range(len(growth)-1):
    if growth[i] == 0:
        gprime.append(0)
    else:
        gprime.append(growth[i+1]/growth[i])

final_growth = []
final_entropy = []
for i in range(len(gprime)):
    if gprime[i] > 0:
        final_growth.append(gprime[i])
        final_entropy.append(entropy[i])

print(final_entropy)
print(final_growth)

plt.grid()
plt.scatter(final_growth, final_entropy)
plt.xlabel("G'")
plt.ylabel("$Entropy(MI_{state})$")
plt.title("Stress: 2.5")
plt.show()
plt.clf()
exit(-1)

print(sum/len(totalCells))
#non  = [11.26043464362008, 12.120183411334366, 34.710449085353424, 25.08943314051025, 20.67708409040647, 23.27401659725306, 26.616194229917873, 44.2524317327881, 21.815797149698632]
#prediciton  = [20.90962342825688, 38.773785241294235, 24.359919334650012, 34.11222834257083, 27.234352259209345, 12.470726076180798, 23.94530478149776, 30.56467034660499, 10.794854901274366]

#non_sum = 0.0
#prediciton_sum = 0.0
#print(len(non))
#print(len(prediciton))
#for i in range(len(non)):
#    non_sum += non[i]/len(non)
#    prediciton_sum += prediciton[i]/len(prediciton)

#print(non_sum)
#print(prediciton_sum)

#TestSuite stuff ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^