import matplotlib.pyplot as plt
import utils

#only combined realocation
path = "data/2022-04-25/test_run_03_21_11.841302tmp.bin"
path = "Data/2022-04-25/5_5_run_0.66Combined_03_41_58.055640tmp.bin"
path = "Data/2022-04-25/5_5_run_0.66Combined_03_52_03.555464tmp.bin"
path = "Data/2022-04-25/5_5_run_0.0Combined_03_55_42.571186tmp.bin"

#divide_combine also
path = "Data/2022-04-25/5_5_run_0.66Combined_both_04_00_53.141518tmp.bin"
path = "Data/2022-04-25/5_5_run_0.33Combined_both_04_04_15.704694tmp.bin"
path = "Data/2022-04-25/5_5_run_0.0Combined_both_04_25_41.387294tmp.bin"

path = "Data/2022-04-25/Time_Step_Test_06_23_46.985354tmp.bin"
path = "Data/2022-04-25/Time_Step_Test_06_39_27.521796tmp.bin"
path = "Data/2022-04-25/Time_Step_Test_06_44_04.716131tmp.bin"
path = "Data/2022-04-25/Time_Step_Test_final_06_49_31.740073tmp.bin"
path = "Data/2022-04-25/Time_Step_Test_final_06_55_07.230066tmp.bin"
path = "Data/2022-04-25/Time_Step_Test_final_07_27_56.616540tmp.bin"
path = "Data/2022-04-25/Time_Step_Test_final_07_38_36.203365tmp.bin"
path = "Data/2022-04-25/Time_Step_Test_final_07_47_15.740381tmp.bin"

#25.5
path = "Data/2022-04-25/Time_Step_Test_final_07_59_28.289724tmp.bin"
#10
#path = "Data/2022-04-25/Time_Step_Test_final_08_05_10.389015tmp.bin"
#5
#path = "Data/2022-04-25/Time_Step_Test_final_08_12_09.786105tmp.bin"
#0
#path = "Data/2022-04-25/Time_Step_Test_final_07_54_30.356977tmp.bin"


#10
#path = "Data/2022-04-25/Time_Step_Test_final_08_37_01.713490tmp.bin"

#path = "Data/2022-04-25/Time_Step_Test_final_08_44_44.283335tmp.bin"

#----------------------------------------------------
#25.5
#path = "Data/2022-04-25/Time_Step_Test_final_08_49_05.422902tmp.bin"
#10
#path = "Data/2022-04-25/Time_Step_Test_final_08_54_30.536347tmp.bin"
#5
#path = "Data/2022-04-25/Time_Step_Test_final_08_31_36.061870tmp.bin"
#0
#path = "Data/2022-04-25/Time_Step_Test_final_08_24_37.989781tmp.bin"

#combined 0.66
#path = "Data/2022-04-25/Time_Step_Test_final_09_21_30.079296tmp.bin"

#combined 0.0
#path = "Data/2022-04-25/Time_Step_Test_final_09_31_29.660279tmp.bin"

#noisless channel
path = "Data/2022-04-25/Time_Step_Test_final_09_37_35.166590tmp.bin"


path = "Data/2022-05-16/Time_Step_Test_final_06_17_02.418992tmp.bin"

path = "Data/2022-05-16/Time_Step_Test_final_06_20_18.103886tmp.bin"
data, path = utils.loadDataDate(path, False)
noise = data[2]
cell_stress_array = data[3]

growth_measured = data[4]
growth_max = data[5]

MI_measured = data[6]
MI_max = data[7]

growth_measured_full = []
growth_max_full = []
MI_measured_full = []
MI_max_full = []
for i in range(len(growth_measured[0])):
    growth_measured_full.append(growth_measured[0][i])
    growth_max_full.append(growth_max[0][i])

    MI_measured_full.append(MI_measured[0][i])
    MI_max_full.append(MI_max[0][i])

plt.plot(cell_stress_array, growth_measured_full,marker = 'o', label = "Adaptive Receptor Growth")
plt.plot(cell_stress_array, growth_max_full,marker = 'o', label = "Equal Receptor Growth")
plt.xlabel("Cell Stress")
plt.ylabel("Cell Growth")
plt.title("Growth")
plt.grid()
plt.legend()
plt.show()
plt.clf()

plt.plot(cell_stress_array, MI_measured_full,marker = 'o', label = "Adaptive Receptor MI")
plt.plot(cell_stress_array, MI_max_full,marker = 'o', label = "Equal Receptor MI")
plt.xlabel("Cell Stress")
plt.ylabel("Cell MI")
plt.title("Information")
plt.grid()
plt.legend()
plt.show()
plt.clf()


max_vars = data[16]
measured_vars = data[17]

max_vars_full = []
measured_vars_full = []
for i in range(len(max_vars[0])):
    max_vars_full.append(max_vars[0][i][0])
    measured_vars_full.append(measured_vars[0][i][0])

plt.plot(cell_stress_array, measured_vars_full,marker = 'o', label = "Adaptive Receptor SI")
plt.plot(cell_stress_array, max_vars_full,marker = 'o', label = "Equal Receptor SI")
plt.xlabel("Cell Stress")
plt.ylabel("Cell SI")
plt.title("Subjective Information")
plt.grid()
plt.legend()
plt.show()
plt.clf()