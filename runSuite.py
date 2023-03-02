import run
import math
import utils
import utils
import log
l = log.log()

#parameters[0] is SimParams
#parameters[1] is ConcParams
#parameters[2] is CellMetaStats
#parameters[3] is CellStats
#parameters[4] is Celllocations
#parameters[5] is EnviornmentParams
#parameters[6] is simrepeat

def simulate_noisy(cell_stress, number, cell_type, location, noise, parameters, bind_type):
    survivaly = []
    MIy = []
    entx = []
    entrxgiveny = []
    MI_vars_all = []
    entr = []
    entr_in = []
    parameters[2][0] = 1.0/cell_stress
    parameters[2][6] = noise
    parameters[2][7] = bind_type

    enviornment_stats = []

    if cell_type == "Measured":
        parameters[2][5] = "measured"
    else:
        parameters[2][5] = "non"
    for i in range(parameters[6] ):
        #print(str(i) + ": Repeat Number")
        timeRun = run.divisonMIRun(parameters[0], parameters[1], parameters[2], parameters[3], parameters[4], parameters[5], False)
        survivaly.append(timeRun[0])
        MIy.append(timeRun[1])
        entr.append(timeRun[2])
        entr_in.append(timeRun[3])
        entx.append(timeRun[4])
        entrxgiveny.append(timeRun[5])
        MI_vars = timeRun[7]
        E_MI_var = 0.0
        enviornment_stats.append(timeRun[8])
        for i in range(len(MI_vars)):
            E_MI_var = E_MI_var + MI_vars[i]
        if len(MI_vars) > 0:
            MI_vars_all.append(E_MI_var/len(MI_vars))

    new_survival = []
    #print(len(survivaly))
    for j in range(0,len(survivaly)):
        l.log(survivaly[j])
        new_survival.append([])
        for i in range(0,len(survivaly[j])-1):
            new_survival[j].append(math.log(survivaly[j][i+1]/survivaly[j][i],2))
            #print(math.log(survivaly[j][i+1]/survivaly[j][i],2))
        #print(len(new_survival[j]))
        new_survival[j].append(new_survival[j][len(new_survival[j]) - 1])
    #exit(-1)
    growth_sum = 0.0
    MI_sum = 0.0
    entr_sum = 0.0
    entr_input = 0.0
    count = 0.0
    l.log(len(MIy))
    for j in range(len(MIy)):
        l.log(len(MIy[j]))
        for i in range(len(MIy[j])):
            count = count + 1
            growth_sum = growth_sum + new_survival[j][i]
            MI_sum = MI_sum + MIy[j][i]
            entr_sum = entr_sum + entr[j][i]
            entr_input = entr_input + entr_in[j][i]

    count1 = 0.0
    for j in range(len(new_survival)):
        for i in range(len(new_survival[j])):
            count1 = count1 + 1
            growth_sum = growth_sum + new_survival[j][i]

    growth_sum = growth_sum/count1
    MI_sum = MI_sum /count
    entr_sum = entr_sum / count
    entr_input = entr_input / count
    return growth_sum, MI_sum, entr_sum, entr_input, MIy, survivaly, entx, entrxgiveny, MI_vars_all, enviornment_stats

#parameters[0] is SimParams
#parameters[1] is ConcParams
#parameters[2] is CellMetaStats
#parameters[3] is CellStats
#parameters[4] is Celllocations
#parameters[5] is EnviornmentParams
#parameters[6] is simrepeat

def simulate(cell_stress, number, cell_type, location, parameters):
    survivaly = []
    MIy = []
    entr = []
    entr_in = []
    parameters[2][0] = 1.0/cell_stress
    if cell_type == "Measured":
        parameters[2][5] = "measured"
    else:
        parameters[2][5] = "non"
    for i in range(parameters[6]):
        l.log(str(i) + ": Repeat Number")
        timeRun = run.divisonMIRun(parameters[0], parameters[1], parameters[2], parameters[3], parameters[4], parameters[5], False)
        survivaly.append(timeRun[0])
        MIy.append(timeRun[1])
        entr.append(timeRun[2])
        entr_in.append(timeRun[3])

    new_survival = []
    l.log(len(survivaly))
    for j in range(0,len(survivaly)):
        new_survival.append([])
        for i in range(0,len(survivaly[j])-1):
            new_survival[j].append(math.log(survivaly[j][i+1]/survivaly[j][i],2))
        l.log(len(new_survival[j]))
        new_survival[j].append(new_survival[j][len(new_survival[j]) - 1])

    growth_sum = 0.0
    MI_sum = 0.0
    entr_sum = 0.0
    entr_input = 0.0
    count = 0.0
    for j in range(len(new_survival)):
        for i in range(len(new_survival[j])):
            count = count + 1
            growth_sum = growth_sum + new_survival[j][i]
            MI_sum = MI_sum + MIy[j][i]
            entr_sum = entr_sum + entr[j][i]
            entr_input = entr_input + entr_in[j][i]

    growth_sum = growth_sum/count
    MI_sum = MI_sum /count
    entr_sum = entr_sum / count
    entr_input = entr_input / count
    return growth_sum, MI_sum, entr_sum, entr_input

#noise_array = [0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
#noise_array = [0.0,0.01, 0.02, 0.03, 0.04, 0.05,0.06, 0.07, 0.08, 0.09, 0.1,0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0]
#noise_array = [0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,2,3,4]
#noise_array = [0.0, 0.2, 0.3]
noise_array = [0.0]
#noise_array = [16, 32 ,64 , 128, 256, 400, 410, 420, 450, 475, 500, 575, 700, 850, 1200, 1400, 1800, 2300, 3000, 4000]
#noise_array = [20000, 30000, 50000]
#noise_array = [4096, 2048, 1024, 512, 256, 128, 64, 32, 16]
#noise_array = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
#noise_array = [0.0, 0.3, 0.7, 1]
#noise_array = [0]

def suite_noise_stress(noise_array, cell_stress_array, parameters, run_name):
    full_measured_growth = []
    full_max_growth = []
    full_max_MI = []
    full_measured_MI = []

    full_max_full_growth = []
    full_measured_full_growth = []
    full_max_full_MI = []
    full_measured_full_MI = []

    full_max_full_Hx = []
    full_max_full_Hxgiveny = []

    full_measured_full_hx = []
    full_measured_full_hxgiveny = []

    enviornment_stats = []
    for noise in noise_array:
        MI_vars_measured_full = []
        MI_vars_max_full = []

        measured_growth = []
        measured_MI = []
        measured_entr = []
        measured_entr_in = []
        max_growth = []
        max_MI = []
        max_entr = []
        max_entr_in = []

        measured_full_growth = []
        max_full_growth = []

        measured_full_MI = []
        max_full_MI = []

        entx_max = []
        entx_measured = []
        entxgiveny_max = []
        entxgiveny_measured = []

        MI_vars_measured = []
        MI_vars_max = []

        for i in range(len(cell_stress_array)):
            l.log("noise ARray: " + str(noise) + "cell stress array: " + str(cell_stress_array[i]))
            l.log("Measured Cell, Stress: " + str(cell_stress_array[i]) + " Percent Done: " + str(i/len(cell_stress_array)) + " ........................................................")
            array = simulate_noisy(cell_stress_array[i], i, "Measured", "10-27-21", noise, parameters, "gaussian")
            measured_growth.append(array[0])
            measured_MI.append(array[1])
            measured_entr.append(array[2])
            measured_entr_in.append(array[3])

            measured_full_growth.append(array[5])
            measured_full_MI.append(array[4])

            entx_measured.append(array[6])
            entxgiveny_measured.append(array[7])

            MI_vars_measured.append(array[8])

            enviornment_stats.append(array[9])
            l.log("Maximum Sensitivity: " + str(cell_stress_array[i]) + " Percent Done: " + str(i/len(cell_stress_array))+ " ........................................................")
            array = simulate_noisy(cell_stress_array[i], i, "Maximum Sensitivity", "10-27-21", noise, parameters, "gaussian")
            max_growth.append(array[0])
            max_MI.append(array[1])
            max_entr.append(array[2])
            max_entr_in.append(array[3])

            max_full_growth.append(array[5])
            max_full_MI.append(array[4])

            entx_max.append(array[6])
            entxgiveny_max.append(array[7])

            MI_vars_max.append(array[8])
            enviornment_stats.append(array[9])

        full_measured_growth.append(measured_growth)
        full_max_growth.append(max_growth)
        full_max_MI.append(max_MI)
        full_measured_MI.append(measured_MI)

        full_max_full_growth.append(max_full_growth)
        full_measured_full_growth.append(measured_full_growth)
        full_max_full_MI.append(max_full_MI)
        full_measured_full_MI.append(measured_full_MI)

        full_max_full_Hx.append(entx_max)
        full_max_full_Hxgiveny.append(entxgiveny_max)

        full_measured_full_hx.append(entx_measured)
        full_measured_full_hxgiveny.append(entxgiveny_measured)

        MI_vars_measured_full.append(MI_vars_measured)
        MI_vars_max_full.append(MI_vars_max)

    title_array = []
    data_array = []

    title_array.append("sim_parameters")
    data_array.append(parameters)

    title_array.append("noise_array")
    data_array.append(noise_array)

    title_array.append("cell_stress_array")
    data_array.append(cell_stress_array)

    title_array.append("meaaured_growth")
    data_array.append(full_measured_growth)

    title_array.append("max_growth")
    data_array.append(full_max_growth)

    title_array.append("meaaured_MI")
    data_array.append(full_measured_MI)

    title_array.append("max_MI")
    data_array.append(full_max_MI)

    title_array.append("max_HX")
    data_array.append(full_max_full_Hx)

    title_array.append("max_HXgivenY")
    data_array.append(full_max_full_Hxgiveny)

    title_array.append("measured_HX")
    data_array.append(full_measured_full_hx)

    title_array.append("measured_HXgivenY")
    data_array.append(full_measured_full_hxgiveny)

    title_array.append("full_max_full_growth")
    data_array.append(full_max_full_growth)

    title_array.append("full_measured_full_growth")
    data_array.append(full_measured_full_growth)

    title_array.append("full_max_full_MI")
    data_array.append(full_max_full_MI)

    title_array.append("full_measured_full_MI")
    data_array.append(full_measured_full_MI)

    title_array.append("full_max_full_vars")
    data_array.append(MI_vars_max_full)

    title_array.append("full_measured_full_vars")
    data_array.append(MI_vars_measured_full)

    title_array.append("envoiornment_stats")
    data_array.append(enviornment_stats)

    data_array.insert(0, title_array)

    l.log("Run Finished Exporting Data...")
    return utils.saveDataDate(run_name, data_array, True)