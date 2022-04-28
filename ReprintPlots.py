"""""
ReprintPlots.py: reprints the plots if they were not originally formatted correctly
"""""

#import matplotlib.pyplot as plt
import utils

#titles are the files with the data from a single simulation run that needs to be reprinted
titles = []
titles.append("Output/Data/qual/Maximum Sensitivity Model Survivability over time0.txt")

finalData = []
for k in range(1):
    f = open(titles[k], "r")
    lines = f.readlines()
    for i in range(2, len(lines)):
        arrline = lines[i].split(", ")
        doubleArrLine = []
        for j in range(len(arrline)-1):
            doubleArrLine.append(float(arrline[j]))
        finalData.append(doubleArrLine)
    f.close()
print(finalData)

fianlTimeDivx = range(101)
titleNum = "0"
newLineTitles = []
newLineTitles.append("Maximum Sensitivity Model Survivability")
newLineTitles.append("Single Cell")
xlabel = "time step"
ylabel = "division probability"

figDIR = "Output/Plots/qual"
dataDIR = "Output/Data/qual"

YMax = 0.5
utils.outputDataError(finalData, fianlTimeDivx, xlabel, ylabel, newLineTitles, YMax, titleNum, figDIR, dataDIR, False)
