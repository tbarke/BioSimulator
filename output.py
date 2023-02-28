import configparser
import math

import outputAnalysis
import utils


class output(object):

    class PrimitiveOutput:
        def __init__(self, totalcellsFile = None, totalReceptorsFile = None, boundReceptorsFile = None, internalABFile = None, totalConcsFile = None, cellLocationsFile = None, cellMovementFile = None, enviornmentConcsFile = None):
            self.totalcellsFile = totalcellsFile
            self.totalReceptorsFile = totalReceptorsFile
            self.boundReceptorsFile = boundReceptorsFile
            self.internalABFile = internalABFile
            self.totalConcsFile = totalConcsFile
            self.cellLocationsFile = cellLocationsFile
            self.cellMovementFile = cellMovementFile
            self.enviornmentConcsFile = enviornmentConcsFile

    class CalcOutput:
        def __init__(self, totalMIFile = None, averageMIFile = None, compositeMIFile = None, weightedMIFile = None, averageGrowthFile = None, SIAsVarFile = None, SIAsEntropyFile = None):
            self.totalMIFile = totalMIFile
            self.averageMIFile = averageMIFile
            self.compositeMIFile = compositeMIFile
            self.weightedMIFile = weightedMIFile
            self.averageGrowthFile = averageGrowthFile
            self.SIAsVarFile = SIAsVarFile
            self.SIAsEntropyFile = SIAsEntropyFile

    class RunClacOut:
        def __init__(self, MItrad = None, MI2D2D = None, MI2D1D = None, MI2D1Din = None, growth = None):
            self.MItrad = MItrad
            self.MI2D2D = MI2D2D
            self.MI2D1D = MI2D1D
            self.MI2D1Din = MI2D1Din
            self.growth = growth

    def __init__(self, runName = None, totalcellsFile = None, totalReceptorsFile = None, boundReceptorsFile = None, internalABFile = None, totalConcsFile = None, cellLocationsFile = None, cellMovementFile = None, enviornmentConcsFile = None, totalMIFile = None, averageMIFile = None, compositeMIFile = None, weightedMIFile = None, averageGrowthFile = None, SIAsVarFile = None, SIAsEntropyFile = None):
        self.runName = runName
        self.fileName = ''
        self.PrimitiveOutput = self.PrimitiveOutput(totalcellsFile,totalReceptorsFile, boundReceptorsFile, internalABFile, totalConcsFile, cellLocationsFile, cellMovementFile, enviornmentConcsFile)
        self.CalcOutput = self.CalcOutput(totalMIFile, averageMIFile, compositeMIFile, weightedMIFile, averageGrowthFile, SIAsVarFile, SIAsEntropyFile)
        self.RunClacOut = self.RunClacOut()

    def calculateMeasures(self, config, bins):
        try:
            MITrad, MI2D2D, MImove, growth = outputAnalysis.CalcData(config, bins, True, True, True, True, self.PrimitiveOutput.totalConcsFile, self.PrimitiveOutput.cellMovementFile, self.PrimitiveOutput.boundReceptorsFile, self.PrimitiveOutput.totalcellsFile)
            #TODO fix below and in write (object should be lower case to use the constructor, but this will break below)
            self.RunClacOut.MItrad = MITrad
            self.RunClacOut.MI2D2D = MI2D2D
            self.RunClacOut.MI2D1D = MImove
            #TODO implement below
            #self.RunClacOut.MI2D1Din = MI2D1Din
            self.RunClacOut.growth = growth
        except Exception as e:
            print("Error: ", e)

    def read(self, filename, supressWarnings = True):
        def parse_string(input_str):
            import re
            input_str = input_str.strip()
            # Check for integer
            if re.match(r'^[-+]?\d+$', input_str):
                return int(input_str)

            # Check for float
            if re.match(r'^[-+]?(\d+(\.\d*)?|\.\d+)$', input_str):
                return float(input_str)

            # Check for list of integers or floats
            if re.match(r'^\[(?:\s*[-+]?(?:\d+(?:\.\d*)?|\.\d+)\s*,?)*\]$', input_str):
                values_str = re.findall(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)', input_str)
                values = [int(x) if x.isdigit() else float(x) for x in values_str]
                return values

            # Check for boolean
            if input_str.lower() == 'true':
                return True
            if input_str.lower() == 'false':
                return False

            # Default to returning the input string as a string
            return input_str

        config = configparser.ConfigParser()
        config.read(filename)

        paramOrg = self.__dict__
        primitive = (int, str, bool, float)
        for key1, value1 in paramOrg.items():
            if isinstance(value1, primitive):
                paramOrg[key1] = config.get('nosection', key1)
                continue
            paramnames = value1.__dict__
            for key2 in paramnames:
                try:
                    paramnames[key2] = parse_string(config.get(key1, key2))
                except configparser.NoOptionError:
                    if not supressWarnings:
                        print("Warning: could not load in: " + key1 + "." + key2 + ", Probably loading in from an old version config. Using default Value.")

    def write(self, path, absolute = False):
        config = configparser.ConfigParser()

        priDict = {}
        calcDict = {}
        runClacDict = {}
        paramOrg = self.__dict__
        primitive = (int, str, bool, float)
        for key1, value1 in paramOrg.items():
            if isinstance(value1, primitive):
                config['nosection'] = {key1: value1}
                continue
            paramnames = value1.__dict__
            for key2, value2 in paramnames.items():
                newValue2 = value2
                if not value2:
                    newValue2 = ''
                if key1 == "PrimitiveOutput":
                    priDict[key2] = newValue2
                elif key1 == "CalcOutput":
                    calcDict[key2] = newValue2
                elif key1 == "RunClacOut":
                    runClacDict[key2] = newValue2

        config['PrimitiveOutput'] = priDict
        config['CalcOutput'] = calcDict

        fullpath = path
        if not absolute:
            nameDir = self.runName.split('/')[0]
            fullpath = path + "/" + nameDir  + "/OutputProfiles/" + self.runName.replace("/","-") + ".cfg"
        with open(fullpath, 'w+') as configfile:
            config.write(configfile)

        return fullpath

