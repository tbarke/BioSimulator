import configparser
import math

class configuration(object):

    class ConcParams:
        def __init__(self, simParams):
            self.concProfile = "vonMises"
            self.constantConc = 10
            self.uphillHigh = 14
            self.uphillLow = 0
            self.diffCoeff = 5
            self.arp = "periodic"
            self.repeatFrequency = 0
            self.magnitude = 100000.0
            self.locationA = 25
            self.locationB = 75
            self.locationC = 50
            self.concDiffuse = "static"
            self.VonMisesVar = 10 * ((2*math.pi)/simParams.length)
            self.VonMisesMagnitude = 175
            self.VonMisesAOffset = -25 * (1/simParams.locationStep)
            self.VonMisesBOffset = 25 * (1/simParams.locationStep)
            self.VonMisesCOffset = 0 * (1 / simParams.locationStep)

    class SimParams:
        def __init__(self):
            self.length = 100
            self.locationStep = 0.1
            self.simTimeStep = 0.01
            self.simLength = 20
            self.presetA_time = []
            self.presetB_time = []
            self.presetA_loc = []
            self.presetB_loc = []
            self.presetSave = False
            self.presetBool = False
            self.uniformConcStart = False
            self.simulationConcStart = True
            self.lam = 0.4
            self.gamma = 0.1
            self.aknoght = 200

    class CellMetaStats:
        def __init__(self):
            self.stress = 0.01
            self.absorptionRate = 1/self.stress
            self.receptorConsumptionRate = 0.0
            self.survivalCost = 250
            self.velocityMultiplier = 1.0
            self.mutate = "non"
            self.decisiontype = "non"
            self.noise = 0.0
            self.receptorMode = "simulate"
            self.combinedPortion = 1
            self.dividePortion = 1
            self.adaptiveRatio = 0
            self.cellLocations = []
            self.fullDivide = True
            self.fullDie = True
            self.dissocociationConstant = 2.0

    class OutputFlags:
        def __init__(self):
            self.totalcells = True
            self.totalReceptors = True
            self.boundReceptors = True
            self.internalAB = True
            self.totalConcs = True
            self.cellLocations = True
            self.cellMovement = True
            self.concProfile = True

    class RunOutputFlags:
        def __init__(self):
            self.save = True
            self.compressSave = True
            self.totalMI = True
            self.averageMI = True
            self.compositeMI = True
            self.weightedMI = True
            self.averageGrowth = True
            self.SIAsVar = True
            self.SIAsEntropy = True

    class CellStats:
        def __init__(self):
            self.Arec = 200
            self.Brec = 200
            self.maxRec = 400
            self.Amol = 0
            self.Bmol = 0
            self.ATP = 0
            self.biomass = 5
            self.generation = 1
            self.distTrav = 0
            self.startID = 0

    class RunStats:
        def __init__(self):
            self.stressArray = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
            self.noiseArray = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
            self.enviornmentArray = ["vonMises", "constUpDown", "constant", "manaFromHeaven"]
            self.cellStrategiesArray = ["non", "measured", "counter", "adjusted", "drastic", "drastic2"]
            self.cellStrategiesArrayFlag = True
            self.enviornmentArrayFlag = True
            self.runStress = True
            self.runNoise = False
            self.runDetermine = False
            self.beginningRandInt = 1
            self.save = True

    def __init__(self):
        self.simParams = self.SimParams()
        self.concParams = self.ConcParams(self.simParams)
        self.cellMetaStats = self.CellMetaStats()
        self.outputFlags = self.OutputFlags()
        self.runOutputFlags = self.RunOutputFlags()
        self.runStats = self.RunStats()
        self.cellStats = self.CellStats()

    def convert_array(self, stri):
        if stri[0] != '[' or stri[len(stri)-1] != ']':
            raise Exception("\"" + stri + "\": is not an array")

        sub_stri = stri[1:len(stri)-1]
        str_arr = sub_stri.split(',')
        ret = []
        if sub_stri.strip() == "":
            return ret
        for s in str_arr:
            try:
                ret.append(float(s.strip()))
            except ValueError:
                ret.append(s)
        return ret

    def readConfig(self, filename):
        config = configparser.ConfigParser()
        config.read(filename)

        self.runStats.stressArray = self.convert_array(config.get('RunStats', 'stressArray'))
        self.runStats.noiseArray = self.convert_array(config.get('RunStats', 'noiseArray'))
        self.runStats.runStress = config.get('RunStats', 'runStress') == "True"
        self.runStats.runDetermine = config.get('RunStats', 'runNoise') == "True"
        self.runStats.cellStrategiesArrayFlag = config.get('RunStats', 'cellStrategiesArrayFlag') == "True"
        self.runStats.enviornmentArrayFlag = config.get('RunStats', 'enviornmentArrayFlag') == "True"
        self.runStats.stressArray = config.get('RunStats', 'runDetermine') == "True"
        self.runStats.beginningRandInt = int(config.get('RunStats', 'beginningRandInt'))
        self.runStats.save = config.get('RunStats', 'save')  == "True"
        self.runStats.enviornmentArray = self.convert_array(config.get('RunStats', 'enviornmentArray'))
        self.runStats.cellStrategiesArray = self.convert_array(config.get('RunStats', 'cellStrategiesArray'))

        self.cellStats.Arec = int(config.get('CellStats', 'Arec'))
        self.cellStats.Brec = int(config.get('CellStats', 'Brec'))
        self.cellStats.maxRec = int(config.get('CellStats', 'maxRec'))
        self.cellStats.Amol = int(config.get('CellStats', 'Amol'))
        self.cellStats.Bmol = int(config.get('CellStats', 'Bmol'))
        self.cellStats.ATP = int(config.get('CellStats', 'ATP'))
        self.cellStats.biomass = int(config.get('CellStats', 'biomass'))
        self.cellStats.generation = int(config.get('CellStats', 'generation'))
        self.cellStats.distTrav = int(config.get('CellStats', 'distTrav'))
        self.cellStats.startID = int(config.get('CellStats', 'startID'))

        self.runOutputFlags.save = config.get('RunOutputFlags', 'save')  == "True"
        self.runOutputFlags.compressSave = config.get('RunOutputFlags', 'compressSave')  == "True"
        self.runOutputFlags.totalMI = config.get('RunOutputFlags', 'totalMI')  == "True"
        self.runOutputFlags.averageMI = config.get('RunOutputFlags', 'averageMI')  == "True"
        self.runOutputFlags.compositeMI = config.get('RunOutputFlags', 'compositeMI')  == "True"
        self.runOutputFlags.weightedMI = config.get('RunOutputFlags', 'weightedMI')  == "True"
        self.runOutputFlags.averageGrowth = config.get('RunOutputFlags', 'averageGrowth')  == "True"
        self.runOutputFlags.SIAsVar = config.get('RunOutputFlags', 'SIAsVar')  == "True"
        self.runOutputFlags.SIAsEntropy = config.get('RunOutputFlags', 'SIAsEntropy')  == "True"

        self.outputFlags.totalcells = config.get('OutputFlags', 'totalcells') == "True"
        self.outputFlags.totalReceptors = config.get('OutputFlags', 'totalReceptors') == "True"
        self.outputFlags.boundReceptors = config.get('OutputFlags', 'boundReceptors') == "True"
        self.outputFlags.internalAB = config.get('OutputFlags', 'internalAB') == "True"
        self.outputFlags.totalConcs = config.get('OutputFlags', 'totalConcs') == "True"
        self.outputFlags.cellLocations = config.get('OutputFlags', 'cellLocations') == "True"
        self.outputFlags.cellMovement = config.get('OutputFlags', 'cellMovement') == "True"
        self.outputFlags.concProfile = config.get('OutputFlags', 'concProfile') == "True"

        self.cellMetaStats.stress = float(config.get('CellMetaStats', 'stress'))
        self.cellMetaStats.absorptionRate = float(config.get('CellMetaStats', 'absorptionRate'))
        self.cellMetaStats.receptorConsumptionRate = float(config.get('CellMetaStats', 'receptorConsumptionRate'))
        self.cellMetaStats.survivalCost = int(config.get('CellMetaStats', 'survivalCost'))
        self.cellMetaStats.velocityMultiplier = float(config.get('CellMetaStats', 'velocityMultiplier'))
        self.cellMetaStats.mutate = config.get('CellMetaStats', 'mutate')
        self.cellMetaStats.decisiontype = config.get('CellMetaStats', 'decisiontype')
        self.cellMetaStats.noise = float(config.get('CellMetaStats', 'noise'))
        self.cellMetaStats.dissocociationConstant = float(config.get('CellMetaStats', 'dissocociationConstant'))
        self.cellMetaStats.receptorMode = config.get('CellMetaStats', 'receptorMode')
        self.cellMetaStats.combinedPortion = int(config.get('CellMetaStats', 'combinedPortion'))
        self.cellMetaStats.dividePortion = int(config.get('CellMetaStats', 'dividePortion'))
        self.cellMetaStats.adaptiveRatio = int(config.get('CellMetaStats', 'adaptiveRatio'))
        self.cellMetaStats.cellLocations = self.convert_array(config.get('CellMetaStats', 'cellLocations'))
        self.cellMetaStats.fullDivide = config.get('CellMetaStats', 'fullDivide') == "True"
        self.cellMetaStats.fullDie = config.get('CellMetaStats', 'fullDie') == "True"

        self.simParams.length = int(config.get('SimParams', 'length'))
        self.simParams.locationStep = float(config.get('SimParams', 'locationStep'))
        self.simParams.simTimeStep = float(config.get('SimParams', 'simTimeStep'))
        self.simParams.simLength = float(config.get('SimParams', 'simLength'))
        self.simParams.presetA_time = self.convert_array(config.get('SimParams', 'presetA_time'))
        self.simParams.presetB_time = self.convert_array(config.get('SimParams', 'presetB_time'))
        self.simParams.presetA_loc = self.convert_array(config.get('SimParams', 'presetA_loc'))
        self.simParams.presetB_loc = self.convert_array(config.get('SimParams', 'presetB_loc'))
        self.simParams.presetBool = config.get('SimParams', 'presetBool') == "True"
        self.simParams.uniformConcStart = config.get('SimParams', 'uniformConcStart') == "True"
        self.simParams.simulationConcStart = config.get('SimParams', 'simulationConcStart') == "True"
        self.simParams.lam = float(config.get('SimParams', 'lam'))
        self.simParams.gamma = float(config.get('SimParams', 'gamma'))
        self.simParams.aknoght = int(config.get('SimParams', 'aknoght'))
        self.simParams.presetSave = config.get('SimParams', 'presetSave') == "True"

        self.concParams.diffCoeff = int(config.get('ConcParams', 'diffCoeff'))
        self.concParams.constantConc = int(config.get('ConcParams', 'constantConc'))
        self.concParams.uphillHigh = int(config.get('ConcParams', 'uphillHigh'))
        self.concParams.uphillLow = int(config.get('ConcParams', 'uphillLow'))
        self.concParams.concProfile = config.get('ConcParams', 'concProfile')
        self.concParams.arp = config.get('ConcParams', 'arp')
        self.concParams.repeatFrequency = float(config.get('ConcParams', 'repeatFrequency'))
        self.concParams.magnitude = float(config.get('ConcParams', 'magnitude'))
        self.concParams.locationA = int(config.get('ConcParams', 'locationA'))
        self.concParams.locationB = int(config.get('ConcParams', 'locationB'))
        self.concParams.locationC = int(config.get('ConcParams', 'locationC'))
        self.concParams.concDiffuse = config.get('ConcParams', 'concDiffuse')
        self.concParams.VonMisesVar = float(config.get('ConcParams', 'VonMisesVar'))
        self.concParams.VonMisesMagnitude = float(config.get('ConcParams', 'VonMisesMagnitude'))
        self.concParams.VonMisesAOffset = float(config.get('ConcParams', 'VonMisesAOffset'))
        self.concParams.VonMisesBOffset = float(config.get('ConcParams', 'VonMisesBOffset'))
        self.concParams.VonMisesCOffset = float(config.get('ConcParams', 'VonMisesCOffset'))

    def writeConfig(self, filename):
        config = configparser.ConfigParser()
        config['RunStats'] = {'stressArray': self.runStats.stressArray,
                              'noiseArray': self.runStats.noiseArray,
                              'enviornmentArray': self.runStats.enviornmentArray,
                              'cellStrategiesArray': self.runStats.cellStrategiesArray,
                              'runStress': self.runStats.runStress,
                              'cellStrategiesArrayFlag': self.runStats.cellStrategiesArrayFlag,
                              'enviornmentArrayFlag': self.runStats.enviornmentArrayFlag,
                              'runNoise': self.runStats.runDetermine,
                              'runDetermine': self.runStats.stressArray,
                              'beginningRandInt': self.runStats.beginningRandInt,
                              'save': self.runStats.save
        }

        config['CellStats'] = {'Arec': self.cellStats.Arec,
                              'Brec': self.cellStats.Brec,
                              'maxRec': self.cellStats.maxRec,
                              'Amol': self.cellStats.Amol,
                              'Bmol': self.cellStats.Bmol,
                              'ATP': self.cellStats.ATP,
                              'biomass': self.cellStats.biomass,
                                'generation': self.cellStats.generation,
                               'distTrav': self.cellStats.distTrav,
                               'startID': self.cellStats.startID
                              }

        config['RunOutputFlags'] = {'save': self.runOutputFlags.save,
                               'compressSave': self.runOutputFlags.compressSave,
                               'totalMI': self.runOutputFlags.totalMI,
                               'averageMI': self.runOutputFlags.averageMI,
                               'compositeMI': self.runOutputFlags.compositeMI,
                               'weightedMI': self.runOutputFlags.weightedMI,
                               'averageGrowth': self.runOutputFlags.averageGrowth,
                               'SIAsVar': self.runOutputFlags.SIAsVar,
                               'SIAsEntropy': self.runOutputFlags.SIAsEntropy
                               }

        config['OutputFlags'] = {'totalcells': self.outputFlags.totalcells,
                                    'totalReceptors': self.outputFlags.totalReceptors,
                                    'boundReceptors': self.outputFlags.boundReceptors,
                                    'internalAB': self.outputFlags.internalAB,
                                    'totalConcs': self.outputFlags.totalConcs,
                                    'cellLocations': self.outputFlags.cellLocations,
                                    'cellMovement': self.outputFlags.cellMovement,
                                    'concProfile': self.outputFlags.concProfile
                                    }

        config['CellMetaStats'] = {'stress': self.cellMetaStats.stress,
                                    'absorptionRate': self.cellMetaStats.absorptionRate,
                                    'receptorConsumptionRate': self.cellMetaStats.receptorConsumptionRate,
                                    'survivalCost': self.cellMetaStats.survivalCost,
                                    'velocityMultiplier': self.cellMetaStats.velocityMultiplier,
                                    'mutate': self.cellMetaStats.mutate,
                                    'decisiontype': self.cellMetaStats.decisiontype,
                                    'noise': self.cellMetaStats.noise,
                                    'receptorMode': self.cellMetaStats.receptorMode,
                                    'combinedPortion': self.cellMetaStats.combinedPortion,
                                    'dividePortion': self.cellMetaStats.dividePortion,
                                    'adaptiveRatio': self.cellMetaStats.adaptiveRatio,
                                    'cellLocations': self.cellMetaStats.cellLocations,
                                    'fullDivide': self.cellMetaStats.fullDivide,
                                    'fullDie': self.cellMetaStats.fullDie,
                                    'dissocociationConstant': self.cellMetaStats.dissocociationConstant
                                    }

        config['SimParams'] = {'length': self.simParams.length,
                                   'locationStep': self.simParams.locationStep,
                                   'simTimeStep': self.simParams.simTimeStep,
                                   'simLength': self.simParams.simLength,
                                   'presetA_time': self.simParams.presetA_time,
                                   'presetB_time': self.simParams.presetB_time,
                                   'presetA_loc': self.simParams.presetA_loc,
                                   'presetB_loc': self.simParams.presetB_loc,
                                   'presetBool': self.simParams.presetBool,
                                   'uniformConcStart': self.simParams.uniformConcStart,
                                   'simulationConcStart': self.simParams.simulationConcStart,
                                   'lam': self.simParams.lam,
                                   'gamma': self.simParams.gamma,
                                   'aknoght': self.simParams.aknoght,
                                   'presetSave': self.simParams.presetSave
                               }

        config['ConcParams'] = {'diffCoeff': self.concParams.diffCoeff,
                               'arp': self.concParams.arp,
                               'concProfile': self.concParams.concProfile,
                               'constantConc': self.concParams.constantConc,
                               'uphillHigh': self.concParams.uphillHigh,
                               'uphillLow': self.concParams.uphillLow,
                               'repeatFrequency': self.concParams.repeatFrequency,
                               'magnitude': self.concParams.magnitude,
                               'locationA': self.concParams.locationA,
                               'locationB': self.concParams.locationB,
                               'locationC': self.concParams.locationC,
                               'concDiffuse': self.concParams.concDiffuse,
                               'VonMisesVar': self.concParams.VonMisesVar,
                               'VonMisesMagnitude': self.concParams.VonMisesMagnitude,
                               'VonMisesAOffset': self.concParams.VonMisesAOffset,
                               'VonMisesBOffset': self.concParams.VonMisesBOffset,
                               'VonMisesCOffset': self.concParams.VonMisesCOffset
                               }
        with open(filename, 'w+') as configfile:
            config.write(configfile)