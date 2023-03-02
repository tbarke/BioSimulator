import random
import log
l = log.log()

class randomClass(object):


    def __init__(self, seed):
        self.currentSeed = seed
        self.allseeds = []

    def getNewInt(self, low, high):
        random.seed(self.currentSeed)
        ret = random.randint(low,high)
        return random.randint(low, high)
        random.seed(self.currentSeed)
        self.currentSeed = random.randint(1,10000) + 1
        return ret

    def appendAllSeeds(self, seed):
        self.allseeds.append(seed)

    def choices(self, arr, probs):
        random.seed(self.currentSeed)
        ret = random.choices(arr,probs)
        random.seed(self.currentSeed)
        self.currentSeed = random.randint(1,10000) + 1
        return ret

    def giveAllSeeds(self):
        return self.allseeds
