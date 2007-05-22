
class Area:
    def __init__(self,index):
        self.index = index
        #add con
        #add ext func
        self.dispersalMap = {}
        self.localextinctionrate = 1.0

    def setLocalDispersalRate(self, area, rate):
        dispersalMap.put(area, rate);
    
    def getLocalDispersalRate(self,area):
        return self.dispersalMap[area]

