
import random, math
import scipy
import nchoosem
from n_node import *
from area import *

ran = random.Random()

class BioGeoTree:
    def __init__(self, extantstop = 0, timestop = 0, birthrate = 0.4, disprate = 0.02,
                 deathrate = 0.02, DEBUG = False, areaspecificextinction = False, areaextrate = [],
                 speciationlinked = False, areas = []):
        self.failures = 0
        self.maxfail = 10000
        self.areas = areas
        self.nareas = len(self.areas)
        self.dists = self.setup_dists()
        self.disprate = disprate # birth rate in an area
        self.birthrate = birthrate # speciation rate
        self.deathrate = deathrate # in an area, which translates to
        self.extantstop = extantstop
        self.timestop = timestop
        self.DEBUG = DEBUG
        self.areaspecificextinction = areaspecificextinction
        self.speciationlinked = speciationlinked
        self.printnumCH = False
        self.accdisp = []
        self.accext = []

    def setup_dists(self):
        # results in somthing like 
        # [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)]
        # for three areas
        dists = [ dist for dist in \
                           nchoosem.iterate_all_bv(self.nareas) ]
        return dists

    def make_tree(self, show_dead = False):
        self.accdisp = []
        self.accext = []
        self.numofchanges = 0
        self.numdisp = 0
        self.numext = 0
        self.numapext = 0
        self.numapdisp = 0
        self.setup_parameters()
        self.root = InternalNode(isroot = 1)
        ar = []
        for i in range(len(self.areas)):
            ar.append(0)
        ar[ran.randint(0,len(self.areas)-1)] = 1
        if self.DEBUG == True:
            for i in range(len(self.areas)):
                print ar[i]
        self.SETAREAS[self.root] = ar
        self.root.label = listToString(self.SETAREAS[self.root])
        self.CURAREAS[self.root] = ar # copy?
        self.BIRTHTIME[self.root] = self.currenttime
        self.extantNodes.append(self.root)
        while self.checkStopConditions():
            #try{
            dt = self.timeToNextSPEvent();
            at = self.timeToNextBioGeoEvent();
            if dt < at:
                if self.DEBUG == True:
                    print "spevent"
                self.currenttime += dt
                self.speciationEvent()
            # else if speciation event is before biogeo event
            else: 
                if self.DEBUG == True:
                    print "biogeoevent"
                self.currenttime += at
                self.biogeoEvent()
            for i in range(len(self.extantNodes)):
                if self.CURAREAS[self.extantNodes[i]].count(1) < 1:
                    self.killNode(self.extantNodes[i])
            if len(self.extantNodes) < 1:
                self.failures += 1
                if self.DEBUG == True:
                    print "failures: "+str(self.failures)
                self.setup_parameters()
                self.root = Node(label="root")
                ar = []
                for i in range(len(self.areas)):
                    ar.append(0)
                ar[ran.randint(0,len(self.areas)-1)] = 1
                if self.DEBUG == True:
                    for i in range(len(self.areas)):
                        print ar[i]
                self.SETAREAS[self.root] = ar
                self.root.label = listToString(self.SETAREAS[self.root])
                self.CURAREAS[self.root] = ar # copy?
                self.BIRTHTIME[self.root] = self.currenttime
                self.extantNodes.append(self.root)
        """
            set internal branch lengths
        """
        tempExtantNodes = self.extantNodes[:]
        for i in range(len(tempExtantNodes)):
            try:
                self.DEATHTIME[tempExtantNodes[i]]
            except KeyError:
                self.killNode(tempExtantNodes[i])
                tempExtantNodes[i].istip = 1
        if show_dead == False:
            self.deleteDeadNodes()
        self.root.length = 0
        totallength = 0 
        for i in self.root.leaves():
            totallength += i.length
        self.accdisp = len(self.accdisp)/totallength
        self.accext = len(self.accext)/totallength
        return self.root

    def setup_parameters(self):
        self.numofchanges = 0
        self.numdisp = []
        self.numext = []
        self.numapext = []
        self.numapdisp= []
        self.currenttime = 0.0
        self.relativedispersalrate = self.disprate / (self.disprate+self.deathrate)
        self.sumrate = self.disprate + self.deathrate
        self.extantNodes = []
        self.SETAREAS = {}
        self.CURAREAS = {}
        self.BIRTHTIME = {}
        self.DEATHTIME = {}

    def checkStopConditions(self):
        keepgoing = True
        if self.extantstop > 0:
            if len(self.extantNodes) >= self.extantstop:
                self.currenttime += self.timeToNextSPEvent()
                keepgoing = False;
        if self.timestop > 0.0:
            if self.currenttime >= self.timestop:
                self.currenttime = self.timestop
                keepgoing = False
        return keepgoing
    
    def speciationEvent(self):
        # pick a random lineage
        extant = self.extantNodes[ran.randint(0,len(self.extantNodes)-1)]
        self.nodeBirth(extant) # real speciation
            
    def biogeoEvent(self):
        # pick a random lineage
        extant = self.extantNodes[ran.randint(0,len(self.extantNodes)-1)]
        if self.eventIsDispersal():
            if self.DEBUG == True:
                    print "dispersal"
            self.nodeDispersal(extant)
        else:
            if self.DEBUG == True:
                    print "extinction"
            self.nodeExtinction(extant)

    def eventIsDispersal(self):
        x = ran.random()
        if x < self.relativedispersalrate:
            return True
        else:
            return False

    def iter_dist_splits(self, dist):
        assert dist in self.dists
        if sum(dist) == 1:
            yield (dist, dist)
        else:
            for i in scipy.nonzero(dist)[0]:
                x = scipy.zeros((len(dist),), dtype="i")
                x[i] = 1
                x = tuple(x)
                if x in self.dists:
                    yield (x, dist)
                    yield (dist, x)
                    y = tuple(scipy.array(scipy.logical_xor(dist, x),
                                          dtype="i"))
                    if y in self.dists:
                        yield (x, y)
                        if sum(y) > 1:
                            yield (y, x)
                        
    def nodeBirth(self,node):
        left = Node()
        right = Node()
        node.add_child(left)
        node.add_child(right)
        self.BIRTHTIME[left] = self.currenttime
        self.BIRTHTIME[right] = self.currenttime
        self.killNode(node)
        self.extantNodes.append(left)
        self.extantNodes.append(right)
        node.label = listToString(self.CURAREAS[node])
        splits = [ s for s in self.iter_dist_splits(tuple(self.CURAREAS[node])) ]
        nsplits = len(splits)
        x = ran.randint(0,nsplits-1)
        self.SETAREAS[left] = list(splits[x][0][:])
        self.SETAREAS[right] = list(splits[x][1][:])
        self.CURAREAS[left] = list(splits[x][0][:])
        self.CURAREAS[right] = list(splits[x][1][:])
        left.label = listToString(self.SETAREAS[left])
        right.label = listToString(self.SETAREAS[right])
    """
        node birth into specific area for speciation linked with dispersal
    """
    def nodeBirth2(self,node,sparea):
        left = Node()
        right = Node()
        node.add_child(left)
        node.add_child(right)
        self.BIRTHTIME[left] = self.currenttime
        self.BIRTHTIME[right] = self.currenttime
        self.killNode(node)
        self.extantNodes.append(left)
        self.extantNodes.append(right)
        node.label = listToString(self.CURAREAS[node])
        x = []
        for i in range(self.nareas):
            x[i] = 0;
            if i == sparea:
                x[i] = 1
        self.SETAREAS[left] = self.CURAREAS[node][:]
        self.SETAREAS[right] = x[:]
        self.CURAREAS[left] = self.CURAREAS[node][:]
        self.CURAREAS[right] = x[:]
        left.label = listToString(self.SETAREAS[left])
        right.label = listToString(self.SETAREAS[right])

    def nodeExtinction(self, node):
        ar = self.CURAREAS[node]
        x = 0
        if self.areaspecificextinction == False:
             x = ran.randint(0,ar.count(1)-1)
        else:
            tx = ran.random()
            start = 0
            end = 0
            for i in range(self.nareas):
                end = start + self.areas[i].getLocalExtinctionRate()
                if tx > start:
                    if tx <= end:
                        x = i
                        break
                start = end
        cur = 0
        for i in range(len(ar)):
            if ar[i] == 1:
                if x == cur:
                    ar[i] = 0
                    self.accext.append(node)
                    break
            if ar[i]==1:
                cur += 1
        if ar.count(1)<1:
            self.killNode(node)
        self.CURAREAS[node] = ar

    def nodeDispersal(self, node):
        ar = self.CURAREAS[node]
        if ar.count(1) != self.nareas:
            x = ran.randint(0,ar.count(1)-1)
            cur = 0
            sa = None
            sarea = 0
            for i in range(len(ar)):
                if ar[i] == 1:
                    if x == cur:
                        sa = self.areas[i]
                        sarea = i
                    cur += 1
            # attempt dispersal
            rt = ran.random()
            start = 0
            for i in range(len(ar)):
                if self.areas[i]!=sa and ar[i] == 0:
                    end = start + sa.getLocalDispersalRate(self.areas[i])
                    if rt > start and rt <= end and sa.getLocalDispersalRate(self.areas[i])!= 0:
                        if self.speciationlinked == False:
                            ar[i] = 1
                        self.accdisp.append(node)
                        if self.speciationlinked == True:
                            self.nodeBirth2(node,i)
                        break
                    start = end
            if self.speciationlinked == False:
                self.CURAREAS[node] = ar

    def killNode(self,node):
        self.DEATHTIME[node] = self.currenttime
        bl = self.DEATHTIME[node] - self.BIRTHTIME[node]
        node.length = bl
        node.label = listToString(self.CURAREAS[node])
        self.extantNodes.remove(node)

    def timeToNextSPEvent(self):
        return (-math.log(ran.random()))/ ( (len(self.extantNodes)) * self.birthrate)

    def timeToNextBioGeoEvent(self):
        return (-math.log(ran.random()))/ ( (len(self.extantNodes)) * self.sumrate)

    def deleteDeadNodes(self):
        kill = []
        setDistanceToTip(self.root)
        leaves = self.root.leaves()
        for i in range(len(leaves)):
            if getDistanceFromTip(leaves[i]) != self.root.distanceToTip:
                #print str(getDistanceFromTip(leaves[i]))+" "+str(self.root.distanceToTip)
                kill.append(leaves[i])
        for i in range(len(kill)):
            for j in self.accdisp:
                if j == kill[i]:
                    self.accdisp.remove(j)
            for j in self.accext:
                if j == kill[i]:
                    self.accext.remove(j)
            self.deleteANode(kill[i])

    def deleteANode(self, node):
        tparent = node.parent
        if tparent != self.root:
            child = None
            for i in range(len(tparent.children())):
                if tparent.children()[i] != node:
                    child = tparent.children()[i]
            pparent = tparent.parent
            tparent.children().remove(node)
            tparent.children().remove(child)
            pparent.children().remove(tparent)
            pparent.add_child(child)
            child.parent = pparent
            child.length = child.length+tparent.length
        else:
            child = None
            for i in range(len(tparent.children())):
                if tparent.children()[i] != node:
                    child = tparent.children()[i]
            tparent.children().remove(node)
            child.parent = None    
            self.root = child

def setDistanceToTip(root):
    for leaf in root.leaves():
        curh = 0.0
        leaf.distanceToTip = curh
        while leaf != None:
            curh += leaf.length
            if leaf.distanceToTip<curh:
                leaf.distanceToTip = curh
            leaf = leaf.parent

def listToString(inlist):
    ret = ""
    for i in inlist:
        ret += str(i)
    return ret

if __name__ == "__main__":
    from tree_reader_n import *
    import tree_printer_n
    ar = []
    ar1 = Area(0)
    ar2 = Area(1)
    ar3 = Area(2)
    ar1.dispersalMap[ar2]=0.5
    ar1.dispersalMap[ar3]=0.5
    ar2.dispersalMap[ar1]=0.5
    ar2.dispersalMap[ar3]=0.5
    ar3.dispersalMap[ar1]=0.5
    ar3.dispersalMap[ar2]=0.5
    ar.append(ar1)
    ar.append(ar2)
    ar.append(ar3)
    bgt = BioGeoTree(DEBUG = False,disprate = 0.02, deathrate = 0.02,extantstop = 20,birthrate = 0.08,areas = ar)
    root = bgt.make_tree(False);
    #for(int j=0;j<tree.getExternalNodeCount();j++){
    #    System.out.print((j+1)+"\t"+tree.getExternalNode(j).getName()+"\n");
    #    tree.getExternalNode(j).setName(String.valueOf(j+1));
    #}
    #System.out.println(tree.getRoot().getNewick(true)+"\n");
    print to_string(root)
    print tree_printer_n.render(root)
