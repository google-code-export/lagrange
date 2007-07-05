import psyco
psyco.profile()
psyco.log()
import tree_reader_n
from tree_reader_n import *
import tree_printer_n
from area import *
from biogeotree import * 
import sys
from pprint import pprint
import model_optimize
import rate_matrix
from rate_model import *
from ancsplit import *
import branch_segment
import ancsplit
import scipy
from scipy import optimize
from marginal_range_reconstructor import *
LARGE = 10e5
PMAX = 10.0



def str_list_to_int_tuple(strlist):
    for i in range(len(strlist)):
        strlist[i] = int(strlist[i])
    tup = tuple(strlist)
    return tup

def comb(ds1, ds2):
    ret = ""
    for i in range(len(ds1)):
        if ds1[i] == '1' or ds2[i] == '1':
            ret += '1'
        else:
            ret += '0'
    return ret

def calculate_global_for_all_nodes_with_internal_node_counts(node, model, d, e, skip=[], _rv={}):
    global internal_correct
    global internal_bigger
    global internal_smaller
    global contained
    if (not node.istip) and (not node.label in skip):
        c1, c2 = node.children()
        calculate_global_for_all_nodes_with_internal_node_counts(c1, model, d, e)
        calculate_global_for_all_nodes_with_internal_node_counts(c2, model, d, e)
        print ", ".join([node.label, c1.label, c2.label])
        results = []
        #node.tree.clear_startdist()
        for disti, dist in model.enumerate_dists():
            S = set()
            for ancsplit in model.iter_ancsplits(dist):
                d1, d2 = ancsplit.descdists
                # calcuate likelihood of this split at this node
                lh = model_optimize.ancsplit_likelihood_de(node, ancsplit, model, d, e)
                
                ds1 = model.diststrings[model.dist2i[d1]]
                ds2 = model.diststrings[model.dist2i[d2]]
                root = zip(scipy.log(node.tree.root.dist_conditionals),
                           model.diststrings)
                root.sort(); root.reverse()
                for i, r in enumerate(root):
                    dc, ds = r
                    if dc < (root[0][0] - 2):
                        break
                root = root[:i]
                try:
                    results.append((-(math.log(lh)), ds1, ds2, root))
                except:
                    pass

        results.sort()
        v = []
        first = True
        got_it = False
        cont = False
        for opt, ds1, ds2, root in results:
            combd = comb(ds1,ds2)
            right = 0
            wrong = 0
            for i in range(len(node.label)):
                if node.label[i] == combd[i]:
                    right += 1
                else:
                    wrong += 1
            if first == True:
                if right == len(node.label) and wrong == 0:
                    internal_correct += 1
                    got_it = True
                elif combd.count('1') < node.label.count('1') and right > 1:
                    internal_smaller += 1
                elif combd.count('1') > node.label.count('1') and right > 1:
                    internal_bigger += 1
            if got_it == False and first == False and right == len(node.label) and wrong == 0 and cont == False:
                contained += 1
                cont = True
            if opt > (results[0][0] + 2):
                break
            v.append((opt, ds1, ds2))
            print "  -lnL %g, %s, %s" % (opt, ds1, ds2)
            first = False
        _rv[node.label] = v

    return _rv

def likelihood_de_multi(params, model, rec):
    d, e = params
    for p in (d, e):
        if (p < 0) or (p > PMAX):
            return LARGE
    model.setup_D(d)
    model.setup_E(e)
    model.setup_Q()
    try:
        lh = 1
        for i in rec:
            lh *= i.eval_likelihood()
        return -(math.log(lh))
    except:
        return LARGE

#rec start
reps = 20
dists = None
periods = [100]
model = RateModel(3, periods=periods, dists=dists)
print "disp\text\tavgaccdisp\tavgaccext\tdisp\text"

#make biogeotree
extant = 100
dispers = 0.2
death = 0.2
brate = 0.4
numberofsimul = 1
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
bgt = BioGeoTree(DEBUG = False,disprate = dispers, deathrate = death,extantstop = extant,birthrate = brate,areas = ar)
for rep in range(reps):
    #simulation
    data = []
    newickstr = []
    accdisp = []
    accext = []
    for i in range(numberofsimul):
        root = bgt.make_tree(False)
        while len(root.leaves()) != extant:
            root = bgt.make_tree(False)
        accdisp.append(bgt.accdisp)
        accext.append(bgt.accext)
        count = 0
        datat = {}
        for j in root.leaves():
            #print str(count+1)+"\t"+j.label
            datat[str(count+1)] = str_list_to_int_tuple(list(j.label))
            j.label = str(count+1)
            count += 1
        newickstrt = to_string(root)+";"
        #print newickstrt
        data.append(datat)
        newickstr.append(newickstrt)
    
    #reconstruction
    rec = []
    for i in range(numberofsimul):
        rect = MarginalRangeReconstructor(model, newickstr[i], periods)
        rect.set_default_model(model)
        rect.set_tip_conditionals(data[i])
        rec.append(rect)
    
    v = optimize.fmin_powell(
        likelihood_de_multi, [dispers,death], args=(model, rec),
        full_output=True,disp=0
        )
    print str(v[:2][0][0])+"\t"+str(v[:2][0][1])+"\t"+str(sum(accdisp) / len(accdisp))+"\t"+str(sum(accext) / len(accext))+"\t"+str(dispers)+"\t"+str(death)

#internal_correct = 0
#internal_bigger = 0
#internal_smaller = 0
#contained = 0
#calculate_global_for_all_nodes_with_internal_node_counts(rec.root, model, v[:2][0][0], v[:2][0][1])
#print "internal_correct "+ str(internal_correct)
#print "internal_bigger "+ str(internal_bigger)
#print "internal_smaller "+ str(internal_smaller)
#print "contained "+ str(contained)

"""private void print_results(){
        double int_cont = rr.get_tree().getInternalNodeCount();
        String print = "";
        if(internal_correct != 0)
            print = print + (internal_correct/int_cont)+"\t";
        else
            print = print + 0+"\t";
        if(internal_bigger != 0)
            print = print + (internal_bigger/int_cont)+"\t";
        else
            print = print + 0+"\t";
        if(internal_smaller != 0)
            print = print + (internal_smaller/int_cont)+"\t";
        else
            print = print + 0+"\t";
        if(contained != 0)
            print = print + (contained/int_cont)+"\t";
        else
            print = print + 0+"\t";
        System.out.println(this.currentnumber+"\t"+this.finish_time+"\t"+print);
        print = print +"\t"+this.opt_vag;
        print = print + "\t"+this.opt_vag2;
        print = print + "\t"+this.opt_vag3;
        print = print + "\t"+this.opt_vag4;
        print = print + "\t"+this.opt_vag5;
        print = print + "\t"+this.opt_vag6;
        print = print + "\t"+this.opt_ext;
        //print = print + "\t"+this.opt_ext2;
        //print = print + "\t"+this.opt_ext3;
        print = print +"\t"+this.apdisp;
        print = print + "\t"+this.apext;
        print = print +"\t"+this.redisp;
        print = print + "\t"+this.reext;
        print = print + "\t"+this.height;
        if(testSig == true){
            print = print + "\t"+this.isSig;
        }
        try {
            fw.write(print+"\n");
            fw.flush();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }"""

sys.exit()
