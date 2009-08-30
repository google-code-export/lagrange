import sys
import random
# the lagrange folder must be in the current folder or in the
# PYTHONPATH environment variable
import lagrange

#outfile = sys.stdout
# To print output to a file, uncomment the following line by removing
# the '#':
#outfile = file("results.txt", "a")


# Read matrix of geographic range data at leaf nodes on the phylogeny.
# The data file is formatted as follows: the first line contains the
# names of areas (short, preferably single-letter words), and
# following lines contain the taxon (single-word name) and its
# presence/absence in those areas, indicated by 1's and 0's.  Any
# amount of whitespace can be used to separate columns.  View the file
# "psychotria.matrix" for an example.
datafile = file("CAPR.matrix")
labels, data = lagrange.input.parse_matrix(datafile)

# Create a DEC RateModel for the set of areas parsed from the range
# data.  In this example, rates do not vary through time.
nareas = len(labels)
model = lagrange.RateModelGE(nareas, labels)

reps = 1000
ntrees = 27002
# The newick-formatted phylogeny with taxon names matching the
# geographic range data
for x in range(reps):
	print x
	handle = open("../STORE/all.tre", "rU")
	newicktree=""
	which = random.randint(0, ntrees)
	count = 0
	for i in handle:
		if count == which:
			newicktree=i
			break
		count += 1
	handle.close()
		
	tree = lagrange.Tree(newicktree)
	tree.set_default_model(model)
	tree.set_tip_conditionals(data)
	
	# Constraints
	# fix root node to Kauai only
	#tree.root.excluded_dists = [ dist for dist in model.dists if dist != (1,0,0,0) ]
	
	# Next, estimate global rates of dispersal and extinction
	mainoutfile = file("1/results.txt", "a")
	d, e = lagrange.output.optimize_dispersal_extinction(mainoutfile, tree, model)
	
	# Use global rates to evaluate ancestral range scenarios at internal nodes
	#lagrange.output.ancsplits(outfile, tree, model, d, e)
	
	#BASE
	outfile = file("1/1.txt", "a")
	names = ["Heptacodium_miconioides","Lonicera_japonica"]
	lagrange.output.ancsplits_with_names(outfile, tree, model, d, e, names)
	outfile.close()
	
	#SYMP LEY
	outfile = file("1/2.txt", "a")
	names = ["Symphoricarpos_albus","Leycesteria_formosa"]
	lagrange.output.ancsplits_with_names(outfile, tree, model, d, e, names)
	outfile.close()
	
	#SYMP
	outfile = file("1/3.txt", "a")
	names = ["Symphoricarpos_albus","Symphoricarpos_hesperius", "Symphoricarpos_occidentalis"]
	lagrange.output.ancsplits_with_names(outfile, tree, model, d, e, names)
	outfile.close()
	
	#TRIO
	outfile = file("1/4.txt", "a")
	names = ["Triosteum_angustifolium","Triosteum_sp", "Triosteum_perfoliatum"]
	lagrange.output.ancsplits_with_names(outfile, tree, model, d, e, names)
	outfile.close()
	
	#LONI CAPR
	outfile = file("1/5.txt", "a")
	names = ["Lonicera_arizonica","Lonicera_etrusca", "Lonicera_pilosa"]
	lagrange.output.ancsplits_with_names(outfile, tree, model, d, e, names)
	outfile.close()
	
	#LONI
	outfile = file("1/6.txt", "a")
	names = ["Lonicera_xylosteum","Lonicera_arizonica"]
	lagrange.output.ancsplits_with_names(outfile, tree, model, d, e, names)
	outfile.close()
	
	#LONI LONI
	outfile = file("1/7.txt", "a")
	names = ["Lonicera_xylosteum","Lonicera_hispida", "Lonicera_orespia"]
	lagrange.output.ancsplits_with_names(outfile, tree, model, d, e, names)
	outfile.close()
	
	#SYMP LONI
	outfile = file("1/8.txt", "a")
	names = ["Symphoricarpos_albus","Lonicera_japonica"]
	lagrange.output.ancsplits_with_names(outfile, tree, model, d, e, names)
	outfile.close()
	
	#TRIO LONI
	outfile = file("1/9.txt", "a")
	names = ["Triosteum_angustifolium","Lonicera_japonica"]
	lagrange.output.ancsplits_with_names(outfile, tree, model, d, e, names)
	outfile.close()
	
	# Run the analysis.  First, print the tree
	#lagrange.output.ascii_tree(outfile, tree, model, data)
	treeoutfile = open("1/sampledtrees.txt", "a")
	treeoutfile.write(lagrange.newick.to_string(tree.root)+";\n")
	treeoutfile.close()
	
