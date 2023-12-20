from mod import *
from grammar import *
from itertools import chain, combinations

# Lattice

# Define a function to generate the powerset of an iterable.
def powerset(iterable):
    s = list(iterable)  # Convert the input iterable to a list.
    return list(chain.from_iterable(combinations(s, r) for r in range(len(s) + 1)))

enzymes = powerset(inputRules)[1:]

good_enzymes = []
bad_enzymes = []

# For NOW:
# (E) Use the same dg for all the solutions (2-3)
# (E) Look at which errors we get
# (L) To access the values of the flow solution use eval() (3) Make the table solution, amount input, amount waste, amount output, num hyperedges in the solution
# Does the quality improve with more enzymes? Is the flow more complex w. more enzymes?
# (Optional): compare multiple solutions w. same objective value.
# (L) Visualize lattice with values (3)

# For LATER:
# Add an explanation for the subsets that worked
# Use another chemistry: knock?
# Design a system w. artificial rules for illustration
# "Waste" is the amount of carbons
# W. more enzymes the flow might be shorter.

# This part is from dg.py
pShortest = lambda d: all(a.vLabelCount("C") <= 8 for a in d.right) # The likelihood of a big molecule undergoing a reaction is very high

strat = (
    addSubset(ribuloseP, water)
    >> rightPredicate[pShortest]( # put pShorter, pEvenShorter or pShortest here, instead of p (they are all equivalent)
            repeat(inputRules)
    )
)
dg = DG(graphDatabase=inputGraphs)
dg.build().execute(strat)

for power_set in enzymes:
    
    # This is from pathway.py
    flow = Flow(dg)
    flow.addSource(ribuloseP)
    flow.addSource(water)
    flow.addConstraint(inFlow[ribuloseP] == 60) # ideal would be to get 50 of fructose
    for e in dg.edges:
        for rule in e.rules:
            if rule not in power_set:
                flow.addConstraint(isEdgeUsed[e] == 0)
    # flow.addConstraint(inFlow[water] == 1)
    flow.addSink(fructoseP)
    flow.addConstraint(outFlow[fructoseP] >= 1)
    flow.objectiveFunction = -outFlow[fructoseP] # dg.findVertex(fructoseP)
    for v in dg.vertices: flow.addSink(v)
    flow.findSolutions(verbosity=0) #, absGap=0, maxNumSolution=2**20) # check absGap
    
    # while True: 
    #   sols = flow.findSolution()
    #   if jfiohfaohaio:
    #   break
    
    if (len(flow.solutions) == 0):
        bad_enzymes.append(power_set)
        continue
    else:
        good_enzymes.append(power_set)
        post.summarySection("Enzymes: %s " % power_set)
        flowPrinter = FlowPrinter()
        flowPrinter.printUnfiltered = False
        flow.solutions.print()       
    
for x in good_enzymes:
    print(x)
    
# First solution
# s0 = flow.solutions[0]
# s0.eval(outFlow[fructoseP])
# Print info about the solution
# s0.list()
# Flow for this solution
# s0_flow = s0.flow
# s0_dg = s0.flow.dg