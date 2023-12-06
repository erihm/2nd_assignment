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


for set in enzymes:
    inputRules = list(set)

    try: 
        # This part is from dg.py
        pShortest = lambda d: all(a.vLabelCount("C") <= 8 for a in d.right) # WHY

        strat = (
            addSubset(ribuloseP, water)
            >> rightPredicate[pShortest]( # put pShorter, pEvenShorter or pShortest here, instead of p (they are all equivalent)
                    repeat(inputRules)
            )
        )
        dg = DG(graphDatabase=inputGraphs)
        dg.build().execute(strat)

        # This is from pathway.py
        flow = Flow(dg)
        flow.addSource(ribuloseP)
        flow.addSource(water)
        # flow.addConstraint(inFlow[ribuloseP] == 6)
        # flow.addConstraint(inFlow[water] == 1)
        flow.addSink(fructoseP)
        flow.addConstraint(outFlow[fructoseP] >= 1)
        flow.objectiveFunction = -outFlow[fructoseP]
        for v in dg.vertices: flow.addSink(v.graph)
        flow.findSolutions(verbosity=0)
        print("For this set of enzymes: ", inputRules)
        flow.solutions.list()
        good_enzymes.append(inputRules)

        print("For this set of enzymes we obtained %s solution/s" % len(flow.solutions))
        
        postSection("Enzymes: %s " % inputRules)
        flowPrinter = FlowPrinter()
        flowPrinter.printUnfiltered = False
        flow.solutions.print()
        
        del dg
        del flow
        del inputRules
        del strat
        
    except:
        bad_enzymes.append(inputRules)
        continue
    
for x in good_enzymes:
    print(x)
    
# First solution
# s0 = flow.solutions[0]
# Print info about the solution
# s0.list()
# Flow for this solution
# s0_flow = s0.flow
# s0_dg = s0.flow.dg