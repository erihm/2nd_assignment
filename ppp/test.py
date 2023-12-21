from mod import *
from grammar import *
from itertools import chain, combinations
import graphviz

# Lattice

# Define a function to generate the powerset of an iterable.
def powerset(iterable):
    s = list(iterable)  # Convert the input iterable to a list.
    return list(chain.from_iterable(combinations(s, r) for r in range(len(s) + 1)))

def power_set(s):
    power_set_list = []
    for i in range(len(s) + 1):
        power_set_list.extend(combinations(s, i))
    return [set(subset) for subset in power_set_list]

def format_label(subset):
    return ',\n'.join(map(str, subset))

def draw_power_set_lattice(subsets_to_plot):
    graph = graphviz.Digraph('PowerSetLattice', format='png', graph_attr={'rankdir': 'BT'})

    for i, subset in enumerate(subsets_to_plot):
        graph.node(str(i), label=format_label(subset))

    for i, subset in enumerate(subsets_to_plot):
        for j, other_subset in enumerate(subsets_to_plot):
            if isinstance(subset, set) and isinstance(other_subset, set) and subset.issubset(other_subset) and i != j and len(subset) + 1 == len(other_subset):  # Check for immediate superset
                graph.edge(str(i), str(j))

    graph.render(filename='power_set_lattice', format='png', cleanup=True)


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

solutions = []

# power_set = enzymes[25]
# for u in range(1):
for power_set in enzymes:
    
    try: 
    
        pShortest = lambda d: all(a.vLabelCount("C") <= 8 for a in d.right) # The likelihood of a big molecule undergoing a reaction is very high

        strat = (
            addSubset(ribuloseP, water)
            >> rightPredicate[pShortest]( # put pShorter, pEvenShorter or pShortest here, instead of p (they are all equivalent)
                    repeat(power_set)
            )
        )
        dg = DG(graphDatabase=inputGraphs)
        dg.build().execute(strat)
            
        # This is from pathway.py
        flow = Flow(dg)
        flow.addSource(ribuloseP)
        flow.addSource(water)
        flow.addConstraint(inFlow[ribuloseP] == 60) # ideal would be to get 50 of fructose
        # for e in dg.edges:
        #     for rule in e.rules:
        #         if rule not in power_set:
        #             flow.addConstraint(isEdgeUsed[e] == 0)
        # flow.addConstraint(inFlow[water] == 1)
        
        flow.addSink(fructoseP)
        
    except libpymod.LogicError as e:
        print(f"Error: {e}")
        continue
    
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
        good_enzymes.append(set(power_set))
        post.summarySection("Enzymes: %s " % list(power_set))
        flowPrinter = FlowPrinter()
        flowPrinter.printUnfiltered = False
        flow.solutions.print()

        output = flow.solutions[0].eval(outFlow[fructoseP])
        waste = flow.solutions[0].eval(outFlow) - output
        reactions = flow.solutions[0].eval(isEdgeUsed)
        solutions.append([power_set, output, waste, reactions])
    
for x in good_enzymes:
    print(x)
    
draw_power_set_lattice(good_enzymes)

latex = (
"\\newpage\n" +
"\\section{Result}\n")

latex += ("\\begin{center}\n" +
"\\begin{longtable}{ |c|c|c|c|c| }\n" +
"\\hline\n" +
"\#enzymes & enzymes & \#Fructose-6-Phosphate & \#waste & \#reactions \\\\\n" +
"\\hline\n")

for sol in solutions:
    i = 0
    for en in sol[0]:
        if i < len(sol[0]) - 1:
            latex += " & " + str(en) + " & & & \\\\\n"
        else:
            latex += str(len(sol[0])) + " & " +str(en)
        i += 1
    latex += " & " + str(sol[1]) + " & " + str(sol[2]) + " & " + str(sol[3]) + "\\\\\n \\hline\n"

latex += ("\\end{longtable}\n" +
"\\end{center}\n")

post.summaryRaw(latex)

# for e in s0.flow.dg.edges:
#     for r in e.rules:
#         print(r.name)

# First solution
# s0 = flow.solutions[0]
# s0.eval(outFlow[fructoseP])
# Print info about the solution
# s0.list()
# Flow for this solution
# s0_flow = s0.flow
# s0_dg = s0.flow.dg