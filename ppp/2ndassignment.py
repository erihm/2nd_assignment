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


# Power set of enzymes
enzymes = powerset(inputRules)[1:]

good_enzymes = []
solutions = []

for subset in enzymes:
    
    try: # Here we want to avoid the error "Can not add null vertex as sink"
    
        pShortest = lambda d: all(a.vLabelCount("C") <= 8 for a in d.right) # The likelihood of a big molecule undergoing a reaction is very high

        strat = (
            addSubset(ribuloseP, water)
            >> rightPredicate[pShortest]( # put pShorter, pEvenShorter or pShortest here, instead of p (they are all equivalent)
                    repeat(subset)
            )
        )
        dg = DG(graphDatabase=inputGraphs)
        dg.build().execute(strat)
            
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
    flow.objectiveFunction = -outFlow[fructoseP] 
    for v in dg.vertices: flow.addSink(v)
    flow.findSolutions(verbosity=0) 
    
    if (len(flow.solutions) > 0):
        
        good_enzymes.append(set(subset))
        post.summarySection("Enzymes: %s " % list(subset))
        flowPrinter = FlowPrinter()
        flowPrinter.printUnfiltered = False
        flow.solutions.print()

        # save data for output and quality measure
        output = flow.solutions[0].eval(outFlow[fructoseP]) # outflow of fructose
        waste = flow.solutions[0].eval(outFlow) - output # outflow of waste = everything but fructose
        reactions = flow.solutions[0].eval(edge) # number of reactions
        solutions.append([subset, output, waste, round(waste / output, 4), reactions])
    
for x in good_enzymes:
    print(x)
    
draw_power_set_lattice(good_enzymes)

# generate table with results in summary
latex = (
#"\\newpage\n" +
"\\section{Results}\n" +
"\\begin{center}\n" +
"\\begin{longtable}{ |c|c|c|c|c|c|c| }\n" +
"\\hline\n" +
"id & \#enzymes & enzymes & \#fructose & \#waste & \#w\\slash\#f & \#reactions \\\\\n" +
"\\hline\n")

j = 0
for sol in solutions[::-1]:
    i = 0
    for en in sol[0]:
        if i < len(sol[0]) - 1:
            latex += " & " + " & " + str(en) + " & & & & \\\\\n"
        else:
            latex += str(j) + " & " + str(len(sol[0])) + " & " + str(en)
        i += 1
    latex += " & " + str(sol[1]) + " & " + str(sol[2]) + " & " + str(sol[3]) + " & " + str(sol[4]) + "\\\\\n \\hline\n"
    j += 1

latex += ("\\end{longtable}\n" +
"\\end{center}\n")

post.summaryRaw(latex)
