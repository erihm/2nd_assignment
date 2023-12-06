postSection("Loaded Graphs")
for a in inputGraphs: a.print()
postSection("Loaded Rules")
for a in inputRules: a.print()

dg.print()
postSection("Product Graphs")
for v in dg.vertices: v.graph.print()
flowPrinter = FlowPrinter()
flowPrinter.printUnfiltered = False
flow.solutions.print()
