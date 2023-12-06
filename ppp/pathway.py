flow = Flow(dg)
flow.addSource(ribuloseP)
flow.addSource(water)
flow.addConstraint(inFlow[ribuloseP] == 6)
flow.addConstraint(inFlow[water] == 1)
flow.addSink(fructoseP)
flow.addConstraint(outFlow[fructoseP] >= 1)
flow.objectiveFunction = -outFlow[fructoseP]
for v in dg.vertices: flow.addSink(v.graph)
flow.findSolutions(
	#verbosity=2
)
flow.solutions.list()
