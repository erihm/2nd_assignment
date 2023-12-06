for s in flow.solutions:
	query = CausalityFlowQuery(flow.dg)
	printer = DGPrinter()
	printer.graphvizPrefix = 'layout = "dot";'
	printer.pushVertexVisible(lambda v: s.eval(vertex[v]) != 0)
	dag = query.findDAG(s)
	if dag:
		data = dag.getPrintData(False)
		dg.print(printer=printer, data=data)
	else:
		print(s, "has no DAG, find catalysts")
		for v in dg.vertices:
			a = v.graph
			if s.eval(vertex[a]) == 0: continue
			query = CausalityFlowQuery(flow.dg)
			query[a] = 1
			dag = query.findDAG(s)
			if dag:
				postSection("Catalyst %s" % a.name)
				data = dag.getPrintData(True)
				dg.print(printer=printer, data=data)
