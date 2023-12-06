include("grammar.py")

dg = DG(graphDatabase=inputGraphs)
dg.build().execute(
	addSubset(ribuloseP, water)
	>> repeat[steps](inputRules)
)
