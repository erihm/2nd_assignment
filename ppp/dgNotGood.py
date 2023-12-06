include("grammar.py")

strat = (addSubset(ribuloseP, water)
        >> repeat(inputRules)
        )
dg = DG(graphDatabase=inputGraphs)
dg.build().execute(strat)
