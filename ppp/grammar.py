from mod import *

ribuloseP = smiles("OCC(=O)C(O)C(O)COP(=O)(O)O", "Ribulose-5-Phosphate")
water = smiles("O", "H2O")
aldoKetoB = ruleGML("aldo_keto_backward.gml")
aldoKetoF = ruleGML("aldo_keto_forward.gml")
transKeto = ruleGML("transketolase.gml")
transAldo = ruleGML("transaldolase.gml")
aldolase = ruleGML("aldolase.gml")
phosphohydro = ruleGML("phosphohydrolase.gml")

fructoseP = smiles("OCC(=O)C(O)C(O)C(O)COP(=O)(O)O", "Fructose-6-Phosphate")
