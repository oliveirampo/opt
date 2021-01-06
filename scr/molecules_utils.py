

def checkAtoms(molecules):
	for cod in molecules:
		molecules[cod].checkAtoms()


def computeEEM(eem, molecules, atomTypes):
	for cod in molecules:
		mol = molecules[cod]

		chargeGroups = mol.CGs
		for cg in chargeGroups:
			eem.computeEEM(cg, atomTypes, mol)


def computeCR(cr, molecules, atomTypes, matrix):
	for cod in molecules:
		mol = molecules[cod]

		for prm in mol.parameters:
			prm.computeCR(cr, atomTypes, matrix)


def createEffectivePrms(atomTypes, molecules, eem):
	for cod in molecules:
		mol = molecules[cod]

		eem.setCG(molecules[cod])
		mol.createLJPairs(atomTypes)


def printPrmValue(idx, cod, molecules):
	prm = molecules[cod].parameters[idx]
	print(prm.nam, prm.cur)

