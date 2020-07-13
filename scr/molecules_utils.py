

def checkAtoms(molecules):
	for cod in molecules:
		molecules[cod].checkAtoms()


def computeEEM(eem, molecules, atomTypes):
	for cod in molecules:
		mol = molecules[cod]
		atoms = mol.atoms

		for prm in mol.parameters:
			prm.computeEEM(eem, atomTypes, atoms)


def computeCR(cr, molecules, atomTypes, matrix):
	for cod in molecules:
		mol = molecules[cod]

		for prm in mol.parameters:
			prm.computeCR(cr, atomTypes, matrix)


def createEffectivePrms(atomTypes, molecules, eem, matrix):
	for cod in molecules:
		mol = molecules[cod]
		eem.setCG(molecules[cod])
		mol.createEffectivePrms(atomTypes, eem, matrix)




