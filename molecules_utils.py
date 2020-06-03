

def checkAtoms(molecules):
	for cod in molecules:
		molecules[cod].checkAtoms()


def setCG(molecules, eem):
	for cod in molecules:
		eem.setCG(molecules[cod])


def computeEEM(eem, molecules, atomTypes):
	for cod in molecules:
		eem.computeEEM(molecules[cod], atomTypes)


def createEffectivePrms(molecules, eem):
	for cod in molecules:
		mol = molecules[cod]
		mol.createEffectivePrms(eem)
