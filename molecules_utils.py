

def checkAtoms(molecules):
	for cod in molecules:
		molecules[cod].checkAtoms()


def createLJPairs(molecules):
	for cod in molecules:
		molecules[cod].createLJPairs()


def setCG(molecules, eem):
	for cod in molecules:
		eem.setCG(molecules[cod])


def computeEEM(eem, molecules, atomTypes):
	for cod in molecules:
		eem.computeEEM(molecules[cod], atomTypes)


