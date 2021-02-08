

def checkAtoms(molecules):
	for cod in molecules:
		molecules[cod].checkAtoms()


def computeChargeDistribution(eem, molecules, atom_types, kap, lam):
	max_order = 2

	for cod in molecules:
		mol = molecules[cod]

		chargeGroups = mol.CGs
		for cg in chargeGroups:
			eem.compute(max_order, mol, cg, atom_types, kap, lam)


def computeCR(cr, molecules, atomTypes, matrix):
	for cod in molecules:
		mol = molecules[cod]

		for prm in mol.parameters:
			prm.computeCR(cr, atomTypes, matrix)


def createEffectivePrms(atomTypes, molecules, charge_group_type):
	for cod in molecules:
		mol = molecules[cod]

		charge_group_type.setCG(mol)
		mol.createChargeTransferTopology()
		mol.createEffectiveAtomicCharges()
		mol.createLJPairs(atomTypes)


def printPrmValue(idx, cod, molecules):
	prm = molecules[cod].parameters[idx]
	print(prm.nam, prm.cur)

