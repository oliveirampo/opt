"""This module provides methods for handling collection of molecules.

Methods:
	checkAtoms(molecules)
	computeChargeDistribution(eem, molecules, atom_types, kap, lam)
	computeCR(cr, molecules, atomTypes, matrix)
	createEffectivePrms(atomTypes, molecules, charge_group_type)
	printPrmValue(idx, cod, molecules)
"""


def checkAtoms(molecules):
	"""Checks atoms in molecules.

	:param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
	"""

	for cod in molecules:
		molecules[cod].checkAtoms()


def computeChargeDistribution(eem, molecules, atom_types, kap, lam):
	"""Computes charge distribution for each molecule.

	:param eem: (ChargeDistributionMethod) Charge distribution method.
	:param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
	:param atom_types: (collections.OrderedDict) Ordered dictionary of atom types.
	:param kap: (float) Kappa parameter.
	:param lam: (float) Lambda parameter
	"""

	max_order = 2

	for cod in molecules:
		mol = molecules[cod]

		chargeGroups = mol.CGs
		for cg in chargeGroups:
			eem.compute(max_order, mol, cg, atom_types, kap, lam)


def computeCR(cr, scl_sig_NEI, scl_eps_NEI, molecules, atomTypes, matrix):
	"""Computes combining rule for each molecule.

	:param cr: (combiningRule) Combining rule.
	:param scl_sig_NEI: (float) Scaling factor for 1-4 sigma.
	:param scl_eps_NEI: (float) Scaling factor for 1-4 epsilon.
	:param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
	:param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
	:param matrix: (Matrix) Matrix with usage of C12(II) parameters.
	:return:
	"""

	for cod in molecules:
		mol = molecules[cod]

		for prm in mol.parameters:
			prm.computeCR(cr, scl_sig_NEI, scl_eps_NEI, atomTypes, matrix)


def createEffectivePrms(atomTypes, molecules, charge_group_type):
	"""Creates parameters and adds to list of parameters for each molecule.

	:param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
	:param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
	:param charge_group_type: (Charge_group_type) Type of charge group.
	"""

	for cod in molecules:
		mol = molecules[cod]

		charge_group_type.setCG(mol)
		mol.createChargeTransferTopology()
		mol.createEffectiveAtomicCharges()
		mol.createLJPairs(atomTypes)


def printPrmValue(idx, cod, molecules):
	"""Prints name and current value of parameter for a given molecule.

	:param idx: (int) Parameter position in list of parameters.
	:param cod: (str) Code of molecule.
	:param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
	"""

	prm = molecules[cod].parameters[idx]
	print(prm.nam, prm.cur)
