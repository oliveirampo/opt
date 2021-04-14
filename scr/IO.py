"""The IO module provides methods for general input handling.

Methods:
	readFile(fileName)
	readConf(conf, fileName)
	readConf_helper(i, lines)
	checkInpFiles(cnf)
	readMolData(conf)
	readPrm(fileName)
	readPrmNei(atomTypes, fileName)
	readMatrix(fileName)
	readListAtom(conf, molecules, fileName)
	readRefBondAndAngle(fileName)
	readRefBondAndAngle(fileName)
	getRef(idx, lines)
	readListBond(bnd, molecules, fileName)
	readListAngle(ang, molecules, fileName)
	readSymmetry(atomTypes, symTyp, fileName)
	readSamTemplateFile(fileName)
"""

from collections import OrderedDict
import pandas as pd
import numpy as np
import math
import sys
import os

import myExceptions
import configuration
import molecules_utils
import effectiveParameter
from iac import IAC
from atom import Atom
from matrix import Matrix
from molecule import Molecule
from property import Dns, Hvp


def readFile(fileName):
	"""Reads lines of input file.

	:param fileName: (str) Name of input file.
	:return: (list) List with rows extracted from input file.
	"""

	with open(fileName, 'r') as f:
		lines = f.readlines()
	lines = [row.strip().split() for row in lines]
	return lines


def readConf(conf, fileName):
	"""Reads configuration file.

	:param conf: (Conf) Configuration object.
	:param fileName: (str) File name.
	:exception
		myExceptions.MissingKeyWordError
	"""

	confDict = {}

	lines = readFile(fileName)

	i = 0
	while i < len(lines):
		row = lines[i]

		if row[0].startswith('#'):
			i = i + 1

		elif row[0].startswith('BEGIN'):
			i, dictName, add_data = readConf_helper(i, lines)
			confDict[dictName] = add_data

		else:
			if len(row) == 1:
				confDict[row[0]] = np.nan

			elif len(row) == 2:
				confDict[row[0]] = row[1]

			if len(row) > 2:
				confDict[row[0]] = row[1:]
			
			i = i + 1

	try:
		key = "nJobs"
		conf.nJobs = confDict[key]

		key = "scl_dns"
		conf.scl_dns = confDict[key]

		key = "scl_hvp"
		conf.scl_hvp = confDict[key]

		key = "opt_nit"
		conf.opt_nit = confDict[key]

		key = "rng_scl"
		conf.rng_scl = confDict[key]

		key = "eq_stp"
		conf.eq_stp = confDict[key]

		key = "eq_frq"
		conf.eq_frq = confDict[key]

		key = "prd_stp"
		conf.prd_stp = confDict[key]

		key = "prd_frq"
		conf.prd_frq = confDict[key]

		key = "wall_time"
		conf.wall_time = confDict[key]

		key = "cr"
		conf.cr = confDict[key]

		key = "charge_group_type"
		conf.charge_group_type = confDict[key]

		key = "charge_distribution"
		conf.charge_distribution_method = confDict[key]

		key = "kappa"
		conf.kappa = confDict[key]

		key = "lam"
		conf.lam = confDict[key]

		key = "scl_sig_NEI"
		conf.scl_sig_NEI = confDict[key]

		key = "scl_eps_NEI"
		conf.scl_eps_NEI = confDict[key]

		key = "ignoreIAC"
		conf.ignoreIAC = confDict[key]

	except KeyError:
		raise myExceptions.MissingKeyWordError(key, "cnfFile")

	# Add additional configuration.
	try:
		plotConf = configuration.PlotConf()
		plotConf.map_iac_name = confDict['iac_name']
		plotConf.map_cod_family = confDict['cod_family']
		plotConf.map_cod_color = confDict['cod_color']
		plotConf.map_cod_marker = confDict['cod_marker']
		plotConf.settings = confDict['plot_settings']

		conf.plotConf = plotConf

	except KeyError:
		conf.plotConf = plotConf
		pass


def readConf_helper(i, lines):
	"""Helper function to read configuration file.

	:param i: (int) Row number.
	:param lines: (list) List of rows.
	:return:
		nextRow: (int) Row number of next row.
		dictName: (str)
		add_data: (OrderedDict)
	"""

	dictName = ''
	add_data = OrderedDict()

	nexRow = i
	while i < len(lines):
		row = lines[i]

		if row[0].startswith('#'):
			i = i + 1

		elif row[0].startswith('END'):
			nextRow = i + 1
			i = len(lines)

		else:
			if row[0].startswith('BEGIN'):
				dictName = row[0].replace('BEGIN_', '')

			else:
				add_data[row[0]] = row[1]

			i = i + 1

	return nextRow, dictName, add_data

	
def checkInpFiles(cnf):
	"""Checks if input files exist.

	:param cnf: (Conf) Configuration object.
	:exception
		myExceptions.NoSuchFile
	"""

	files = cnf.inpFiles
	for fileName in files:
		if not os.path.exists(fileName):
			raise myExceptions.NoSuchFile(fileName)


def readMolData(conf):
	"""Reads file with information about molecules.

	:param conf: (Conf) Configuration object.
	:return:
		molecules: (collections.OrderedDict) Ordered dictionary of molecules.
	"""

	fileName = conf.molDataFile
	scl_dns = conf.scl_dns
	scl_hvp = conf.scl_hvp

	lines = readFile(fileName)
	molecules = OrderedDict()

	for line in lines:
		if line[0][0] != '#':
			cod = line[0]

			# check if cod is already in molData
			if cod in molecules:
				sys.exit(cod + ' is already in molData')

			run = float(line[2])
			if run != 0.0:
				frm = line[1]
				pre_sim = float(line[3])
				tem_sim = float(line[4])
				dns_wei = float(line[5])
				dns_ref = float(line[6])
				hvp_wei = float(line[7])
				hvp_ref = float(line[8])
				mlp_ref = line[9]
				blp_ref = line[10]
				tem_cri = line[11]
				eps_ref = line[12]

				dns = Dns(scl_dns, dns_wei, dns_ref)
				hvp = Hvp(scl_hvp, hvp_wei, hvp_ref)
				properties = [dns, hvp]

				mol = Molecule(cod, frm, run, pre_sim, tem_sim, properties, mlp_ref, blp_ref, tem_cri, eps_ref)
				molecules[cod] = mol

	return molecules


def readPrm(fileName):
	"""Reads parameter file.

	:param fileName: (str) Input file name.
	:return:
		atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
	"""

	atomTypes = OrderedDict()
	lines = readFile(fileName)

	for line in lines:
		if line[0][0] != '#':
			iac = int(line[0])

			if iac in atomTypes:
				sys.exit(iac + ' is already in prm_IT.dat')

			typ = line[1]
			nam = line[2]
			sig = float(line[3])
			rng_sig = float(line[4])
			eps = float(line[5])
			rng_eps = float(line[6])
			hrd = float(line[7])
			rng_hrd = float(line[8])
			eln = float(line[9])
			rng_eln = float(line[10])
			sig_2 = float(line[11])
			rng_sig_2 = float(line[12])
			eps_2 = float(line[13])
			rng_eps_2 = float(line[14])

			atomType = IAC(iac, typ, nam, sig, rng_sig, eps, rng_eps, hrd, rng_hrd, eln, rng_eln, sig_2, rng_sig_2,
						   eps_2, rng_eps_2)
			atomTypes[iac] = atomType

	return atomTypes


def readPrmNei(atomTypes, fileName):
	"""Reads file with 1-4 parameters and adds data to atom types.

	:param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
	:param fileName: (str) Input file name.
	"""

	df = pd.read_csv(fileName, sep='\s+', comment='#', names=['iac', 'name', 'sig_nei', 'eps_nei'])

	for idx, row in df.iterrows():
		iac = str(row['iac'])
		iac = int(iac)

		if iac not in atomTypes:
			continue

		sig_nei = row['sig_nei']
		eps_nei = row['eps_nei']

		atomType = atomTypes[iac]
		atomType.sig_nei = sig_nei
		atomType.eps_nei = eps_nei
		atomType.fixed_nei = True


def readMatrix(fileName):
	"""Reads matrix file.

	:param fileName: (str) Input file name.
	:return:
		matrix: (Matrix) Matrix with information about C12(II) usage.
	"""

	matrix = Matrix()
	lines = readFile(fileName)

	for line in lines:
		if line[0].startswith('#'):
			continue

		iac1 = line[0]
		iac2 = line[1]

		matrix.add(iac1, iac2)

	return matrix


def readListAtom(conf, molecules, fileName):
	"""Reads list of atoms and add to molecules.

	:param conf: (Conf) Configuration object.
	:param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
	:param fileName: (str) Input file name.
	:return:
	"""

	lines = readFile(fileName)

	for line in lines:
		if line[0][0] != '#':
			cod = line[0]
			idx = line[1]
			nam = line[2]
			iac = line[3]
			val = line[4]

			charge = effectiveParameter.createEffectiveParameterFactory('Charge', idx, '', '', iac, '', nam, '', val)
			atom = Atom(idx, nam, iac, charge)

			if cod in molecules:
				molecules[cod].addAtom(atom, conf)

	molecules_utils.checkAtoms(molecules)


def readRefBondAndAngle(fileName):
	"""Reads reference bond and bond angles.

	:param fileName: (str) Input file name.
	:return:
		bnd: (dict) Dictionary of reference bonds.
		ang: (dict) Dictionary of reference bond angles.
	"""

	lines = readFile(fileName)

	idx = 0
	bnd = {}
	for i in range(len(lines)):
		if lines[i][0][0] != '#':
			if lines[i][0] == 'BONDSTRETCHTYPECODE':
				idx, bnd = getRef(i + 1, lines)

	ang = {}
	for i in range(idx, len(lines)):
		if lines[i][0][0] != '#':
			if lines[i][0] == 'BONDANGLEBENDTYPECODE':
				idx, ang = getRef(i + 1, lines)

	return bnd, ang


def getRef(idx, lines):
	"""Extracts reference values of bond or bond angle from row.

	:param idx: (int) Row number.
	:param lines: (list) Array of rows.
	"""

	nBnd = 0
	bnd = {}

	for i in range(idx, len(lines)):
		if lines[i][0][0] != '#':
			if lines[i][0] == 'END':
				return i, bnd

			if len(lines[i]) == 2:
				nBnd = lines[i][1]

			else:
				bond_type = lines[i][0]
				bond_distance = lines[i][3]
				bnd[bond_type] = bond_distance

	if len(bnd) != nBnd:
		sys.exit('NBTY({}) != {}'.format(nBnd, len(bnd)))


def readListBond(bnd, molecules, fileName):
	"""Reads list of bonds of molecules,
	and creates connectivity and distance matrices.

	:param bnd: (dict) Reference bonds.
	:param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
	:param fileName: (str) Input file name.
	"""

	lines = readFile(fileName)
	lines = np.asarray(lines)

	for cod in molecules:
		mol = molecules[cod]

		if mol.run:
			mol_rows = np.where(lines == cod)
			mol_bonds = lines[mol_rows[0], 1:]

			mol_bonds = mol_bonds.astype(object)
			for i in range(mol_bonds.shape[0]):
				bond_type = mol_bonds[i, -1]
				bond_dist = bnd[bond_type]
				mol_bonds[i, -1] = bond_dist

			mol_bonds[:, 0:2] = mol_bonds[:, 0:2].astype(int) - 1
			mol_bonds[:, -1] = mol_bonds[:, -1].astype(float)

			n_atoms = mol.nAtoms()
			connectivity = np.zeros([n_atoms, n_atoms], dtype=bool)
			distance_matrix = np.zeros([n_atoms, n_atoms])

			for i in range(mol_bonds.shape[0]):
				atom_pos_1 = mol_bonds[i, 0]
				atom_pos_2 = mol_bonds[i, 1]
				distance = mol_bonds[i, 2]

				connectivity[atom_pos_1, atom_pos_2] = True
				connectivity[atom_pos_2, atom_pos_1] = True

				distance_matrix[atom_pos_1, atom_pos_2] = distance
				distance_matrix[atom_pos_2, atom_pos_1] = distance

			mol.connectivity = connectivity
			mol.distance_matrix = distance_matrix


def readListAngle(ang, molecules, fileName):
	"""Reads list of angles
	and updates distance matrix.

	:param ang: (dict) Reference bond angles.
	:param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
	:param fileName: (str) Input file name.
	"""

	lines = readFile(fileName)

	for line in lines:
		if line[0][0] != '#':
			cod = line[0]
			if cod in molecules:
				idx1 = int(line[1])
				idx2 = int(line[2])
				idx3 = int(line[3])
				theta = float(ang[line[4]])

				mol = molecules[cod]

				atm1 = mol.getAtom(idx1)
				atm2 = mol.getAtom(idx2)
				atm3 = mol.getAtom(idx3)

				if (atm1.ignore is False) and (atm2.ignore is False) and (atm3.ignore is False):

					bond_1_2 = mol.are_bonded(idx1, idx2)
					bond_2_3 = mol.are_bonded(idx2, idx3)

					if (bond_1_2 is False) or (bond_2_3 is False):
						sys.exit('Problem while reading list of angles: Atom {} is not a central atom'.format(atm2.nam))

					d1 = mol.get_bond_distance(idx1, idx2)
					d3 = mol.get_bond_distance(idx2, idx3)

					dist = math.sqrt(d1 * d1 + d3 * d3 - 2 * d1 * d3 * math.cos(math.radians(theta)))

					mol.add_bond_distance(idx1, idx3, dist)


def readSymmetry(atomTypes, symTyp, fileName):
	"""Reads symmetry file
	and mark parameters as symmetric.
	Symmetric parameters are optimized as one.

	:param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
	:param symTyp: (str) Code for type of symmetry (sig or eps).
	:param fileName: (str) Input file name.
	:return:
	"""

	if not os.path.exists(fileName):
		return

	lines = readFile(fileName)

	for line in lines:
		if line[0][0] != '#':
			iac = int(line[0])
			nam = line[1]
			sym = line[2]

			if not sym.isupper():
				sys.exit('\n\tWrong value for symmetry: {} \
				\n\tIt should be uppercase\n'.format(sym))

			if iac not in atomTypes:
				sys.exit('\n\tIAC = {} not found in prms.dat\n'.format(iac))

			if nam != atomTypes[iac].nam:
				sys.exit('\n\tWrong typ for IAC={}: {} != {}\n'.format(iac, nam, atomTypes[iac].typ))

			atomTypes[iac].addSymmetry(symTyp, sym)


def readSamTemplateFile(fileName):
	"""Reads sam template file.

	:param fileName: (str) Input file name.
	:return:
		lines: (list) List with rows.
	"""

	lines = readFile(fileName)
	return lines


def readVdWRadii(fileName):
	"""Reads list of van der Waals radii.

	:param fileName: (str) Input file name.
	:return
		vdw: (dict) Dictionary of van der Waals radii.
	"""

	lines = readFile(fileName)
	vdw = {}

	for line in lines:
		if line[0][0] != '#':
			atom = line[0]
			radius = line[1]

			vdw[atom] = radius
	return vdw


def readPrmCrFile(fileName):
	"""Reads file with parameter for linear combination of combining rules.

	:param fileName: (str) Name of file.
	:return:
	"""

	lines = readFile(fileName)
	crPrms = {}

	for line in lines:
		if line[0][0] != '#':
			prm = line[0]
			value = line[1]
			val_range = line[2]

			value = float(value)
			val_range = float(val_range)

			crPrms[prm] = {'val': value, 'rng': val_range}
	return crPrms
