from collections import OrderedDict
import pandas as pd
import numpy as np
import math
import sys
import os


from scr import myExceptions, molecules_utils
from scr.molecule import Molecule
import effectiveParameter
from scr.iac import IAC
from scr.matrix import Matrix
from scr.atom import Atom
from property import Dns, Hvp


def readFile(fileName):
	# print('Reading ' + fileName)
	with open(fileName, 'r') as f:
		lines = f.readlines()
	lines = [row.strip().split() for row in lines]
	return lines


def readConf(conf, fileName):
	# print('reading ' + fileName)
	confDict = {}

	lines = readFile(fileName)

	i = 0
	while i < len(lines):
		row = lines[i]

		if row[0].startswith('#'):
			i = i + 1

		elif row[0].startswith('BEGIN'):
			i, dictName, add_data = readConf_helper(i, lines)
			#print(dictName)
			#print(add_data)
			# TODO

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

		key = "EEM"
		conf.eem = confDict[key]

		key = "scl_sig_NEI"
		conf.scl_sig_NEI = confDict[key]

		key = "scl_eps_NEI"
		conf.scl_eps_NEI = confDict[key]

		key = "ignoreIAC"
		conf.ignoreIAC = confDict[key]

	except KeyError:
		raise myExceptions.MissingKeyWordError(key, "cnfFile")


def readConf_helper(i, lines):
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
	files = cnf.inpFiles
	for fileName in files:
		if not os.path.exists(fileName):
			raise myExceptions.NoSuchFile(fileName)


def readMolData(conf):
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
				eps_ref = line[11]

				dns = Dns(scl_dns, dns_wei, dns_ref)
				hvp = Hvp(scl_hvp, hvp_wei, hvp_ref)
				properties = [dns, hvp]

				mol = Molecule(cod, frm, run, pre_sim, tem_sim, properties, mlp_ref, blp_ref, eps_ref)
				molecules[cod] = mol

	return molecules


def readPrm(fileName):
	atomTypes = OrderedDict()
	lines = readFile(fileName)

	for line in lines:
		if line[0][0] != '#':
			iac = int(line[0])

			if iac in atomTypes:
				sys.exit(iac + ' is already in prm_IT.dat')

			typ       = line[1]
			sig       = float(line[2])
			rng_sig   = float(line[3])
			eps       = float(line[4])
			rng_eps   = float(line[5])
			hrd       = float(line[6])
			rng_hrd   = float(line[7])
			eln       = float(line[8])
			rng_eln   = float(line[9])
			sig_2     = float(line[10])
			rng_sig_2 = float(line[11])
			eps_2     = float(line[12])
			rng_eps_2 = float(line[13])

			atomType = IAC(iac, typ, sig, rng_sig, eps, rng_eps, hrd, rng_hrd, eln, rng_eln, sig_2, rng_sig_2, eps_2, rng_eps_2)
			atomTypes[iac] = atomType

	return atomTypes


def readPrmNei(atomTypes, fileName):
	# print('Reading ' + fileName)
	df = pd.read_csv(fileName, sep='\s+', comment='#', names=['iac', 'name', 'sig_nei', 'eps_nei'])

	for idx, row in df.iterrows():
		iac = str(row['iac'])
		iac = int(iac)

		if not iac in atomTypes:
			continue

		sig_nei = row['sig_nei']
		eps_nei = row['eps_nei']

		atomType = atomTypes[iac]
		atomType.sig_nei = sig_nei
		atomType.eps_nei = eps_nei
		atomType.fixed_nei = True


def readMatrix(fileName):
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
	nBnd = 0
	bnd = {}

	for i in range(idx, len(lines)):
		if lines[i][0][0] != '#':
			if lines[i][0] == 'END':
				return i, bnd

			if len(lines[i]) == 2:
				nBnd = lines[i][1]

			else:
				bnd[lines[i][0]] = lines[i][3]

	if len(bnd) != nBnd:
		sys.exit('NBTY({}) != {}'.format(nBnd, len(bnd)))


def readListBond(bnd, molecules, fileName):
	lines = readFile(fileName)

	for line in lines:
		if line[0][0] != '#':
			cod = line[0]
			if cod in molecules:
				idx1 = line[1]
				idx2 = line[2]
				dist = bnd[line[3]]

				atm1 = molecules[cod].getAtom(idx1)
				atm2 = molecules[cod].getAtom(idx2)

				if (not atm1.ignore) and (not atm2.ignore):
					atm1.addBndNrm(idx2, dist)
					atm2.addBndNrm(idx1, dist)


def readListAngle(ang, molecules, fileName):
	lines = readFile(fileName)

	for line in lines:
		if line[0][0] != '#':
			cod = line[0]
			if cod in molecules:
				idx1 = line[1]
				idx2 = line[2]
				idx3 = line[3]
				theta = ang[line[4]]

				atm1 = molecules[cod].getAtom(idx1)
				atm2 = molecules[cod].getAtom(idx2)
				atm3 = molecules[cod].getAtom(idx3)

				if (atm1.ignore == False) and (atm2.ignore == False) and (atm2.ignore == False):
					if (not atm1.isNrmNB(idx2)) and (not atm3.isNrmNB(idx2)):
						sys.exit('{} is not the central atom'.format(atm2.nam))

					d1 = atm1.getNrmBndDist(idx2)
					d3 = atm3.getNrmBndDist(idx2)
					theta = float(theta)
					dist = math.sqrt(d1 * d1 + d3 * d3 - 2 * d1 * d3 * math.cos(math.radians(theta)))

					atm1.addBndNei(idx3, dist)
					atm3.addBndNei(idx1, dist)


def readSymmetry(atomTypes, symTyp, fileName):
	if not os.path.exists(fileName):
		return

	lines = readFile(fileName)

	for line in lines:
		if line[0][0] != '#':
			iac = int(line[0])
			typ = line[1]
			sym = line[2]

			if not sym.isupper():
				sys.exit('\n\tWrong value for symmetry: {} \
				\n\tIt should be uppercase\n'.format(sym))

			if not iac in atomTypes:
				sys.exit('\n\tIAC = {} not found in prms.dat\n'.format(iac))

			if typ != atomTypes[iac].typ:
				sys.exit('\n\tWrong typ for IAC={}: {} != {}\n'.format(iac, typ, atomTypes[iac].typ))

			atomTypes[iac].addSymmetry(symTyp, sym)


def readSamTemplateFile(fileName):
	lines = readFile(fileName)
	return lines


