import numpy as np
import sys
import os


import myExceptions


def runAna(conf, molecules, anaDir):
    it = conf.it
    allOutFile = anaDir + '/all_' + str(it) + '.out'
    allSumFile = anaDir + '/all_sum_' + str(it) + '.out'
    resFile = anaDir + '/res_' + str(it) + '.dat'
    newMolFile = anaDir + '/mol.dat'
    rmsdFile = anaDir + '/rmsd_' + str(it) + '.dat'

    with open(allOutFile, 'w') as allOut:
        for cod in molecules:
            mol = molecules[cod]

            if mol.run == True:
                # add running average of each job
                addSimprop(conf, mol)

                writeAllFile(mol, allOut)


def addSimprop(conf, mol):
    nJobs = conf.nJobs
    it = conf.it
    cod = mol.cod

    properties = {}
    startJob = 1
    for i in range(startJob, nJobs):
        outFileName  = 'dat_{}/o_{}_{}'.format(it, cod, i)
        if not os.path.exists(outFileName):
            raise myExceptions.NoSuchFile(outFileName)

        nRow = os.popen('wc ' + outFileName).read().split()[0]
        # if nRow is not correct ignore molecule (set run to 0)
        if not nRow in ['100', '120', '240', '1000']:
            mol.run = 0.0
            return

        with open(outFileName, 'r') as f:
            lines = f.readlines()

        lines = [row.strip().split() for row in lines]

        for row in lines:
            letter = row[0]
            tim = row[1]
            val = row[2]
            tim = float(tim)
            val = float(val)
            if not letter in properties:
                properties[letter] = {'tim': [tim]}
                for j in range(startJob, nJobs):
                    properties[letter][j] = []
                properties[letter][i].append(val)

            else:
                properties[letter]['tim'].append(tim)
                properties[letter][i].append(val)


    # do not add instantaneous values
    for letter in properties:
        tim = np.unique(properties[letter]['tim'])

        traj = np.zeros((tim.shape[0], nJobs), dtype=np.float)
        traj[:,0] = tim
        for j in range(startJob, nJobs):
            traj[:,j] = np.array(properties[letter][j])

        # mol.addTrajectory(letter, traj)
        avgs = traj[-1,1:]
        mol.addRunningAverages(letter, avgs)


def getMaxDev(prop):
    maxDev = 0.0
    for i in range(len(prop)):
        for j in range(i, len(prop)):
            iVal = prop[i]
            jVal = prop[j]
            dev = abs(iVal - jVal)
            if dev > maxDev:
                maxDev = dev
    return maxDev


def writeAllFile(mol, out):
    cod = mol.cod
    print(cod)

    frm = mol.frm
    run = '1.0'
    pre_sim = mol.pre_sim
    tem_sim = mol.tem_sim

    s = '{:6} {:8} {:3} {:6} {:6}'.format(cod, frm, run, pre_sim, tem_sim)

    properties = mol.properties
    for prop in properties:
        wei = prop.wei
        ref = prop.ref
        runningAverages = prop.runningAverages

        s = '{} {:3} {:6}'.format(s, wei, ref)
        for val in runningAverages:
            s = '{} {:6.1f}'.format(s, val)

    out.write(s)
