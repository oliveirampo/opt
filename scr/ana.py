import pandas as pd
import numpy as np
import sys
import os


import myExceptions
from sensitivity import Sensitivity


def runAna(conf, molecules, anaDir):
    it = conf.it
    allOutFile = anaDir + '/all_' + str(it) + '.out'
    allSumFile = anaDir + '/all_sum_' + str(it) + '.out'
    resFile = anaDir + '/res_' + str(it) + '.dat'
    newMolFile = anaDir + '/mol.dat'
    rmsdFile = anaDir + '/rmsd_' + str(it) + '.dat'

    with open(allOutFile, 'w') as allOut, open(allSumFile, 'w') as allSum:
        for cod in molecules:
            mol = molecules[cod]

            # add running average of each job
            addSimprop(conf, mol)

            if mol.run == True:
                writeAllFile(mol, allOut)

                writeAllSum(mol, allSum)

                addSens(conf, mol)
                print(mol.sens)

                print('TODO: writeResFile')


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


def addSens(conf, mol):
    cod = mol.cod
    it = conf.it
    nJobs = conf.nJobs

    sens = Sensitivity(pd.DataFrame(columns=['typ', 'idx1', 'idx2']))

    sensitivities = {}
    startJob = 1
    for i in range(startJob, nJobs):
        outName = 'dat_{}/s_{}_{}'.format(it, cod, i)

        df = pd.read_csv(outName, sep='\s+',
        names=['letter', 'tim', 'typ', 'idx1', 'idx2', 'nam1', 'nam2', 'val', 'sens'])

        lastTim = df['tim'].unique()[-1]
        df = df.loc[df['tim'] == lastTim]

        for letter in df['letter'].unique():
            dat = df.loc[df['letter'] == letter]

            if not letter in sensitivities:
                sensitivities[letter] = {}
                sensitivities[letter][i] = dat

                if len(sensitivities) == 1:
                    sens.addPrmInfo(dat['typ'], dat['idx1'], dat['idx2'])

            else:
                sensitivities[letter][i] = dat

    # get avg and maxDev
    for letter in sensitivities:
        nRows = sensitivities[letter][startJob].shape[0]
        nReplicas = len(sensitivities[letter].keys())
        data = np.zeros(shape=(nRows, nReplicas))

        # save data into nd-array
        col = 0
        for job in sensitivities[letter]:
            dat = sensitivities[letter][job]

            if dat.shape[0] != nRows:
                sys.exit('ERROR: Wrong number of rows')

            dat = dat.set_index(['typ', 'idx1', 'idx2'])
            dat = dat['sens'].values

            data[:, col] = dat
            col += 1

        # get avg
        avg = np.mean(data, axis=1)

        # get max dev
        dev = np.zeros(shape=(nRows,1))
        for i in range(data.shape[1]):
            for j in range(i + 1, data.shape[1]):
                dat = abs(data[:,i] - data[:,j])
                dat = np.reshape(dat, (dat.shape[0], 1))

                dev = np.concatenate((dev, dat), axis=1)

        maxDev = np.amax(dev, axis=1)

        sens.addDerivative(letter, avg, maxDev)

    # sens.getAllDerivatives('D')
    mol.sens = sens


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


def writeAllSum(mol, allSum):
    cod = mol.cod
    frm = mol.frm
    run = mol.run
    pre_sim = mol.pre_sim
    tem_sim = mol.tem_sim

    allSum.write('{:6} {:12} {:3} {:7.3f} {:5.0f} {} '
    .format(cod, frm, run, pre_sim, tem_sim, '/'))

    properties = mol.properties
    for prop in properties:
        wei = prop.wei
        ref = prop.ref
        sim = prop.sim
        unit = prop.unit
        dd = prop.maxDev

        dev = sim - ref
        err = sim
        if ref != 0.0:
            err = 100 * dev / ref

        if wei == 0.0:
            wei = '*'
        else:
            ref = '{:.1f}'.format(ref)
            sim = '{:.1f}'.format(sim)
            dev = '{:.1f}'.format(dev)
            err = '{:.1f}'.format(err)
            dd  = '{:.1f}'.format(dd)

        allSum.write('{:3} {:>8} {:6} {:>6} {:>6} ( {:>4} {:6} ) {} '
        .format(wei, ref, sim, dev, err, dd, unit, '/'))

    allSum.write('\n')














