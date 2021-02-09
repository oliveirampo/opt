"""Performs analysis of the simulation results.

Methods:
--------
    runAna(conf, molecules, anaDir):
    addSimprop(conf, mol):
    addSens(conf, mol):
    getMaxDev(prop):
    writeAllFile(mol, out):
    writeAllSum(mol, allSum):
    writeResFile(mol, out):
    writeMolData(mol, out):
    getExpSimData(molecules):
    writeRmsd(df, out):
    getRmsd(prop, df):
"""

import pandas as pd
import numpy as np
import sys
import os


import myExceptions
from sensitivity import Sensitivity


def runAna(conf, molecules, anaDir):
    """Runs analysis of simulation results and writes to txt file.

    :param conf: (configuration.Conf) Configuration object.
    :param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
    :param anaDir: (float) Directory for analysis results.

    The files created are:
        - allOutFile: Summary of ensemble averages.
        - allSumFile: Summary of simulation results.
        - resFile: Summary of sensitivity results.
        - newMolFile: List of molecules that did not vaporize.
        - rmsdFile: Root-mean-square deviations.
    """

    it = conf.it
    allOutFile = anaDir + '/all_' + str(it) + '.out'
    allSumFile = anaDir + '/all_sum_' + str(it) + '.out'
    resFile = anaDir + '/res_' + str(it) + '.dat'
    newMolFile = anaDir + '/mol.dat'
    rmsdFile = anaDir + '/rmsd_' + str(it) + '.dat'

    with open(allOutFile, 'w') as allOut, open(allSumFile, 'w') as allSum,\
            open(resFile, 'w') as res, open(newMolFile, 'w') as newMolData,\
            open(rmsdFile, 'w') as rmsdOut:

        for cod in molecules:
            mol = molecules[cod]

            # add running average of each job
            addSimprop(conf, mol)

            if mol.run:
                writeAllFile(mol, allOut)

                writeAllSum(mol, allSum)

                addSens(conf, mol)
                # print(mol.sens)

                writeResFile(mol, res)

                writeMolData(mol, newMolData)

        df = getExpSimData(molecules)
        writeRmsd(df, rmsdOut)


def addSimprop(conf, mol):
    """Extracts ensemble averages from simulation results and saves into molecule object.

    :param conf: (configuration.Conf) Configuration object.
    :param mol: (Molecule) Molecule.
    """

    nJobs = conf.nJobs
    it = conf.it
    cod = mol.cod

    properties = {}
    startJob = 1
    for i in range(startJob, nJobs):
        outFileName = 'dat_{}/o_{}_{}'.format(it, cod, i)
        if not os.path.exists(outFileName):
            raise myExceptions.NoSuchFile(outFileName)

        nRow = os.popen('wc ' + outFileName).read().split()[0]
        # if nRow is not correct ignore molecule (set run to 0)
        if nRow not in ['100', '120', '240', '1000']:
            mol.run = 0.0
            return

        with open(outFileName, 'r') as f:
            lines = f.readlines()

        lines = [row.strip().split() for row in lines]

        for row in lines:
            propCode = row[0]
            tim = row[1]
            val = row[2]
            tim = float(tim)
            val = float(val)
            if propCode not in properties:
                properties[propCode] = {'tim': [tim]}
                for j in range(startJob, nJobs):
                    properties[propCode][j] = []
                properties[propCode][i].append(val)

            else:
                properties[propCode]['tim'].append(tim)
                properties[propCode][i].append(val)

    # do not add instantaneous values
    for propCode in properties:
        tim = np.unique(properties[propCode]['tim'])

        traj = np.zeros((tim.shape[0], nJobs), dtype=np.float)
        traj[:, 0] = tim
        for j in range(startJob, nJobs):
            traj[:, j] = np.array(properties[propCode][j])

        # mol.addTrajectory(letter, traj)
        avgs = traj[-1, 1:]
        mol.addRunningAverages(propCode, avgs)


def addSens(conf, mol):
    """Extracts sensitivities from simulation results and saves into molecule object.

    :param conf: (configuration.Conf) Configuration object.
    :param mol: (Molecule) Molecule.
    """

    cod = mol.cod
    it = conf.it
    nJobs = conf.nJobs

    sens = Sensitivity(pd.DataFrame(columns=['typ', 'idx1', 'idx2', 'nam1', 'nam2', 'val']))

    sensitivities = {}
    startJob = 1
    for i in range(startJob, nJobs):
        outName = 'dat_{}/s_{}_{}'.format(it, cod, i)

        df = pd.read_csv(outName, sep='\s+',
        names=['propCod', 'tim', 'typ', 'idx1', 'idx2', 'nam1', 'nam2', 'val', 'sens'])

        lastTim = df['tim'].unique()[-1]
        df = df.loc[df['tim'] == lastTim]

        for propCod in df['propCod'].unique():
            dat = df.loc[df['propCod'] == propCod]

            if propCod not in sensitivities:
                sensitivities[propCod] = {}
                sensitivities[propCod][i] = dat

                if len(sensitivities) == 1:
                    sens.addPrmInfo(dat['typ'], dat['idx1'], dat['idx2'], dat['nam1'], dat['nam2'], dat['val'])

            else:
                sensitivities[propCod][i] = dat

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
        minVal = np.min(data, axis=1)
        maxVal = np.max(data, axis=1)
        maxDev = np.abs(minVal - maxVal)

        sens.addDerivative(letter, avg, maxDev)

    # sens.getAllDerivatives('D')
    mol.sens = sens


def getMaxDev(prop):
    """Gets maximum deviation from array of values.

    :param prop: (list, arr) Array of values.
    :return: (float) Maximum deviation.
    """

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
    """Writes summary of ensemble averages.

    :param mol: (Molecule) Molecule.
    :param out: (output object) Output file.
    :return:
    """

    cod = mol.cod

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
    """Writes summary of simulation results.

    :param mol: (Molecule) Molecule.
    :param allSum: (output object) Output file.
    """

    cod = mol.cod
    frm = mol.frm
    run = mol.run
    pre_sim = mol.pre_sim
    tem_sim = mol.tem_sim

    allSum.write('{:6} {:12} {:3} {:7.3f} {:5.0f} {} '.format(cod, frm, run, pre_sim, tem_sim, '/'))

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
            dd = '{:.1f}'.format(dd)

        allSum.write('{:3} {:>8} {:6} {:>6} {:>6} ( {:>4} {:6} ) {} '
        .format(wei, ref, sim, dev, err, dd, unit, '/'))

    allSum.write('\n')


def writeResFile(mol, out):
    """Writes summary of sensitivity results.

    :param mol: (Molecule) Molecule.
    :param out: (output object) Output file.
    """

    cod = mol.cod
    sens = mol.sens
    sens.writeResFile(cod, out)


def writeMolData(mol, out):
    """Writes list of molecules that did not vaporize.

    :param mol: (Molecule) Molecule.
    :param out: (output object) Output file.
    """

    cod = mol.cod
    frm = mol.frm
    run = mol.run
    pre_sim = mol.pre_sim
    tem_sim = mol.tem_sim
    mlp_ref = mol.mlp_ref
    blp_ref = mol.blp_ref
    eps_ref = mol.eps_ref

    s = ''
    properties = mol.properties
    for prop in properties:
        letter = prop.letter
        wei = prop.wei
        ref = prop.ref
        sim = prop.sim
        maxDev = prop.maxDev

        if letter == 'D' and sim < 400.0:
            run = False
            wei = 0.0

        elif letter == 'D' and maxDev > 100.0:
            run = False
            wei = 0.0

        s = '{} {:4} {:6.1f}'.format(s, wei, ref)

    if run:
        run = 1

    out.write('{:8} {:10} {:2} {:6} {:7}'
    .format(cod, frm, run, pre_sim, tem_sim))
    out.write(s)
    out.write('{:7} {:7} {:6}\n'.format(mlp_ref, blp_ref, eps_ref))


def getExpSimData(molecules):
    """Extracts experimental and simulated data and save into DataFrame.

    :param molecules: (Molecule) Molecule.
    :return:
        df: (DataFrame) Table with experimental and simulated results.
    """
    # get number of molecules
    # get number of data points for each property

    colNames = ['cod']
    df = pd.DataFrame(columns=colNames)

    for cod in molecules:
        mol = molecules[cod]
        run = mol.run
        if run:
            idx = df.shape[0]
            df.loc[idx, 'cod'] = cod

            properties = mol.properties
            for prop in properties:
                wei = prop.wei
                letter = prop.letter
                ref = np.nan
                sim = np.nan
                if wei != 0.0:
                    ref = prop.ref
                    sim = prop.sim

                df.loc[idx, 'ref_{}'.format(letter)] = ref
                df.loc[idx, 'sim_{}'.format(letter)] = sim
                df.loc[idx, 'diff_{}'.format(letter)] = sim - ref
                df.loc[idx, 'diff2_{}'.format(letter)] = (sim - ref) * (sim - ref)
                df.loc[idx, 'err_{}'.format(letter)] = 100 * (sim - ref) / ref

    return df


def writeRmsd(df, out):
    """Writes root-mean-square deviations to txt file.

    :param df: (DataFrame) Table with main results.
    :param out: (output object) Output file.
    """

    propLetters = []
    for col in df.columns[1:]:
        letter = col.split('_')[-1]
        if letter not in propLetters:
            propLetters.append(letter)

    df['short_cod'] = df['cod'].map(lambda x: x[0] + x[2])

    # get statistics per group of compounds (based on short_cod)
    df_grouped = df.groupby('short_cod')
    codes = list(df_grouped.groups.keys())
    for i in range(len(codes)):
        cod = codes[i]
        data = df_grouped.get_group(cod)

        # print('{:3} {:4}'.format(cod, 1), end=' ')
        out.write('{:3} {:4}'.format(cod, 1))

        for letter in propLetters:
            dat = data.dropna(subset=['ref_{}'.format(letter)])
            nData = dat.shape[0]

            if nData == 0:
                # print('{:4} {:>6} {:>6} {:>6} {:>7}'.format(nData, '-', '-', '-', '-', '-'), end=' ')
                out.write('{:4} {:>6} {:>6} {:>6} {:>7}'.format(nData, '-', '-', '-', '-', '-'))

            else:
                prop_avg, rmsd, aved, mean_err = getRmsd(letter, dat)

                # print('{:4} {:6.1f} {:6.1f} {:6.1f} {:7.1f}'.format(nData, rmsd, aved, mean_err, prop_avg), end=' ')
                out.write('{:4} {:6.1f} {:6.1f} {:6.1f} {:7.1f}'.format(nData, rmsd, aved, mean_err, prop_avg))

        # print('')
        out.write('\n')

    # get total statistics
    nMol = df['cod'].map(lambda x: x[:5]).unique().shape[0]

    # print('{:3} {:4}'.format('T', nMol), end=' ')
    out.write('{:3} {:4}'.format('T', nMol))

    for letter in propLetters:
        dat = df.dropna(subset=['ref_{}'.format(letter)])
        nData = dat.shape[0]
        prop_avg, rmsd, aved, mean_err = getRmsd(letter, dat)

        # print('{:4} {:6.1f} {:6.1f} {:6.1f} {:7.1f}'.format(nData, rmsd, aved, mean_err, prop_avg), end=' ')
        out.write('{:4} {:6.1f} {:6.1f} {:6.1f} {:7.1f}'.format(nData, rmsd, aved, mean_err, prop_avg))
    # print('')
    out.write('\n')


def getRmsd(prop, df):
    """Computes root-mean-square deviations.

    :param prop:
    :param df: (DataFrame) Table with main results.
    :return:
    """

    # prop_avg = np.mean(df[['ref_{}'.format(prop)]]).values[0]
    # rmsd = np.sqrt(np.mean(df[['diff2_{}'.format(prop)]], axis=1)).values[0]
    # aved = np.mean(df[['diff_{}'.format(prop)]], axis=1).values[0]
    # mean_err = np.mean(df[['err_{}'.format(prop)]], axis=1).values[0]
    # return prop_avg, rmsd, aved, mean_err
    return 0, 0, 0, 0
