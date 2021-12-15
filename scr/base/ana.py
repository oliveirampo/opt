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
import math
import sys
import os


# import scr.base.myExceptions
from scr.base.sensitivity import Sensitivity


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

    for cod in molecules:
        # print(cod)
        extractData(conf, molecules[cod])

    df = getExpSimData(molecules)

    with open(allOutFile, 'w') as allOut, open(allSumFile, 'w') as allSum,\
            open(resFile, 'w') as res, open(newMolFile, 'w') as newMolData:
        for cod in molecules:
            writeData(molecules[cod], allOut, allSum, res, newMolData)

    with open(rmsdFile, 'w') as rmsdOut:
        writeRmsd(df, conf, rmsdOut)


def extractData(conf, mol):
    # add running average of each job
    addSimprop(conf, mol)

    if mol.run:
        # print(mol.cod)
        addSens(conf, mol)


def writeData(mol, allOut, allSum, res, newMolData):
    if mol.run:
        writeAllFile(mol, allOut)
        writeResFile(mol, res)

    writeMolData(mol, newMolData)
    writeAllSum(mol, allSum)


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

    for i in range(startJob, nJobs + 1):
        addSimData(i, it, nJobs, startJob, properties, cod, mol)

    for propCode in properties:
        addSimRunningAverage(nJobs, startJob, propCode, properties, mol)


def addSimData(i, it, nJobs, startJob, properties, cod, mol):
    outFileName = 'dat_{}/o_{}_{}'.format(it, cod, i)
    if not os.path.exists(outFileName):
        # raise myExceptions.NoSuchFile(outFileName)
        mol.run = False
        return

    nRow = os.popen('wc ' + outFileName).read().split()[0]
    # if nRow is not correct ignore molecule (set run to 0)
    if nRow not in ['100', '120', '200', '240']:
        mol.run = False
        return

    with open(outFileName, 'r') as f:
        lines = f.readlines()

    lines = [row.strip().split() for row in lines]
    for row in lines:
        addSimInstantaneousValues(i, nJobs, startJob, row, properties)


def addSimInstantaneousValues(i, nJobs, startJob, row, properties):
    propCode = row[0]
    tim = row[1]
    val = row[2]
    tim = float(tim)
    val = float(val)
    if propCode not in properties:
        properties[propCode] = {'tim': [tim]}
        for j in range(startJob, nJobs + 1):
            properties[propCode][j] = []
        properties[propCode][i].append(val)

    else:
        properties[propCode]['tim'].append(tim)
        properties[propCode][i].append(val)


def addSimRunningAverage(nJobs, startJob, propCode, properties, mol):
    if mol.run:
        tim = np.unique(properties[propCode]['tim'])

        traj = np.zeros((tim.shape[0], nJobs + 1), dtype=np.float)
        traj[:, 0] = tim
        for j in range(startJob, nJobs + 1):
            traj[:, j] = np.array(properties[propCode][j])

        # mol.addTrajectory(letter, traj)
        avgs = traj[-1, 1:]

        if propCode == 'D':
            if avgs[-1] < 340:
                mol.run = False

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
        prepareSens(it, i, cod, sensitivities, sens)

    for letter in sensitivities:
        addSensData(startJob, letter, sensitivities, sens)

    mol.sens = sens


def prepareSens(it, i, cod, sensitivities, sens):
    outName = 'dat_{}/s_{}_{}'.format(it, cod, i)

    df = pd.read_csv(outName, sep='\s+',
                     names=['propCod', 'tim', 'typ', 'idx1', 'idx2', 'nam1', 'nam2', 'val', 'sens'])

    lastTim = df['tim'].unique()[-1]
    df = df.loc[df['tim'] == lastTim]

    for propCod in df['propCod'].unique():
        prepareSensHelper(i, propCod, df, sensitivities, sens)


def prepareSensHelper(i, propCod, df, sensitivities, sens):
    dat = df.loc[df['propCod'] == propCod]

    if propCod not in sensitivities:
        sensitivities[propCod] = {}
        sensitivities[propCod][i] = dat

        if len(sensitivities) == 1:
            # Adds information about parameters to DataFrame.
            sens.addPrmInfo(dat['typ'], dat['idx1'], dat['idx2'], dat['nam1'], dat['nam2'], dat['val'])

    else:
        sensitivities[propCod][i] = dat


def addSensData(startJob, letter, sensitivities, sens):
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


def getErr(prop):
    """Returns error on the mean at 95% confidence interval over repetitions.

    :param prop: (list, arr) Array of values.
    :return:
        err: (float) Error.
    """

    std = np.std(prop, ddof=1)
    N = len(prop)
    err = 1.96 * std / math.sqrt(N)
    return err


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
    run = mol.run
    pre_sim = mol.pre_sim
    tem_sim = mol.tem_sim

    s = '{:6} {:12} {:3} {:6} {:6}'.format(cod, frm, run, pre_sim, tem_sim)

    properties = mol.properties
    for prop in properties:
        wei = prop.wei
        ref = prop.ref
        runningAverages = prop.runningAverages

        s = '{} {:>8} {:7}'.format(s, wei, ref)
        for val in runningAverages:
            s = '{} {:6.1f}'.format(s, val)

    s = s + '\n'

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

    allSum.write('{:6} {:20} {:3} {:7.3f} {:5.0f} {} '.format(cod, frm, run, pre_sim, tem_sim, '/'))

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

        if np.isnan(sim):
            ref = '{:.1f}'.format(ref)
            sim = '{:.1f}'.format(0.0)
            dev = '-'
            err = '-'
            dd = '-'

        else:
            if wei == 0.0:
                wei = '*'
                ref = '-'
                sim = '{:.1f}'.format(sim)
                dev = '-'
                err = '-'
                dd = '-'
                # dd = '{:.1f}'.format(dd)

            else:
                ref = '{:.1f}'.format(ref)
                sim = '{:.1f}'.format(sim)
                dev = '{:.1f}'.format(dev)
                err = '{:.1f}'.format(err)
                dd = '{:.1f}'.format(dd)

        allSum.write('{:3} {:>8} {:6} {:>7} {:>6} ( {:>6} {:6} ) {} '
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
    tem_cri = mol.tem_cri
    eps_ref = mol.eps_ref

    s = ''
    properties = mol.properties
    for prop in properties:
        letter = prop.code
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

    out.write('{:8} {:20} {:2} {:6} {:7}'
    .format(cod, frm, run, pre_sim, tem_sim))
    out.write(s)
    out.write('{:7} {:7} {:6} {:6}\n'.format(mlp_ref, blp_ref, tem_cri, eps_ref))


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
                letter = prop.code
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


def writeRmsd(df, conf, out):
    """Writes root-mean-square deviations to txt file.

    :param df: (DataFrame) Table with main results.
    :param conf: (configuration.Conf) Configuration object.
    :param out: (output object) Output file.
    """

    propLetters = []
    for col in df.columns[1:]:
        letter = col.split('_')[-1]
        if letter not in propLetters:
            propLetters.append(letter)

    map_cod_family = conf.plotConf.map_cod_family
    fam_ordered, nhb_group, hbd_group = get_ordered_family_code(map_cod_family)

    df = df.apply(assign_family, axis=1, map_cod_family=map_cod_family, nhb_group=nhb_group, hbd_group=hbd_group)
    ref_columns = ['ref_' + prop for prop in propLetters]

    # get statistics per group of compounds (based on short_cod)
    df_grouped = df.groupby(['fam', 'n_fg'])

    # oder groups
    codes = list(df_grouped.groups.keys())
    codes = sorted(codes, key=lambda x: fam_ordered.index(x[0]))

    for i in range(len(codes)):
        cod = codes[i]
        data = df_grouped.get_group(cod)

        n_fg = cod[1]
        if n_fg == '0':
            n_fg = '-'

        # For each fam + n_fg
        write_rmsd_helper(data, n_fg, cod, propLetters, ref_columns, out)

        # For each fam
        n_fgs = df[df['fam'] == cod[0]]['n_fg'].unique()
        last_n_fg = n_fgs[-1]
        if n_fg != '0' and n_fg == last_n_fg:
            data = df[df['fam'] == cod[0]]
            write_rmsd_helper(data, n_fgs[0] + '-' + n_fgs[-1], cod, propLetters, ref_columns, out)

    # get total statistics
    write_rmsd_helper(df[df['group'] == 'HAL'], '-', ['HAL'], propLetters, ref_columns, out)
    write_rmsd_helper(df[df['group'] == 'NHB'], '-', ['NHB'], propLetters, ref_columns, out)
    write_rmsd_helper(df[df['group'] == 'HBD'], '-', ['HBD'], propLetters, ref_columns, out)
    write_rmsd_helper(df, '-', ['T'], propLetters, ref_columns, out)


def get_ordered_family_code(map_cod_family):
    fam_ordered = []
    for key, value in map_cod_family.items():
        fam = value
        if value == 'HAL':
            fam = key.split('_')[0]
        if fam not in fam_ordered:
            fam_ordered.append(fam)

    nhb_group = ['ROR', 'RCOH', 'RCOR', 'RCOOR']
    hbd_group = ['ROH', 'RCOOH', 'RN', 'RN$_2$', 'RCON']
    return fam_ordered, nhb_group, hbd_group


def assign_family(s,  map_cod_family, nhb_group, hbd_group):
    cod = s['cod'][0] + '_' + s['cod'][2]
    fam = map_cod_family[cod]
    s['fam'] = fam

    s['n_fg'] = s['cod'][2]
    if fam == 'MIX':
        s['n_fg'] = '0'

    s['group'] = s['fam']
    if fam in nhb_group:
        s['group'] = 'NHB'
    elif fam in hbd_group:
        s['group'] = 'HBD'

    if fam == 'HAL':
        s['fam'] = s['cod'][0]

    return s


def write_rmsd_helper(data, n_fg, cod, propLetters, ref_columns, out):
    data_non_empty = data.dropna(subset=ref_columns, how='all')
    nMol = getNumberOfMolecules(data_non_empty)

    # print('{:6} {:4}'.format(cod, nMol), end=' ')
    out.write('{:6} {:4} {:4}'.format(cod[0], n_fg, nMol))

    for letter in propLetters:
        dat = data.dropna(subset=['ref_{}'.format(letter)])
        nData = dat.shape[0]

        if nData == 0:
            # print('{:6} {:>6} {:>6} {:>6} {:>7}'.format(nData, '-', '-', '-', '-', '-'), end=' ')
            out.write('{:6} {:>6} {:>6} {:>6} {:>7}'.format(nData, '-', '-', '-', '-', '-'))

        else:
            prop_avg, rmsd, aved, mad = getRmsd(letter, dat)

            # print('{:6} {:6.1f} {:6.1f} {:6.1f} {:7.1f}'.format(nData, rmsd, aved, mad, prop_avg), end=' ')
            out.write('{:6} {:6.1f} {:6.1f} {:6.1f} {:7.1f}'.format(nData, rmsd, aved, mad, prop_avg))

        # print('')
    out.write('\n')


def getRmsd(prop, df):
    """Computes root-mean-square deviations.

    :param prop:
    :param df: (DataFrame) Table with main results.
    :return:
    """

    prop_avg = np.mean(df[['ref_{}'.format(prop)]]).values[0]
    rmsd = np.sqrt(np.mean(df[['diff2_{}'.format(prop)]], axis=0)).values[0]
    aved = np.mean(df[['diff_{}'.format(prop)]], axis=0).values[0]
    mad = np.mean(df[['diff_{}'.format(prop)]].abs(), axis=0).values[0]
    # mean_err = np.mean(df[['err_{}'.format(prop)]], axis=0).values[0]

    return prop_avg, rmsd, aved, mad


def getNumberOfMolecules(df):
    """Returns number of unique molecules.

    :param df: (pandas DataFrame) Table with data.
    :return: nMol (int) Number of unique molecules.
    """

    nMol = df['cod'].map(lambda x: x[:5]).unique().shape[0]
    return nMol
