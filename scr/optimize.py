"""Optimization module based on scipy module.

SciPy ``optimize`` provides functions for minimizing (or maximizing)
objective functions, possibly subject to constraints.

Methods:
    runOptimization(conf, molecules, atomTypes)
    testParameter(prmsToOptmize)
    updateOirginalParameterValues(molecules)
    addSimProp(conf, molecules)
    getParametersToBeOptimized(atomTypes)
    setMinMax(conf, prmsToOptmize)
    minimize(conf, init, prmsToOptmize, atomTypes, molecules)
    targetFunction(init, prmsToOptmize, atomTypes, cr, eem, matrix, molecules, info)
    updatePrm(init, prmsToOptmize)
    updateEffectivePrm(cr, eem, matrix, molecules, atomTypes)
"""

from scipy import optimize
import shutil
import sys
import os

import molecules_utils
import ana


def runOptimization(conf, molecules, atomTypes):
    """Runs optimization of force-field parameters.

    :param conf: (configuration.Conf) Configuration object.
    :param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
    :param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.

    - Computes charge distribution for each molecule.
    - Updates values of effective parameters for each molecule.
    - Extracts ensemble averages ans sensitivities from simulation results.
    - Minimizes parameters.
    """

    charge_method = conf.charge_distribution_method

    molecules_utils.computeChargeDistribution(charge_method, molecules, atomTypes, conf.kappa, conf.lam)

    updateEffectivePrm(conf.cr, conf.scl_sig_NEI, conf.scl_eps_NEI, charge_method, conf.kappa, conf.lam, conf.matrix,
                       molecules, atomTypes)

    updateOirginalParameterValues(molecules)

    # get simulated results
    addSimProp(conf, molecules)

    prmsToOptmize = getParametersToBeOptimized(atomTypes)

    # set min/max values and return initial values
    init = setMinMax(conf, prmsToOptmize)

    optDir = conf.optDir
    # remove optDir/ if it already exists, and create a new one.
    if os.path.exists(optDir):
        shutil.rmtree(optDir)
    os.makedirs(optDir)

    optOutFile = conf.optOutFile
    with open(optOutFile, 'w') as optOut:
        minimize(conf, init, prmsToOptmize, atomTypes, molecules, optOut)


def testParameter(prmsToOptmize):
    """Test function. Checks if parameters are changed when one of its copies if changed."""

    prm2 = prmsToOptmize[0][0]
    print('prm2: ', prm2.iac, prm2.cur)

    prm2.cur = 1.0
    print('prm2: ', prm2.iac, prm2.cur)


def updateOirginalParameterValues(molecules):
    """Updates original values of parameters to their respective current values.

    :param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
    """

    for cod in molecules:
        mol = molecules[cod]
        parameters = mol.parameters
        for prm in parameters:
            prm.ori = prm.cur


def addSimProp(conf, molecules):
    """Extract ensemble averages and sensitivities from simulation results and save into molecule.

    :param conf: (configuration.Conf) Configuration object.
    :param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
    """

    for cod in molecules:
        mol = molecules[cod]

        # add running average of each job
        ana.addSimprop(conf, mol)

        # add sensitivity matrix
        ana.addSens(conf, mol)


def getParametersToBeOptimized(atomTypes):
    """Returns parameters to be optimized.

    :param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
    :return:
        prmsToOptmize: (list) List of parameters to be optimized.
    """

    prmsToOptmize = []

    for iac in atomTypes:
        atmTyp = atomTypes[iac]

        sig = atmTyp.sig
        eps = atmTyp.eps
        sig_2 = atmTyp.sig_2
        eps_2 = atmTyp.eps_2
        hrd = atmTyp.hrd
        eln = atmTyp.eln

        allParameters = [sig, eps, sig_2, eps_2, hrd, eln]
        parameters = []
        for prm in allParameters:
            if prm.rng != 0.0:
                parameters.append(prm)

        if len(parameters) == 0:
            continue

        if len(prmsToOptmize) == 0:
            prmsToOptmize.append([parameters[0]])
            parameters = parameters[1:]

        # check for symmetry
        indexes = []
        for i in range(len(parameters)):
            prm1 = parameters[i]

            j = 0
            while j < len(prmsToOptmize):
                prmList = prmsToOptmize[j]

                k = 0
                while k < len(prmList):
                    prm2 = prmsToOptmize[j][k]

                    if (prm1.hasSymmetricIAC()) and (prm1.typ == prm2.typ) and (prm1.symmetry == prm2.symmetry):
                        prmsToOptmize[j].append(prm1)
                        indexes.append(i)
                        j = len(prmsToOptmize)
                        k = len(prmsToOptmize)

                    k += 1
                j += 1

        for i in range(len(parameters)):
            if i not in indexes:
                prm = parameters[i]
                prmsToOptmize.append([prm])

    return prmsToOptmize


def setMinMax(conf, prmsToOptmize):
    """Sets minimal amd maximal range of parameter variation.

    :param conf: (configuration.Conf) Configuration object.
    :param prmsToOptmize: (list) List of parameters to be optimized.
    :return:
        init: (list) List with initial guess of parameter values.
    """

    rng_scl = conf.rng_scl
    init = []

    for prmList in prmsToOptmize:
        minVal = []
        maxVal = []
        iniValues = []
        iac = []

        for prm in prmList:
            val = prm.ori * (1.0 - prm.rng * rng_scl / 100.0)
            minVal.append(val)
            prm.min = val

            val = prm.ori * (1.0 + prm.rng * rng_scl / 100.0)
            maxVal.append(val)
            prm.max = val

            iniValues.append(prm.ori)

            iac.append(prm.iac)

        minVal = list(set(minVal))
        maxVal = list(set(maxVal))

        if len(minVal) != 1 and len(maxVal) != 1:
            print('Symmetrical parameters with different values: {}'.format(iac))
            sys.exit(1)

        iniValues = list(set(iniValues))
        init.append(iniValues[0])

    if len(init) != len(prmsToOptmize):
        print('\n\tWrong number of parameters to be optimized.')
        sys.exit(1)

    return init


def minimize(conf, init, prmsToOptmize, atomTypes, molecules, optOut):
    """Minimizes target function.

    :param conf: (configuration.Conf) Configuration object.
    :param init: (list) List with initial guess of parameter values.
    :param prmsToOptmize: (list) List of parameters to be optimized.
    :param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
    :param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
    :param optOut: (output object) Output for optimization results.
    """
    kappa = conf.kappa
    lam = conf.lam
    eem = conf.charge_distribution_method
    matrix = conf.matrix
    cr = conf.cr
    opt_nit = conf.opt_nit
    scl_sig_NEI = conf.scl_sig_NEI
    scl_eps_NEI = conf.scl_eps_NEI

    minimizer_kwargs = (prmsToOptmize, atomTypes, cr, scl_sig_NEI, scl_eps_NEI, eem, kappa, lam, matrix, molecules,
                        {'Nfeval': 0}, optOut)
    option = {'disp': True, 'maxiter': opt_nit}

    # x_fnc = targetFunction(init, prmsToOptmize, atomTypes, cr, , scl_sig_NEI, scl_eps_NEI eem, matrix, molecules,
    # {'Nfeval':0})
    ret = optimize.minimize(targetFunction, init, args=minimizer_kwargs, options=option, method='Nelder-Mead')
    print(ret)

    outPrmFile = conf.outPrmFile

    writePrm(ret.x, prmsToOptmize, atomTypes, matrix, cr, scl_sig_NEI, scl_eps_NEI, eem, kappa, lam, molecules,
             outPrmFile)


def targetFunction(init, prmsToOptmize, atomTypes, cr, scl_sig_NEI, scl_eps_NEI, eem, kappa, lam, matrix, molecules,
                   info, optOut):
    """Computes target function given value of parameters.

    :param init: (list) List with initial guess of parameter values.
    :param prmsToOptmize: (list) List of parameters to be optimized.
    :param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
    :param cr: (combiningRule) Combining rule.
    :param scl_sig_NEI: (float) Scaling factor for 1-4 sigma.
    :param scl_eps_NEI: (float) Scaling factor for 1-4 epsilon.
    :param eem: (ChargeDistributionMethod) Charge distribution method.
    :param kappa: (float) Kappa parameter used in charge distribution method.
    :param lam (float) Lambda parameter used in charge distribution method.
    :param matrix: (Matrix) Matrix with usage of C12(II) parameter.
    :param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
    :param info: Dictionary used to define frequency that the target function value is printed.
    :param optOut: (output object) Output for optimization results.
    :return:
        x_fnc: (float) Value of target function.
    """

    updatePrm(init, prmsToOptmize)
    updateEffectivePrm(cr, scl_sig_NEI, scl_eps_NEI, eem, kappa, lam, matrix, molecules, atomTypes)

    x_fnc = 0.0
    x_sum = 0.0

    for cod in molecules:
        mol = molecules[cod]

        run = mol.run
        if run:

            sens = mol.sens
            props = mol.properties
            parameters = mol.parameters

            for prop in props:
                ref = prop.ref
                wei = prop.wei
                sim = prop.sim
                propCode = prop.code
                propScale = prop.scale

                linearApproximationX = sim

                for prm in parameters:
                    ori = prm.ori
                    cur = prm.cur
                    prmName = prm.nam

                    derivative = sens.getDerivative(prmName, propCode)
                    linearApproximationX += derivative * (cur - ori)

                term = abs(linearApproximationX - ref)
                term *= wei / propScale

                x_fnc += term
                x_sum += wei

    x_fnc = x_fnc / x_sum

    if info['Nfeval'] % 10 == 0:
        optOut.write('{0:4d} {1: 3.6f}\n'.format(info['Nfeval'], x_fnc))
        # print(info['Nfeval'], x_fnc)

    info['Nfeval'] += 1
    return x_fnc


def updatePrm(init, prmsToOptmize):
    """Updates parameter values, and adjusts current value to [min, max].

    :param init: (list) List with initial guess of parameter values.
    :param prmsToOptmize: (list) List of parameters to be optimized.
    """

    for i in range(len(prmsToOptmize)):
        prmList = prmsToOptmize[i]

        for prm in prmList:
            newVal = init[i]
            prm.cur = newVal

            if newVal < prm.min:
                prm.cur = prm.min
                init[i] = prm.min
            elif newVal > prm.max:
                prm.cur = newVal
                init[i] = prm.max


def updateEffectivePrm(cr, scl_sig_NEI, scl_eps_NEI, eem, kappa, lam, matrix, molecules, atomTypes):
    """Updated values of effective parameters for each molecule.

    :param cr: (combiningRule) Combining rule.
    :param scl_sig_NEI: (float) Scaling factor for 1-4 sigma.
    :param scl_eps_NEI: (float) Scaling factor for 1-4 epsilon.
    :param eem: (ChargeDistributionMethod) Charge distribution method.
    :param kappa: (float) Kappa parameter used in charge distribution method.
    :param lam (float) Lambda parameter used in charge distribution method.
    :param matrix: (Matrix) Matrix with usage of C12(II) parameter.
    :param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
    :param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
    :return:
    """

    molecules_utils.computeChargeDistribution(eem, molecules, atomTypes, kappa, lam)
    molecules_utils.computeCR(cr, scl_sig_NEI, scl_eps_NEI, molecules, atomTypes, matrix)


def writePrm(val, prmsOptmized, atomTypes, matrix, cr, scl_sig_NEI, scl_eps_NEI, eem, kappa, lam, molecules, fileName):
    """Writes optimized parameters to txt file.

    :param val: (list) List with final optimized parameter values.
    :param prmsOptmized: (list) List of parameters to be optimized.
    :param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
    :param matrix: (Matrix) Matrix with usage of C12(II) parameter.
    :param cr: (combiningRule) Combining rule.
    :param scl_sig_NEI: (float) Scaling factor for 1-4 sigma.
    :param scl_eps_NEI: (float) Scaling factor for 1-4 epsilon.
    :param eem: (ChargeDistributionMethod) Charge distribution method.
    :param kappa: (float) Kappa parameter used in charge distribution method.
    :param lam (float) Lambda parameter used in charge distribution method.
    :param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
    :param fileName: (str) Name out output file.
    :return:
    """

    updatePrm(val, prmsOptmized)
    updateEffectivePrm(cr, scl_sig_NEI, scl_eps_NEI, eem, kappa, lam, matrix, molecules, atomTypes)

    with open(fileName, 'w') as out:
        out.write('{0:<5} {1:11} {2:13} {3:5} {4:13} {5:5} {6:8} {7:6} {8:8} {9:5} {10:13} {11:5} {12:13} {13}\n'
            .format('#iac', 'atnm', 'sig', 'rng', 'eps', 'rng', 'hrd', 'rng', 'eln', 'rng', 'sig_2', 'rng', 'eps_2',
                     'rng'))

        for iac in atomTypes:
            atomTyp = atomTypes[iac]
            typ = atomTyp.typ
            sig = atomTyp.sig.cur
            eps = atomTyp.eps.cur
            hrd = atomTyp.hrd.cur
            eln = atomTyp.eln.cur

            sig_rng = atomTyp.sig.rng
            eps_rng = atomTyp.eps.rng
            hrd_rng = atomTyp.hrd.rng
            eln_rng = atomTyp.eln.rng

            sig_2 = atomTyp.sig_2.cur
            eps_2 = atomTyp.eps_2.cur
            sig_rng_2 = atomTyp.sig_2.rng
            eps_rng_2 = atomTyp.eps_2.rng

            out.write('{0:<5} {1:10} {2:13.6E} {3:5.2f} {4:13.6E} {5:5.2f} {6:8.5f} {7:5.2f} {8:9.5f} {9:5.2f} '
                      '{10:13.6E} {11:5.2f} {12:13.6E} {13:5.2f}\n'
                      .format(iac, typ, sig, sig_rng, eps, eps_rng, hrd, hrd_rng, eln, eln_rng, sig_2, sig_rng_2, eps_2,
                              eps_rng_2))
