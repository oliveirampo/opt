"""Optimization module based on scipy module.

SciPy ``optimize`` provides functions for minimizing (or maximizing)
objective functions, possibly subject to constraints.

Methods:
    runOptimization(conf, molecules, atomTypes)
    testParameter(prmsToOptimize)
    updateOriginalParameterValues(molecules)
    addSimProp(conf, molecules)
    getParametersToBeOptimized(atomTypes)
    setMinMax(conf, prmsToOptimize)
    minimize(conf, init, prmsToOptimize, atomTypes, molecules)
    targetFunction(init, prmsToOptimize, atomTypes, cr, eem, matrix, molecules, info)
    updatePrm(init, prmsToOptimize)
    updateEffectivePrm(cr, eem, matrix, molecules, atomTypes)
"""

# from datetime import datetime
from scipy import optimize
import shutil
import math
import sys
import os

from scr.base import molecules_utils
from scr.base import ana


def runOptimization(conf, crPrms, molecules, atomTypes):
    """Runs optimization of force-field parameters.

    :param conf: (configuration.Conf) Configuration object.
    :param crPrms: (dict) Parameters to do linear combination of combining rules.
    :param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
    :param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.

    - Computes charge distribution for each molecule.
    - Updates values of effective parameters for each molecule.
    - Extracts ensemble averages ans sensitivities from simulation results.
    - Minimizes parameters.
    """

    print("Computing charges")
    charge_method = conf.charge_distribution_method

    molecules_utils.computeChargeDistribution(charge_method, molecules, atomTypes, conf.kappa, conf.lam)

    print("Updating effective parameters")
    updateEffectivePrm(conf.cr, crPrms, conf.scl_sig_NEI, conf.scl_eps_NEI, charge_method, conf.kappa, conf.lam,
                       conf.matrix, molecules, atomTypes)

    updateOriginalParameterValues(molecules)

    print("Adding simulated results")
    # get simulated results
    addSimProp(conf, molecules)

    print("Selecting parameters to be optimized")
    prmsToOptimize = getParametersToBeOptimized(atomTypes)

    # set min/max values and return initial values
    init = setMinMax(conf, prmsToOptimize)

    optDir = conf.optDir
    # remove optDir/ if it already exists, and create a new one.
    if os.path.exists(optDir):
        shutil.rmtree(optDir)
    os.makedirs(optDir)

    print("Optimizing")
    optOutFile = conf.optOutFile
    with open(optOutFile, 'w') as optOut:
        minimize(conf, crPrms, init, prmsToOptimize, atomTypes, molecules, optOut)


def testParameter(prmsToOptimize):
    """Test function. Checks if parameters are changed when one of its copies if changed."""

    prm2 = prmsToOptimize[0][0]
    print('prm2: ', prm2.iac, prm2.cur)

    prm2.cur = 1.0
    print('prm2: ', prm2.iac, prm2.cur)


def updateOriginalParameterValues(molecules):
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

    molecules = [addSimPropHelper(conf, molecules[cod]) for _, cod in enumerate(molecules)]
    return molecules


def addSimPropHelper(conf, mol):
    # add running average of each job
    ana.addSimprop(conf, mol)

    # add sensitivity matrix
    ana.addSens(conf, mol)


def getParametersToBeOptimized(atomTypes):
    """Returns parameters to be optimized.

    :param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
    :return:
        prmsToOptimize: (list) List of parameters to be optimized.
    """

    prmsToOptimize = []

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

        if len(prmsToOptimize) == 0:
            prmsToOptimize.append([parameters[0]])
            parameters = parameters[1:]

        # check for symmetry
        indexes = []
        for i in range(len(parameters)):
            prm1 = parameters[i]

            j = 0
            while j < len(prmsToOptimize):
                prmList = prmsToOptimize[j]

                k = 0
                while k < len(prmList):
                    prm2 = prmsToOptimize[j][k]

                    if (prm1.hasSymmetricIAC()) and (prm1.typ == prm2.typ) and (prm1.symmetry == prm2.symmetry):
                        prmsToOptimize[j].append(prm1)
                        indexes.append(i)
                        j = len(prmsToOptimize)
                        k = len(prmsToOptimize)

                    k += 1
                j += 1

        for i in range(len(parameters)):
            if i not in indexes:
                prm = parameters[i]
                prmsToOptimize.append([prm])

    return prmsToOptimize


def setMinMax(conf, prmsToOptimize):
    """Sets minimal amd maximal range of parameter variation.

    :param conf: (configuration.Conf) Configuration object.
    :param prmsToOptimize: (list) List of parameters to be optimized.
    :return:
        init: (list) List with initial guess of parameter values.
    """

    rng_scl = conf.rng_scl
    init = []

    for prmList in prmsToOptimize:
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

    if len(init) != len(prmsToOptimize):
        print('\n\tWrong number of parameters to be optimized.')
        sys.exit(1)

    return init


def minimize(conf, crPrms, init, prmsToOptimize, atomTypes, molecules, optOut):
    """Minimizes target function.

    :param conf: (configuration.Conf) Configuration object.
    :param crPrms: (dict) Parameters to do linear combination of combining rules.
    :param init: (list) List with initial guess of parameter values.
    :param prmsToOptimize: (list) List of parameters to be optimized.
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

    minimizer_kwargs = (prmsToOptimize, atomTypes, cr, crPrms, scl_sig_NEI, scl_eps_NEI, eem, kappa, lam, matrix,
                        molecules, {'Nfeval': 0}, optOut)
    option = {'disp': True, 'maxiter': opt_nit}

    # x_fnc = targetFunction(init, prmsToOptimize, atomTypes, cr, , scl_sig_NEI, scl_eps_NEI eem, matrix, molecules,
    # {'Nfeval':0})
    ret = optimize.minimize(targetFunction, init, args=minimizer_kwargs, options=option, method='Nelder-Mead')
    # print(ret)

    outPrmFile = conf.outPrmFile

    updatePrm(ret.x, prmsToOptimize, atomTypes, cr, crPrms)
    updateEffectivePrm(cr, crPrms, scl_sig_NEI, scl_eps_NEI, eem, kappa, lam, matrix, molecules, atomTypes)

    writePrm(atomTypes, outPrmFile)

    fileName = '00_inp/prm_cr_{}.dat'.format(conf.it + 1)
    write_cr_prm(crPrms, fileName)


def targetFunction(init, prmsToOptimize, atomTypes, cr, crPrms, scl_sig_NEI, scl_eps_NEI, eem, kappa, lam, matrix,
                   molecules, info, optOut):
    """Computes target function given value of parameters.

    :param init: (list) List with initial guess of parameter values.
    :param prmsToOptimize: (list) List of parameters to be optimized.
    :param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
    :param cr: (combiningRule) Combining rule.
    :param crPrms: (dict) Parameters to do linear combination of combining rules.
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

    updatePrm(init, prmsToOptimize, atomTypes, cr, crPrms)
    updateEffectivePrm(cr, crPrms, scl_sig_NEI, scl_eps_NEI, eem, kappa, lam, matrix, molecules, atomTypes)

    x_fnc = [0.0]
    x_sum = [0.0]

    [getX(x_fnc, x_sum, molecules[mol]) for mol in molecules]

    x_fnc = x_fnc[0] / x_sum[0]

    if info['Nfeval'] % 10 == 0:
        optOut.write('{0:4d} {1: 3.6f}\n'.format(info['Nfeval'], x_fnc))
        print(info['Nfeval'], x_fnc)

    info['Nfeval'] += 1
    return x_fnc


def getX(x_fnc, x_sum, mol):
    if mol.run:
        sens = mol.sens
        props = mol.properties
        parameters = mol.parameters

        # for prop in props:
        [getLinearApproximation(prop, sens, parameters, x_fnc, x_sum) for prop in props]


def getLinearApproximation(prop, sens, parameters, x_fnc, x_sum):
    ref = prop.ref
    wei = prop.wei
    sim = prop.sim
    propCode = prop.code
    propScale = prop.scale

    linearApproximationX = [sim]

    # for prm in parameters:
    [addDerivative(prm, propCode, sens, linearApproximationX) for prm in parameters]

    term = abs(linearApproximationX[0] - ref)
    term *= wei / propScale

    x_fnc[0] += term
    x_sum[0] += wei


def addDerivative(prm, propCode, sens, linearApproximationX):
    ori = prm.ori
    cur = prm.cur
    prmName = prm.nam

    derivative = sens.getDerivative(prmName, propCode)
    linearApproximationX[0] += derivative * (cur - ori)


def updatePrm(init, prmsToOptimize, atomTypes, cr, crPrms):
    """Updates parameter values, and adjusts current value to [min, max].

    :param init: (list) List with initial guess of parameter values.
    :param prmsToOptimize: (list) List of parameters to be optimized.
    :param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
    :param cr: (combiningRule) Combining rule.
    :param crPrms: (dict) Parameters to do linear combination of combining rules.
    """

    for i in range(len(prmsToOptimize)):
        prmList = prmsToOptimize[i]

        for prm in prmList:
            newVal = init[i]
            prm.cur = newVal

            if newVal < prm.min:
                prm.cur = prm.min
                init[i] = prm.min
            elif newVal > prm.max:
                prm.cur = newVal
                init[i] = prm.max

    updateIAC(atomTypes, crPrms, cr)


def updateIAC(atomTypes, crPrms, cr):
    """Update sig-II/eps-II values.

    :param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
    :param cr: (combiningRule) Combining rule.
    :param crPrms: (dict) Parameters to do linear combination of combining rules.
    """

    alpha = crPrms['alpha']['val']

    # update sig_2 and eps_2 values
    for iac in atomTypes:
        atmTyp = atomTypes[iac]

        sig = atmTyp.sig
        eps = atmTyp.eps
        sig_2 = atmTyp.sig_2
        eps_2 = atmTyp.eps_2

        sigi_1 = sig.cur
        sigj_1 = sig.cur
        epsi_1 = eps.cur
        epsj_1 = eps.cur

        sigi_2 = sig_2.cur
        sigj_2 = sig_2.cur
        epsi_2 = eps_2.cur
        epsj_2 = eps_2.cur

        # do not update
        if sig_2.rng == 0.0 and eps_2.rng == 0.0:
            continue

        # update eps_2
        elif sig_2.rng == 0.0:
            sigij_1 = cr.getSigma(sigi_1, sigj_1, alpha)
            epsij_1 = cr.getEpsilon(epsi_1, epsj_1, sigi_1, sigj_1, alpha)
            epsij_2 = cr.getEpsilon(epsi_2, epsj_2, sigi_2, sigj_2, alpha)

            c6 = 4.0 * epsij_1 * math.exp(6.0 * math.log(sigij_1))

            sigij_2 = math.exp((1.0 / 6.0) * math.log(c6 / (4 * epsij_2)))

            atmTyp.sig_2.cur = sigij_2

        # update sig_2
        elif epsi_2.rng == 0.0:
            sigij_1 = cr.getSigma(sigi_1, sigj_1, alpha)
            epsij_1 = cr.getEpsilon(epsi_1, epsj_1, sigi_1, sigj_1, alpha)
            sigij_2 = cr.getSigma(sigi_2, sigj_2, alpha)

            c6 = 4.0 * epsij_1 * math.exp(6.0 * math.log(sigij_1))
            epsij_2 = c6 / (4 * math.exp(6 * math.log(sigij_2)))

            atmTyp.eps_2.cur = epsij_2


def updateEffectivePrm(cr, crPrms, scl_sig_NEI, scl_eps_NEI, eem, kappa, lam, matrix, molecules, atomTypes):
    """Updated values of effective parameters for each molecule.

    :param cr: (combiningRule) Combining rule.
    :param crPrms: (dict) Parameters to do linear combination of combining rules.
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
    molecules_utils.computeCR(cr, crPrms, scl_sig_NEI, scl_eps_NEI, molecules, atomTypes, matrix)


def writePrm(atomTypes, fileName):
    """Writes optimized parameters to txt file.

    :param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
    :param fileName: (str) Name out output file.
    :return:
    """

    with open(fileName, 'w') as out:
        out.write('{0:<5} {1:4} {2:11} {3:13} {4:5} {5:13} {6:5} {7:8} {8:6} {9:8} {10:5} {11:13} {12:5} {13:13} {14}\n'
                  .format('#iac', 'typ', 'nam', 'sig', 'rng', 'eps', 'rng', 'hrd', 'rng', 'eln', 'rng', 'sig_2', 'rng',
                          'eps_2', 'rng'))

        for iac in atomTypes:
            atomTyp = atomTypes[iac]
            typ = atomTyp.typ
            nam = atomTyp.nam
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

            out.write('{0:<5} {1:4} {2:10} {3:13.6E} {4:5.2f} {5:13.6E} {6:5.2f} {7:8.5f} {8:5.2f} {9:9.5f} {10:5.2f} '
                      '{11:13.6E} {12:5.2f} {13:13.6E} {14:5.2f}\n'
                      .format(iac, typ, nam, sig, sig_rng, eps, eps_rng, hrd, hrd_rng, eln, eln_rng, sig_2, sig_rng_2,
                              eps_2, eps_rng_2))


def write_cr_prm(crPrms, fileName):
    with open(fileName, 'w') as out:
        out.write('{:6} {:4} {:4}\n'.format("#name", "val", "rng"))

        for key, val in crPrms.items():
            out.write('{:5} {:4} {:4}\n'.format(key, val['val'], val['rng']))
