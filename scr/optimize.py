from scipy import optimize
import sys


import molecules_utils
import ana


def runOptimization(conf, molecules, atomTypes):
    molecules_utils.computeEEM(conf.eem, molecules, atomTypes)
    updateEffectivePrm(conf.cr, conf.eem, conf.matrix, molecules, atomTypes)

    updateOirginalParameterValues(molecules)

    # get simulated results
    addSimProp(conf, molecules)

    prmsToOptmize = getParametersToBeOptimized(atomTypes)
    # set min/max values and return initial values
    init = setMinMax(conf, prmsToOptmize)

    # prm1 = atomTypes[12].sig
    # print('prm1: ', prm1.iac, prm1.cur)
    # testParameter(prmsToOptmize)
    # print('prm1: ', prm1.iac, prm1.cur)

    minimize(conf, init, prmsToOptmize, atomTypes, molecules)


def testParameter(prmsToOptmize):
    prm2 = prmsToOptmize[0][0]
    print('prm2: ', prm2.iac, prm2.cur)

    prm2.cur = 1.0
    print('prm2: ', prm2.iac, prm2.cur)


def updateOirginalParameterValues(molecules):
    for cod in molecules:
        mol = molecules[cod]
        parameters = mol.parameters
        for prm in parameters:
            prm.ori = prm.cur


def addSimProp(conf, molecules):
    for cod in molecules:
        mol = molecules[cod]

        # add running average of each job
        ana.addSimprop(conf, mol)

        # add sensitivity matrix
        ana.addSens(conf, mol)


def getParametersToBeOptimized(atomTypes):
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


def minimize(conf, init, prmsToOptmize, atomTypes, molecules):
    eem = conf.eem
    matrix = conf.matrix
    cr = conf.cr
    opt_nit = conf.opt_nit
    scl_sig_NEI = conf.scl_sig_NEI
    scl_eps_NEI = conf.scl_eps_NEI

    minimizer_kwargs = (prmsToOptmize, atomTypes, cr, eem, matrix, molecules, {'Nfeval': 0})
    option = {'disp': True, 'maxiter': opt_nit}

    # x_fnc = targetFunction(init, prmsToOptmize, atomTypes, cr, eem, matrix, molecules, {'Nfeval':0})
    ret = optimize.minimize(targetFunction, init, args=minimizer_kwargs, options=option, method='Nelder-Mead')
    print(ret)


def targetFunction(init, prmsToOptmize, atomTypes, cr, eem, matrix, molecules, info):
    updatePrm(init, prmsToOptmize)
    updateEffectivePrm(cr, eem, matrix, molecules, atomTypes)

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

    if info['Nfeval'] % 100 == 0:
        # TODO - write value of target function
        print(x_fnc)
        # sys.exit(123)

    info['Nfeval'] += 1
    return x_fnc


# adjust current value to [min, max]
def updatePrm(init, prmsToOptmize):
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


def updateEffectivePrm(cr, eem, matrix, molecules, atomTypes):
    molecules_utils.computeEEM(eem, molecules, atomTypes)
    molecules_utils.computeCR(cr, molecules, atomTypes, matrix)
