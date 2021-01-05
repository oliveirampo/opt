from scipy import optimize
import sys


import ana


def runOptimization(conf, molecules, atomTypes):
    # get simulated results
    addSimProp(conf, molecules)

    prmsToOptmize = getParametersToBeOptimized(atomTypes)
    setMinMax(conf, prmsToOptmize)

    minimize()


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

    for prmList in prmsToOptmize:
        minVal = []
        maxVal = []
        iac = []

        for prm in prmList:
            val = prm.ori * (1.0 - prm.rng * rng_scl / 100.0)
            minVal.append(val)
            val = prm.ori * (1.0 + prm.rng * rng_scl / 100.0)
            maxVal.append(val)

            prm.min = minVal
            prm.max = maxVal

            iac.append(prm.iac)

        minVal = list(set(minVal))
        maxVal = list(set(maxVal))

        if len(minVal) != 1 and len(maxVal) != 1:
            print('Symmetrical parameters with different values: {}'.format(iac))
            sys.exit(1)


def minimize():
    print('TODO::minimize()')







