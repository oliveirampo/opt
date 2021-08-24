"""Module for sensitivity object.

Classes:
    Sensitivity
"""

class Sensitivity():
    """Base class for sensitivity object.

    Attributes:
        _df: (pandas DataFrame) Table with sensitivities.
        _properties: (list) List of properties.

    Methods:
        addPrmInfo(self, typ, idx1, idx2, nam1, nam2, val)
        setPrmName(self)
        addDerivative(self, propCod, avg, maxDev)
        printAllDerivatives(self, propCod)
        getDerivative(self, prmName, propCod)
        writeResFile(self, cod, out)
    """

    def __init__(self, df):
        """Constructs all the necessary attributes for the sensitivity.

        :param df: (pandas DataFrame) Table with sensitivities.
        """

        self._df = df.copy()
        self._properties = []

    def __str__(self):
        print(self._df.head(2))
        s = '({}, {})'.format(self._df.shape[0], self._df.shape[1])
        return s

    def addPrmInfo(self, typ, idx1, idx2, nam1, nam2, val):
        """Adds information about parameters to DataFrame
        and adds column with parameter name.

        :param typ: (pandas Series) Parameter type: C06_NRM, C12_NRM, C06_NEI, C12_NEI or CHG_ATM.
        :param idx1: (pandas Series) Index of atom 1.
        :param idx2: (pandas Series) Index of atom 2.
        :param nam1: (pandas Series) Name of atom 1.
        :param nam2: (pandas Series) Name of atom 2.
        :param val: (pandas Series) Parameter value.
        """

        self._df['typ'] = typ
        self._df['idx1'] = idx1
        self._df['idx2'] = idx2
        self._df['nam1'] = nam1
        self._df['nam2'] = nam2
        self._df['val'] = val

        self.setPrmName()

        self._df.set_index('prmName', verify_integrity=True, inplace=True)

    def setPrmName(self):
        """Creates column with new name of parameters."""

        self._df['prmName'] = self._df[['typ', 'idx1', 'idx2']].apply(lambda x: '_'.join(x.map(str)), axis=1)

    def addDerivative(self, propCod, avg, maxDev):
        """Adds derivatives to DataFrame.

        :param propCod: (str) Property code.
        :param avg: (numpy.ndarray) Average of derivatives per parameter.
        :param maxDev: (numpy.ndarray) Maximum deviation of mean derivative per parameter.
        """

        self._properties.append(propCod)
        self._df[propCod + '_avg'] = avg
        self._df[propCod + '_maxDev'] = maxDev

    def printAllDerivatives(self, propCod):
        """Prints DataFrame.

        :param propCod: (str) Property code.
        """

        print(self._df[['typ', 'idx1', 'idx2', propCod + '_avg']])

    def getDerivative(self, prmName, propCod):
        """Returns the derivative of property with respect to parameter for the current molecule.

        :param prmName: (str) Parameter name.
        :param propCod: (str) Property code.
        :return: (float) Derivative of property given parameter for the current molecule.
        """

        col = propCod + '_avg'
        val = self._df.loc[prmName, col]
        return val

    def writeResFile(self, cod, out):
        """Writes average and maximal deviation of derivatives to txt file.

        :param cod: (str) Molecule code.
        :param out: (output object) Output file.
        """

        for idx, row in self._df.iterrows():
            typ = row['typ']
            idx1 = row['idx1']
            idx2 = row['idx2']
            nam1 = row['nam1']
            nam2 = row['nam2']

            s = ''
            for letter in self._properties:
                avg = row[letter + '_avg']
                maxDev = row[letter + '_maxDev']
                err = 0
                if avg != 0.0:
                    err = 100.0 * maxDev / avg

                s = '{} {:>13.3e} ( {:10.3e} ) [ {:>7.2f} ]'.format(s, avg, maxDev, err)

            val = row['val']
            out.write('{:6} {:8} {:3} {:3} {:5} {:5} {:>10.3e} '
            .format(cod, typ, idx1, idx2, nam1, nam2, val))
            out.write(s)
            out.write('\n')
