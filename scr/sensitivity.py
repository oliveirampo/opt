import pandas as pd


class Sensitivity():
    def __init__(self, df):
        self._df = df.copy()
        self._properties = []

    def __str__(self):
        print(self._df.head(2))
        s = '({}, {})'.format(self._df.shape[0], self._df.shape[1])
        return s

    def addPrmInfo(self, typ, idx1, idx2, nam1, nam2, val):
        self._df['typ'] = typ
        self._df['idx1'] = idx1
        self._df['idx2'] = idx2
        self._df['nam1'] = nam1
        self._df['nam2'] = nam2
        self._df['val'] = val

        self.setPrmName()

        self._df.set_index('prmName', verify_integrity=True, inplace=True)

    def setPrmName(self):
        self._df['prmName'] = self._df[['typ', 'idx1', 'idx2']].apply(lambda x: '_'.join(x.map(str)), axis=1)

    def addDerivative(self, propCod, avg, maxDev):
        self._properties.append(propCod)
        self._df[propCod + '_avg'] = avg
        self._df[propCod + '_maxDev'] = maxDev

    def getAllDerivatives(self, propCod):
        print(self._df[['typ', 'idx1', 'idx2', propCod + '_avg']])

    def getDerivative(self, prmName, propCod):
        col = propCod + '_avg'
        val = self._df.loc[prmName, col]
        return val

    # write file with average of derivatives
    def writeResFile(self, cod, out):
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
            out.write('{0:6} {1:8} {2:2} {3:2} {4:5} {5:3} {6:>10.3e} '
            .format(cod, typ, idx1, idx2, nam1, nam2, val))
            out.write(s)
            out.write('\n')
