import pandas as pd


class Sensitivity():
    def __init__(self, df):
        self._df = df.copy()

    def __str__(self):
        s = '({}, {})'.format(self._df.shape[0], self._df.shape[1])
        return s

    def addPrmInfo(self, typ, idx1, idx2):
        self._df['typ'] = typ
        self._df['idx1'] = idx1
        self._df['idx2'] = idx2

    def addDerivative(self, letter, avg, maxDev):
        self._df[letter + '_avg'] = avg
        self._df[letter + '_maxDev'] = maxDev

    def getAllDerivatives(self, letter):
        print(self._df[['typ', 'idx1', 'idx2', letter + '_avg']])



