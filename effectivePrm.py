from abc import ABC, abstractmethod


class EffectivePrm(ABC):
    @abstractmethod
    def updateValue(self):
        pass


class EffectivePrmLJ(EffectivePrm):
    def updateValue(self):
        print('TODO')


class EffecivePrmLJC12(EffectivePrmLJ):
    def updateValue(self):
        print('TODO-C12')


class EffecivePrmLJC06(EffectivePrmLJ):
    def updateValue(self):
        print('TODO-C06')
