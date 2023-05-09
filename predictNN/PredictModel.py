from tkinter import N
import numpy
import torch
import pandas
from catboost import CatBoostRegressor, Pool
from dgllife.model import load_pretrained
from rdkit import Chem
from dgllife.utils import AttentiveFPAtomFeaturizer, AttentiveFPBondFeaturizer, smiles_to_bigraph

class PredicModel:
    def __init__(self):
        self._mpnnModel = load_pretrained('MPNN_attentivefp_PCBA')
    
    def predict(self, dataList):
        if len(dataList) != 1:
            predictedDataframe = pandas.DataFrame({'Drug': [], 'Prediction' : []})
            incorrectDataframe = pandas.DataFrame({'Drug' : []})
            for drug in dataList:
                featuresDict = self._getDictFeatures(drug)
                if featuresDict is None:
                    incorrectDataframe = incorrectDataframe.append({'Drug' : drug }, ignore_index=True)
                else:

                    dfToPredict = self._dictToDf(featuresDict)
                    catBoost = CatBoostRegressor()
                    catBoostModel = catBoost.load_model('predictModel.cbm')
                    preditedValue = catBoostModel.predict(Pool(dfToPredict))[0]
                    predictedDataframe = predictedDataframe.append({'Drug': drug, 'Prediction' : preditedValue}, ignore_index=True)
            
            predictedDataframe.to_csv('predicted.csv')
            incorrectDataframe.to_csv('incorrectMolecles.csv')
        else:
            inputSmiles = dataList[0]
            print(f'inputSmiles: {inputSmiles}')
            featuresDict = self._getDictFeatures(inputSmiles)
            if featuresDict is None:
                return None
            dfToPredict = self._dictToDf(featuresDict)
            catBoost = CatBoostRegressor()
            catBoostModel = catBoost.load_model('predictModel.cbm')
            return catBoostModel.predict(Pool(dfToPredict))[0]

    def _extractEmbedding(self, smilesString):
        drugMol = Chem.MolFromSmiles(smilesString)

        atomFeaturizer = AttentiveFPAtomFeaturizer()
        bondFeaturizer = AttentiveFPBondFeaturizer()

        drugAtomsFeatures = atomFeaturizer(drugMol)
        drugBondFeatures = bondFeaturizer(drugMol)
        if(len(drugBondFeatures) == 0):
            return None

        temporaryList = []
        for edge in drugBondFeatures['e'].numpy():
            temporaryList.append(numpy.append(edge, [0.]))
  
        drugBondFeaturesExtend = torch.tensor(temporaryList, dtype=torch.float)

        drugGraph = smiles_to_bigraph(smilesString)
        self._mpnnModel.eval()
        embedding = self._mpnnModel(drugGraph, drugAtomsFeatures['h'], drugBondFeaturesExtend)

        return embedding
    
    def _getDictFeatures(self, smiles):
        embedding = self._extractEmbedding(smiles)
        if(embedding != None):
            featuresDict = {}
            for feature in embedding[0].detach().numpy():
                index = len(featuresDict)
                featuresDict[f'feature_{index}'] = feature
            return featuresDict
        return None

    def _dictToDf(self, featuresDict):
        moleculeFeaturesDF = None
        for key in featuresDict:
            featuresDict[key] = [featuresDict.get(key)]
            moleculeFeaturesDF = pandas.DataFrame.from_dict(featuresDict)
        return moleculeFeaturesDF
