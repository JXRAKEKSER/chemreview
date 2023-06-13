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
    
    def predict(self, smiles):

        featuresDict = self._getDictFeatures(smiles)
        if featuresDict is None:
            return None
            
        dfToPredict = self._dictToDf(featuresDict)
        catBoost = CatBoostRegressor()
        catBoostModel = catBoost.load_model('predictModel.cbm')
        preditedValue = catBoostModel.predict(Pool(dfToPredict))[0]

        return preditedValue

    def _extractEmbedding(self, smiles):
        # получение объекта молекулы Rdkit
        drugMol = Chem.MolFromSmiles(smiles)

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

        # получение графового представления на основе smiles
        drugGraph = smiles_to_bigraph(smiles)
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
