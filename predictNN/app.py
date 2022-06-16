from matplotlib import type1font
import pandas
import torch
import dgllife
import numpy
import tdc
import sklearn
import catboost
from tdc.multi_pred import DrugRes
from dgllife.utils import AttentiveFPAtomFeaturizer, AttentiveFPBondFeaturizer, smiles_to_bigraph, CanonicalBondFeaturizer
from dgllife.model import load_pretrained
from rdkit import Chem
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error
from catboost import CatBoostRegressor, Pool
from PredictModel import PredicModel

#data = DrugRes(name = 'GDSC1')
#df = data.get_data()

def getDataFrame(filePath):
  try:
    df = pandas.read_csv(filePath)
    if 'SMILES' not in df.columns: return None

    return df
  except Exception as error:
    return f'{error}'




def hasNull(pdList):
  for item in pdList:
    if item:
      return item
  return False

def findOptimCellLineID(dataFrame):
  cellLineIdList = dataFrame['Cell Line_ID'].unique()
  uniqueDrugNumber = len(dataFrame['Drug_ID'].unique())
  for cellItem in cellLineIdList:
    uniqueDrugsList = (dataFrame[dataFrame['Cell Line_ID'] == cellItem])['Drug_ID'].unique()
    if len(uniqueDrugsList) == uniqueDrugNumber:
      return cellItem

def extractEmbeddings(smiles):
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
  numberAtomFeatures = drugAtomsFeatures['h'].size(dim=1)
  numberBondFeatures = drugBondFeatures['e'].size(dim=1)

  drugGraph = smiles_to_bigraph(smiles)
  graphModel = load_pretrained('MPNN_attentivefp_PCBA')
  graphModel.eval()
  embedding = graphModel(drugGraph, drugAtomsFeatures['h'], drugBondFeaturesExtend)

  return embedding

def getDictFeatures(drugRow):
  smiles = drugRow['Drug']
  embedding = extractEmbeddings(smiles)
  if(embedding != None):
    featuresDict = {}
    for feature in embedding[0].detach().numpy():
      index = len(featuresDict)
      featuresDict[f'feature_{index}'] = feature
    featuresDict['Y'] = drugRow['Y']
    return featuresDict
  return None

def mape(Y_valid, Y_predict):
  return numpy.mean(numpy.abs((Y_predict - Y_valid)/Y_valid))

""" optimalCellLineID = findOptimCellLineID(df)
fullDrugDF = df[df['Cell Line_ID'] == optimalCellLineID]
fullDrugDF.drop_duplicates(subset='Drug_ID', keep='first', inplace=True)
fullDrugDF.drop(labels=['Cell Line_ID', 'Cell Line'], axis=1, inplace=True)
updatedIndexes = []
for i in range(len(fullDrugDF['Drug'])):
  updatedIndexes.append(i)
fullDrugDF.index = updatedIndexes

moleculeFeaturesDF = None
for index, drugRow in fullDrugDF.iterrows():
  featureRow = getDictFeatures(drugRow)
  if(index):
    if featureRow != None:
      moleculeFeaturesDF = moleculeFeaturesDF.append(featureRow, ignore_index=True)
  else:
    if featureRow != None:
      for key in featureRow:
        featureRow[key] = [featureRow.get(key)]
      moleculeFeaturesDF = pandas.DataFrame.from_dict(featureRow)

Y = moleculeFeaturesDF['Y']
X = moleculeFeaturesDF.drop(['Y'], axis=1)

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.1, random_state=10)




cat = CatBoostRegressor(loss_function='MAPE', iterations=100)
trainPool = Pool(data=X_train, label=Y_train)
evalPool = Pool(data=X_test, label=Y_test)
cat.fit(trainPool, use_best_model=True, verbose=0, eval_set=evalPool)
cat.save_model('predictModel.cbm', 'cbm')

Y_predict = cat.predict(Pool(data=X_test))
print(f'MAPE: {mape(Y_test, Y_predict)}')
print(f'MAE: {mean_absolute_error(Y_test, Y_predict)}') """

if __name__ == '__main__':
  data = DrugRes(name = 'GDSC1')
  df = data.get_data()
  featuresRow = getDictFeatures(df.iloc[0])

  predictModel = PredicModel(prepareDataCallback=None)
  predictModel.predict(list(df.iloc[0]))

  """ moleculeFeaturesDF = None
  for key in featuresRow:
    featuresRow[key] = [featuresRow.get(key)]
    moleculeFeaturesDF = pandas.DataFrame.from_dict(featuresRow)
  dataToPredict = moleculeFeaturesDF.drop(['Y'], axis=1)

  catboostObj = CatBoostRegressor()
  catboostModel = catboostObj.load_model('predictModel.cbm')
  
  prediction = catboostModel.predict(Pool(dataToPredict))

  print(prediction[0]) """