from services.AbstractPredictService import AbstractPredictService

class FilePredictService(AbstractPredictService):
   
  def __init__(self, state, predictorCallback):
    self.state = state
    self.predictorCallback = predictorCallback

  def process(self):
    predictedData = {'Drug': [], 'Prediction' : []}
    incorrectDrugs = {'Drug' : []}
    for drug in self.state:
      predicedValue = self.predictorCallback(drug)
      if predicedValue is None:
        incorrectDrugs['Drug'].append(drug)
      else:
        predictedData['Drug'].append(drug)
        predictedData['Prediction'].append(predicedValue)
    
    return ( predictedData, incorrectDrugs )