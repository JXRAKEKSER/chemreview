import pandas

def getDataFrame(filePath):
  try:
    df = pandas.read_csv(filePath)
    if 'Drug' not in df.columns: raise Exception('Отсутствует колонка Drug')
    return df
  except Exception as error:
    raise Exception(f'{error}')
  
def getDataFrameFromDict(dict):
  try:
    return pandas.DataFrame(dict)
  except Exception as error:
    raise Exception(f'{error}')
  
def createCsvFile(dataFrame, fileName):
  return dataFrame.to_csv(f'{fileName}')
