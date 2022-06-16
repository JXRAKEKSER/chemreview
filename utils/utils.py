import pandas

def getDataFrame(filePath):
  try:
    df = pandas.read_csv(filePath)
    if 'Drug' not in df.columns: raise Exception('Отсутствует колонка Drug')
    return df
  except Exception as error:
    raise Exception(f'{error}')
