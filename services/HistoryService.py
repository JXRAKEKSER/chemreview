"""
HistoryService - сервис для работы с коллекцией Predictions из базы данных
"""
class HistoryService():
    COLLECTION = 'Predictions'

    def __init__(self, dbInstanse):
        self._collection = dbInstanse[self.COLLECTION]
        
    
    def addRecord(self, predictionRecord):
        '''predictionRecord = {
          prediction: Float64,
          drug: String
        }'''
        return self._collection.insert_one(predictionRecord)
    
    def getAllRecords(self):
        return list(self._collection.find({}))
    
    def removeRecordById(self, id):
        self._collection.delete_one({'_id': id})

    def removeAllRecords(self):
        self._collection.delete_many({})