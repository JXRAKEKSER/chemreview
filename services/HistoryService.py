""" from pymongo import * """
class HistoryService():
    COLLECTION = 'Predictions'

    def __init__(self, dbInstanse):
        self._collection = dbInstanse[self.COLLECTION]
        
    
    def addRecord(self, predictionRecord):
        '''predictionRecord = {
          prediction: Float64,
          drug: String
        }'''
        self._collection.insert_one(predictionRecord)
    
    def getAllRecords(self):
        return list(self._collection.find({}))
    
    def removeRecordById(self, id):
        self._collection.delete_one({'_id': id})


""" client = MongoClient('localhost', 27017)

db = client['ChemReview']
HistoryService(db) """

""" predictionCollection = db['Predictions']


def addRecord(collection, predictionRecord):
    collection.insert_one(predictionRecord)

def removeRecordById(collection, id):
    collection.delete_one({'_id': id})

def removeAllRecords(collection):
    collection.delete_many({})

print(list(predictionCollection.find())) """