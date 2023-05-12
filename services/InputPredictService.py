from services.AbstractPredictService import AbstractPredictService

class InputPredictService(AbstractPredictService):
    
    def __init__(self, state, predictorCallback):
        self.state = state
        self.predictorCallback = predictorCallback

    def process(self):
        return self.predictorCallback(self.state)
    