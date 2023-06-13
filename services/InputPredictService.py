from services.AbstractPredictService import AbstractPredictService
""" 
InputPredictService - сервис, реализующий поведение AbstractPredictService
необходим для разделения поведения в случае ввода входных данных через строку ввода
 """
class InputPredictService(AbstractPredictService):
    
    def __init__(self, state, predictorCallback):
        self.state = state
        self.predictorCallback = predictorCallback

    def process(self):
        return self.predictorCallback(self.state)
    