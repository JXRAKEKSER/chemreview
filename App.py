# импорты фреймворка
import tkinter.ttk as tk
import tkinter
from tkinter.filedialog import askopenfile, asksaveasfile
from tkinter.messagebox import showerror
from PIL import Image, ImageTk

# импорты библиотек для работы с веществами
from rdkit import Chem
from rdkit.Chem import Draw

# импорты утилитарных функций
from utils.utils import getDataFrame, createCsvFile, getDataFrameFromDict

# модель и её хелперы
from predictNN.PredictModel import PredicModel

# сервисы
from services.InputPredictService import InputPredictService
from services.FilePredictService import FilePredictService
from services.HistoryService import HistoryService

# инфраструктура
import os
from pymongo import *

# компоненты
from components.RecordTable import RecordTable

class App(tk.Frame):
    def __init__(self, master=None, sizes='1200x600'):
        super().__init__(master)
        self.pack()
        self.master.title('ChemReview')
        self.master.geometry(sizes)
        self._state = {
            'formState': {
                'fileInput': None,
                'smilesInput': None,
            },
            'predictionValue': None,
        }
        self._dbInstanse = MongoClient(host=os.getenv('DB_HOST'), port=int(os.getenv('DB_PORT')))[os.getenv('DB_NAME')]

        

    def render(self):
        self._notebookWidget = tk.Notebook(master=self.master, width=900)
        self._notebookWidget.pack()

        self._predictionFrame = tk.Frame(self._notebookWidget)
        self._historyFrame = tk.Frame(self._notebookWidget)

        self._notebookWidget.add(self._predictionFrame, text='Расчёт')
        self._notebookWidget.add(self._historyFrame, text='История')
        # инициализация виджетов и их слушателей для окна расчёта
        self._smilesInput = tk.Entry(self._predictionFrame, width=60, justify=tkinter.CENTER)
        self._smilesInput.grid(row=0, column=0)

        self._inputCaption = tk.Label(self._predictionFrame, text='Неправильная строка SMILES')
        self._fileButtonCaption = tk.Label(self._predictionFrame, text='Неверные столбцы в файле')
        self._submitButton = tk.Button(self._predictionFrame, state=tkinter.DISABLED, text='Start Review', command=self._handleSumbitForm)
        self._submitButton.grid(row=0, column=1, pady=(10, 5))

        self._fileButton = tk.Button(self._predictionFrame, text='Открыть файл', command=self._handleOpenFile)
        self._fileButton.grid(row=1, column=1, pady=(0, 10))

        self._answerPredict = tk.Label(self._predictionFrame, text='Predicted value: ')
        self._answerPredict.grid(row=2, column=2)
        self._answerPredict.grid_forget()

        self._smilesInput.bind('<KeyRelease>', self._handleValidateSmiles)

        # инициализация виджетов и их слушателей для окна истории

        historyService = HistoryService(self._dbInstanse)
        self._historyRecords = historyService.getAllRecords()
        RecordTable(self._historyFrame,
                    recordsList=self._historyRecords,
                    childComponentService=historyService).render()


    def _handleSumbitForm(self):
        self._answerPredict.grid_forget()
        predictModel = PredicModel()

        if self._state['formState']['smilesInput']:
            predictService = InputPredictService(self._state['formState']['smilesInput'], predictModel.predict)
            self._state['predictionValue'] = predictService.process()
            
            if self._state['predictionValue'] is None:
                showerror('Ошибка преобразования', 'Невозможно распознать молекулу')
                return
            self._answerPredict.config(text=f'Predicted value: {str(self._state["predictionValue"])}')
            mol = Chem.MolFromSmiles(self._state['formState']['smilesInput'])
            molImage = Draw.MolToImage(mol, size=(300, 300))
            img = ImageTk.PhotoImage(molImage)
            labl = tk.Label(self._predictionFrame, image=img)
            labl.image = img
            labl.grid(row=3, column=0)
            self._answerPredict.grid(row=6, column=0)
            historyServise = HistoryService(self._dbInstanse)
            historyServise.addRecord({
                'prediction': self._state['predictionValue'],
                'drug': self._state['formState']['smilesInput'] })
        else:
            drugList = list(self._state['formState']['fileInput']['Drug'])
            predictService = FilePredictService(drugList, predictModel.predict)
        
            predictedDict, incorrectDict = predictService.process()
            path = asksaveasfile(defaultextension='.csv')
            predictedDataFrame = createCsvFile(getDataFrameFromDict(predictedDict), path.name)

    def _handleOpenFile(self):
        filePath = askopenfile(filetypes=[("CSV Files", "*.csv")])
        if not filePath:
            return
        try:
            self._state['formState']['fileInput'] = getDataFrame(filePath.name)
            
            self._fileButtonCaption.grid(row=2, column=1)
            self._fileButtonCaption.config(text=filePath.name)

            self._submitButton['state'] = tkinter.ACTIVE
        except Exception as error:
            showerror('Ошибка чтения файла', message=error)

    def _handleValidateSmiles(self, event):
        self._state['formState']['smilesInput'] = event.widget.get()
        if ' ' in self._state['formState']['smilesInput']:
            self._submitButton['state'] = tkinter.DISABLED
            self._inputCaption.grid(row=1, column=0)
            return

        mol = Chem.MolFromSmiles(self._state['formState']['smilesInput'])
        if mol is None or self._state['formState']['smilesInput'] == '':
            self._submitButton['state'] = tkinter.DISABLED
            if self._state['formState']['smilesInput'] == '':
                self._inputCaption.grid_forget()
                self._fileButton['state'] = tkinter.ACTIVE
            else:
                self._inputCaption.grid(row=1, column=0)
                self._fileButton['state'] = tkinter.DISABLED
        else:
            self._submitButton['state'] = tkinter.ACTIVE
            self._inputCaption.grid_forget()
            self._fileButton['state'] = tkinter.DISABLED
