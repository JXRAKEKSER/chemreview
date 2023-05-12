from importlib.util import module_for_loader
# сервисы
from services.InputPredictService import InputPredictService
from services.FilePredictService import FilePredictService

import tkinter as tk
from predictNN.PredictModel import PredicModel
from utils.utils import getDataFrame, getDataFrameFromDict, createCsvFile
from fileinput import filename
from tkinter.filedialog import askopenfile, asksaveasfile
from tkinter.messagebox import showerror
from tkinter import Canvas, PhotoImage
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageTk

inputCaption = None
dataFrame = None
img = None
formState = {
    'fileInput' : None,
    'smilesInput' : None
}

def validateSmiles(event):
    formState['smilesInput'] = event.widget.get()
    if ' ' in formState['smilesInput']:
        submitButton['state'] = tk.DISABLED
        inputCaption.grid(row=1, column=0)
        return

    mol = Chem.MolFromSmiles(formState['smilesInput'])
    if mol is None or formState['smilesInput'] == '':
        submitButton['state'] = tk.DISABLED
        if formState['smilesInput'] == '':
            inputCaption.grid_forget()
            fileButton['state'] = tk.ACTIVE
        else:
            inputCaption.grid(row=1, column=0)
            fileButton['state'] = tk.DISABLED
    else:
        submitButton['state'] = tk.ACTIVE
        inputCaption.grid_forget()
        fileButton['state'] = tk.DISABLED

def submitForm():
    predictModel = PredicModel()
    if formState['smilesInput']:
        predictService = InputPredictService(formState['smilesInput'], predictModel.predict)
        predictedValue = predictService.process()
        #predictedValue = predictModel.predict(formState['smilesInput'])
        if predictedValue is None:
            showerror('Ошибка преобразования', 'Невозможно распознать молекулу')
            return
        answerPredict.config(text=f'Predicted value: {str(predictedValue)}')
        mol = Chem.MolFromSmiles(formState['smilesInput'])
        molImage = Draw.MolToImage(mol, size=(300, 300))
        img = ImageTk.PhotoImage(molImage)
        labl = tk.Label(window, image=img)
        labl.image = img
        labl.grid(row=3, column=0)
        

    else:
        drugList = list(formState['fileInput']['Drug'])
        predictService = FilePredictService(drugList, predictModel.predict)
        
        predictedDict, incorrectDict = predictService.process()
        path = asksaveasfile(defaultextension='.csv')
        predictedDataFrame = createCsvFile(getDataFrameFromDict(predictedDict), path.name)

def openFile():
    filePath = askopenfile(filetypes=[("CSV Files", "*.csv")])
    if not filePath:
        return
    try:
        formState['fileInput'] = getDataFrame(filePath.name)
        fileButtonCaption.grid(row=2, column=1)
        fileButtonCaption.config(text=filePath.name)
        submitButton['state'] = tk.ACTIVE
    except Exception as error:
        showerror('Ошибка чтения файла', message=error)
    

if __name__ == '__main__':
    window = tk.Tk()
    window.title('ChemReview')
    window.geometry('800x600')

    

    smilesInput = tk.Entry(window, width=60, justify=tk.CENTER)
    smilesInput.grid(row=0, column=0)

    inputCaption = tk.Label(window, text='Неправильная строка SMILES')
    fileButtonCaption = tk.Label(window, text='Неверные столбцы в файле')

    submitButton = tk.Button(window, state=tk.DISABLED, text='Start Review', command=submitForm)
    submitButton.grid(row=0, column=1)

    fileButton = tk.Button(window, text='Открыть файл', command=openFile)
    fileButton.grid(row=1, column=1)


    answerPredict = tk.Label(window, text='Predicted value: ')
    answerPredict.grid(row=6, column=0)

    smilesInput.bind('<KeyRelease>', validateSmiles)

    window.mainloop()