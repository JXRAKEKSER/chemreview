import tkinter as tk
from predictNN.PredictModel import PredicModel
from utils.utils import getDataFrame
from fileinput import filename
from tkinter.filedialog import askopenfile
from tkinter.messagebox import showerror
from rdkit import Chem

inputCaption = None
dataFrame = None
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
        predictedValue = predictModel.predict([formState['smilesInput']])
        if predictedValue is None:
            showerror('Ошибка преобразования', 'Невозможно распознать молекулу')
            return
        answerPredict.config(text=str(predictedValue))
    else:
        drugList = list(formState['fileInput']['Drug'])
        predictedValue = predictModel.predict(drugList)
        

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
    
    
    """ if formState['fileInput'] is None:
        fileButtonCaption.config(text='Не содержит столбца SMILES')
    elif formState['fileInput'] == 'No columns to parse from file':
        fileButtonCaption.config(text=dataFrame)
    else: """
    

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
    answerPredict.grid(row=2, column=2)

    smilesInput.bind('<KeyRelease>', validateSmiles)
    window.mainloop()