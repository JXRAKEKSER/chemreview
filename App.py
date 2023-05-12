import tkinter.ttk as tk

class App(tk.Frame):
    def __init__(self, window=None, sizes='800x600'):
        super().__init__(window)
        self.pack()
        self.master.title('ChemReview')
        self.master.geometry(sizes)

    def render(self):
        self._smilesInput = tk.Entry(self, width=60, justify=tk.CENTER)
        self._smilesInput.grid(row=0, column=0)

        self._inputCaption = tk.Label(self, text='Неправильная строка SMILES')
        self._fileButtonCaption = tk.Label(self, text='Неверные столбцы в файле')
        self._submitButton = tk.Button(self, state=tk.DISABLED, text='Start Review', command=submitForm)

        self._fileButton = tk.Button(self, text='Открыть файл', command=openFile)
        self._fileButton.grid(row=1, column=1)

        self._answerPredict = tk.Label(window, text='Predicted value: ')
        self._answerPredict.grid(row=2, column=2)


        