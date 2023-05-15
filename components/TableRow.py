from tkinter import ttk

class TableRow(ttk.Frame):
    
    def __init__(self, master=None, state={}, service=None):
        super().__init__(master=master)
        self.grid()
        self._state = state
        self._service = service

    def render(self):
        drug = self._state['drug']
        prediction = self._state['prediction']

        self._drugLabel = ttk.Label(self, text=drug)
        self._drugLabel.grid(row=0,
                             column=0,
                             padx=(10, 10),
                             pady=(5, 5))

        self._predictionLabel = ttk.Label(self, text=prediction)
        self._predictionLabel.grid(row=0, column=1, padx=(10, 10), pady=(5, 5))

        self._deleteRecordButton = ttk.Button(self, text=' X ',
                                              command=self._handleDeleteRecord)
        self._deleteRecordButton.grid(row=0, column=2)

    def _handleDeleteRecord(self):
        self._service.removeRecordById(self._state['_id'])
        self.destroy()
