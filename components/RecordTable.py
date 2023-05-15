import tkinter
from tkinter import ttk
from components.TableRow import TableRow

class RecordTable(ttk.Frame):
    def __init__(self,
                 master=None,
                 recordsList=[],
                 childComponentService=None):
        super().__init__(master=master)
        self.pack()
        self._recordsList = recordsList
        self._childComponentService = childComponentService

    def render(self):
        if (len(self._recordsList) == 0):
            emptyCaption = ttk.Label(self,
                                     text='История пуста',
                                     justify=tkinter.CENTER)
            emptyCaption.pack()
            return
        for record in self._recordsList:
            TableRow(self,
                     state=record,
                     service=self._childComponentService).render()
        
