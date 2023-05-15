from importlib.util import module_for_loader
# сервисы
from services.InputPredictService import InputPredictService
from services.FilePredictService import FilePredictService

# инфраструктура
from dotenv import load_dotenv
import os

import tkinter as tk
from App import App

if __name__ == '__main__':
    load_dotenv(dotenv_path=f'{os.getcwd()}/.env.dev')
    
    window = tk.Tk()
    app = App(master=window)
    app.render()

    window.mainloop()