import pandas as pd 
import requests 

class Request():
    def __init__(self, code:int, start:str, end:str):
        self.code = code
        self.start = start
        self.end = end 
    def req_data(self):
        url = f"https://api.bcb.gov.br/dados/serie/bcdata.sgs.{self.code}/dados?formato=json&dataInicial={self.start}&dataFinal={self.end}"
        response = requests.get(url=url)
        if response.status_code == 200:
            print(f'{self.code} - requisição bem sucedida!')
            print()
            print("transformando em dataframe")
            df = pd.DataFrame(response.json())
            df['data'] = pd.to_datetime(df['data'], errors='coerce', dayfirst=True)
            df = df.set_index('data')
            df['valor'] = pd.to_numeric(df['valor'], errors='coerce')
            df = df.rename(columns={'valor':self.code})
            if df is not None:
                print("tranformaçao bem sucedida!")
                return df 
            else:
                print('erro ao transformar em dataframe.')       