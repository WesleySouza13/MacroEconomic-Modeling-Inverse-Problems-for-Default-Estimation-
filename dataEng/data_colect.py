from request_data import Request
import pandas as pd 
import os 
lista_de_codes = [4393, 21084, 24369, 27640]
datas = []
for i in lista_de_codes:
    data_path = os.path.join('data', 'sgs_data', f'{i}.csv')
    req = Request(i, '01/01/2010', '01/01/2026').req_data()
    req.to_csv(data_path)
    datas.append(req)
df = pd.concat(datas, axis=1)
data_concat = os.path.join('data', 'raw', 'dados_completos.csv')
df.to_csv(data_concat)