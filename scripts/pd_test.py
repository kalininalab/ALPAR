import pandas as pd
import numpy as np

df_dict = {"a": {"pa": 1, "wf" : 0, "wag": 1}, "b": {"pa": 0, "wf" : 0, "wag": 0}, "c": {"pa": 1, "wf" : 1, "wag": 1}}

df = pd.DataFrame.from_dict(df_dict, orient='index')

df.index.name = 'panacota_renamed'
df.reset_index(inplace=True)

print(df)