import pandas as pd
import GEOparse
import os
import glob

from src.preprocessing.preprocess import getProbeID

def gse130284(raw_counts_file: str) -> pd.DataFrame:
    gse130284 = GEOparse.get_GEO(geo="GSE130284", silent=True)
    os.remove(glob.glob("GSE130284*.soft*")[0])

    df = pd.read_exceL(raw_counts_file)
    df = df.rename(columns={"Id": "Name", "KNalone_2": "NKalone_2"}) # fix typo
    df = df.set_index("Name").iloc[:, :12]

    col_names = []
    # Align sample titles with GSM samples
    for col in df.columns:
        col = col.split("_")
        col = " donor ".join(col)
        col = col[:2] + " " + col[2:] + "_RNA-seq"
        if "IL" in col:
            col = col.replace("IL", "IL-")
        if "INB" in col:
            col = col.replace("INB16", "CTV-1")
        col_names.append(col)
    df.columns = col_names

    col_rename = {}
    for name, gsm in gse130284.gse.gsms.items():
        # map gsm name to sample title
        col_rename[gsm.metadata["title"][0]] = name
    df = df.rename(columns=col_rename).reset_index()
    df = getProbeID(df)
    return df