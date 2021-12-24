from dataclasses import dataclass
import pandas as pd
from io import BytesIO

@dataclass
class Hegemon:
    expr: pd.DataFrame
    survival: pd.DataFrame

    def idx(self) -> pd.DataFrame:
        df_bin = BytesIO()
        self.expr.to_csv(df_bin, index=False, mode="wb", sep="\t")
        df_bin.seek(0)
        
        ptr = []
        pos = 0
        for line in df_bin:
            if pos == 0:
                pos += len(line)
            else:
                ptr.append(pos)
                pos += len(line)
        
        idx = self.expr[["ProbeID", "Name"]]
        idx.insert(1, "Ptr", ptr)
        return idx

    def ih(self) -> pd.DataFrame:
        ih_df = self.survival.rename(columns={"c title": "Title"})
        ih_df = ih_df[["ArrayID", "Title"]]
        ih_df.insert(1, "ArrayHeader", ih_df["ArrayID"])
        return ih_df

if __name__ == "__main__":
    from preprocessing.gse import GSE
    my_gse = GSE("GSE129800")
    expr = my_gse.expr(my_gse.gpls[0])
    survival = my_gse.survival(my_gse.gpls[0])
    my_hegemon = Hegemon(expr, survival)
    print(my_hegemon.ih())
        
