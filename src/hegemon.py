from dataclasses import dataclass
import pandas as pd
from io import BytesIO


from perl_subprocess import PerlFiles

@dataclass
class Hegemon:
    expr: pd.DataFrame
    survival: pd.DataFrame

    def __post_init__(self):
        self.perl = PerlFiles(self.expr, "filler")

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
    import sys
    expr = sys.argv[1]
    expr = pd.read_csv(expr, sep="\t")
    survival = sys.argv[2]
    survival = pd.read_csv(survival, sep="\t")

    my_hegemon = Hegemon(expr, survival)
    print(my_hegemon.perl.bv())
        
