from dataclasses import dataclass
import pandas as pd
import subprocess
from io import StringIO

@dataclass
class PerlFiles:
    expr_file: str

    def thr(self) -> pd.DataFrame:
        result = subprocess.run(["perl", "-I", 
                                "/booleanfs/sahoo/scripts", 
                                "/booleanfs/sahoo/scripts/absoluteInfo.pl", "thr", 
                                self.expr_file, 
                                "2", "70000", "0.5"], 
                                capture_output=True,
                                text=True)
        thr_df = pd.read_csv(StringIO(result.stdout))
        return thr_df

    def info(self) -> pd.DataFrame:
        result = subprocess.run(["perl", "-I",
                                "/booleanfs/sahoo/scripts",
                                "/booleanfs/sahoo/scripts/hegemonutils.pl",
                                "Info",
                                self.expr_file[:-9]], 
                                capture_output=True,
                                text=True)
        info_df = pd.read_csv(StringIO(result.stdout))
        return info_df


    def vinfo(self) -> pd.DataFrame: #parses thr file
        result = subprocess.run(["perl", "-I",
                                "/booleanfs/sahoo/scripts",
                                "/booleanfs/sahoo/scripts/hegemonutils.pl",
                                "VInfo",
                                self.expr_file[:-9]], 
                                capture_output=True,
                                text=True)
        vinfo_df = pd.read_csv(StringIO(result.stdout))
        return vinfo_df
        
    def bv(self) -> pd.DataFrame: #parses thr file
        result = subprocess.run(["perl", "-I",
                                "/booleanfs/sahoo/scripts",
                                "/booleanfs/sahoo/scripts/hegemonutils.pl",
                                "bv",
                                self.expr_file[:-9]], 
                                capture_output=True,
                                text=True)
        bv_df = pd.read_csv(StringIO(result.stdout))
        return bv_df

if __name__ == "__main__":
    import sys
    expr_file = sys.argv[1]
    obj = PerlFiles(expr_file)
    obj.bv().to_csv(expr_file[:-9]+"-bv.txt", sep="\t")
