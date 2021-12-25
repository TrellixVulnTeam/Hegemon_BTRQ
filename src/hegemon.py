from dataclasses import dataclass
import pandas as pd
from io import BytesIO, StringIO
import subprocess
import os


@dataclass
class Hegemon:
    expr: pd.DataFrame
    survival: pd.DataFrame
    filebase: str

    def __post_init__(self):
        self.survival.to_csv(self.filebase+"-expr.txt", sep="\t", index=False)
        self.expr.to_csv(self.filebase+"-expr.txt", sep="\t", index=False)
        self.idx().to_csv(self.filebase+"-idx.txt", sep="\t", index=False)
        self.ih().to_csv(self.filebase+"-ih.txt", sep="\t", index=False)
        self.thr().to_csv(self.filebase+"-thr.txt", sep="\t", index=False)
        self.vinfo().to_csv(self.filebase+"-vinfo.txt", sep="\t", index=False)
        self.info().to_csv(self.filebase+"-info.txt", sep="\t", index=False)
        self.bv().to_csv(self.filebase+"-bv.txt", sep="\t", index=False)
        self.conf_file()

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
        ih_df = self.survival.rename(columns={"c title": "Title", "ArrayId": "ArrayID"})
        ih_df = ih_df[["ArrayID", "Title"]]
        ih_df.insert(1, "ArrayHeader", ih_df["ArrayID"])
        return ih_df
    
    def thr(self) -> pd.DataFrame: # parses expr file    
        result = subprocess.run(["perl", "-I", 
                                "/booleanfs/sahoo/scripts", 
                                "/booleanfs/sahoo/scripts/absoluteInfo.pl", "thr", 
                                self.filebase+"-expr.txt", 
                                "2", "70000", "0.5"], 
                                capture_output=True,
                                text=True)
        thr_df = pd.read_csv(StringIO(result.stdout), sep="\t")
        return thr_df

    def info(self) -> pd.DataFrame: # parses thr file
        result = subprocess.run(["perl", "-I",
                                 "/booleanfs/sahoo/scripts",
                                 "/booleanfs/sahoo/scripts/hegemonutils.pl",
                                 "Info",
                                 self.filebase], 
                                 capture_output=True,
                                 text=True)
        info_df = pd.read_csv(StringIO(result.stdout), sep="\t")
        return info_df

    def vinfo(self) -> pd.DataFrame: # parses thr file
        result = subprocess.run(["perl", "-I",
                                "/booleanfs/sahoo/scripts",
                                "/booleanfs/sahoo/scripts/hegemonutils.pl",
                                "VInfo",
                                self.filebase], 
                                capture_output=True,
                                text=True)
        vinfo_df = pd.read_csv(StringIO(result.stdout), sep="\t")
        return vinfo_df
        
    def bv(self) -> pd.DataFrame: # parses thr file and idx file
        result = subprocess.run(["perl", "-I",
                                "/booleanfs/sahoo/scripts",
                                "/booleanfs/sahoo/scripts/hegemonutils.pl",
                                "bv",
                                self.filebase], 
                                capture_output=True,
                                text=True)
        bv_df = pd.read_csv(StringIO(result.stdout), sep="\t")
        return bv_df
    
    def conf_file(self) -> None:
        export = os.path.join(os.getcwd(), "conf_file.txt")
        with open(export, "w") as file_out:
            file_out.write("[]\n")
            file_out.write("name=\n")

            names = ["expr", "index", "survival", "indexHeader", "info"]
            types = ["expr", "idx", "survival", "ih", "info"]
            for name, type in zip(names, types):
                file = f"{self.filebase}-{type}.txt"
                filepath = os.path.join(os.getcwd(), file)
                file_out.write(f"{name}={filepath}\n")

            file_out.write("key=")

if __name__ == "__main__":
    import sys
    expr = pd.read_csv(sys.argv[1], sep="\t")
    survival = pd.read_csv(sys.argv[2], sep="\t")

    my_hegemon = Hegemon(expr, survival, filebase="filler")
        
