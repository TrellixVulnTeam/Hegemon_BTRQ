import argparse
import pandas as pd


from src.hegemon import Hegemon
from src.preprocessing.gse import GSE

parser = argparse.ArgumentParser(description="CLI For Hegemon Files")
parser.add_argument(
    "-g",
    "--gse",
    help="NCBI GEO GSE Accession ID to parse",
)
args = parser.parse_args()

def parse_gse(accessionID: str) -> dict:
    my_gse = GSE(accessionID)
    gpls = {}
    for gpl in my_gse.gpls:
        expr = my_gse.expr(gpl)
        survival = my_gse.survival(gpl)
        gpls[gpl] = {"expr": expr, "survival": survival}
    return gpls

def main(filebase: str, expr: pd.DataFrame, survival: pd.DataFrame) -> None:
    hegemon = Hegemon(expr, survival)
    idx = hegemon.idx()
    ih = hegemon.ih()

    dfs = [expr, survival, idx, ih]
    names = ["expr", "survival", "idx", "ih"]
    for df, name in zip(dfs, names):
        df.to_csv(f"{filebase}-{name}.txt", sep="\t")

if __name__ == "__main__":
    gpls = parse_gse(args.gse)
    for gpl_name, hegemon_files in gpls.items():
        expr = hegemon_files["expr"]
        survival = hegemon_files["survival"]
        filebase = f"{args.gse}-{gpl_name}"
        main(filebase=filebase, expr=expr, survival=survival)