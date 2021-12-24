import argparse
from typing import Dict
import pandas as pd
import subprocess
import os

from src.hegemon import Hegemon
from src.preprocessing.gse import GSE

parser = argparse.ArgumentParser(description="CLI For Hegemon Files")
parser.add_argument(
    "-g",
    "--gse",
    help="NCBI GEO GSE Accession ID to parse",
)
# parser.add_argument(
#     "-e",
#     "--expr",
#     help="Expression file",
# )
# parser.add_argument(
#     "-s",
#     "--survival",
#     help="Survival file",
# )
# parser.add_argument(
#     "-f",
#     "--filebase",
#     help="Basename for all exported files",
# )
args = parser.parse_args()

def parse_gse(accessionID: str) -> Dict:
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

def run_perl(expr_file: str) -> None:
    thr_file = expr_file[:-9] + "-thr.txt"
    subprocess.run(f"perl -I /booleanfs/sahoo/scripts /booleanfs/sahoo/scripts/absoluteInfo.pl thr {expr_file} 2 70000 0.5 > {thr_file}")

if __name__ == "__main__":
    gpls = parse_gse(args.gse)
    for gpl_name, hegemon_files in gpls.items():
        expr = hegemon_files["expr"]
        survival = hegemon_files["survival"]
        filebase = f"{args.gse}-{gpl_name}"
        os.mkdir(filebase)
        os.chdir(filebase)
        main(filebase=filebase, expr=expr, survival=survival)
        run_perl(filebase+"-expr.txt")
