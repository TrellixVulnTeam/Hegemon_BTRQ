#! /usr/bin/env python
import argparse
import pandas as pd
import os

from src.preprocessing.gse import GSE
from src.hegemon import Hegemon

parser = argparse.ArgumentParser(description="CLI For GSE ID")
parser.add_argument(
    "-g",
    "--gse",
    help="NCBI GEO GSE Accession ID to parse",
)
parser.add_argument(
    "-e",
    "--expr",
    help="expr text file tab seperated",
)
parser.add_argument(
    "-s",
    "--survival",
    help="survival text file tab seperated",
)
args = parser.parse_args()

def parse_gse(accessionID: str) -> None:
    my_gse = GSE(accessionID)
    for gpl_name in my_gse.gpls:
        expr = my_gse.expr(gpl_name)
        survival = my_gse.survival(gpl_name)
        filebase = f"{accessionID}-{gpl_name}"
        Hegemon(expr=expr, survival=survival, filebase=filebase)

def parse_expr_survival(expr: str, survival: str) -> None:
    expr = pd.read_csv(expr, sep="\t")
    survival = pd.read_csv(survival, sep="\t")
    filebase = os.path.basename(expr)[:-9]
    Hegemon(expr=expr, survival=survival, filebase=filebase)

def main():
    if args.gse != None:
        parse_gse(args.gse)
    elif args.expr != None and args.survival != None:
        parse_expr_survival(args.expr, args.survival)

if __name__ == "__main__":
    main()
    
