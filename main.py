import argparse

from src.preprocessing.gse import GSE
from src.hegemon import Hegemon

parser = argparse.ArgumentParser(description="CLI For GSE ID")
parser.add_argument(
    "-g",
    "--gse",
    help="NCBI GEO GSE Accession ID to parse",
)
args = parser.parse_args()

def parse_gse(accessionID: str) -> None:
    my_gse = GSE(accessionID)
    for gpl_name in my_gse.gpls:
        expr = my_gse.expr(gpl_name)
        survival = my_gse.survival(gpl_name)
        filebase = f"{accessionID}-{gpl_name}"
        Hegemon(expr=expr, survival=survival, filebase=filebase)

if __name__ == "__main__":
    gpls = parse_gse(args.gse)
    