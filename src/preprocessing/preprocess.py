import pandas as pd
import scanpy as sc
import tarfile
import glob
import gzip
import pathlib
import shutil
import os


def getProbeID(df):
    if "ProbeID" not in df.columns:
        # merge with reference file containing gene names and ENSG
        path = os.path.dirname(os.path.realpath(__file__))
        ensg_file = os.path.join(path, "ensg.txt")
        ensg_df = pd.read_csv(ensg_file, sep="\t", header=0)
        ensg_df = ensg_df.drop("ENST", axis=1).drop_duplicates("Name")

        df = pd.merge(df, ensg_df, on="Name", how="left")
        df = df.set_index(["ProbeID", "Name"]).reset_index()
    return df


def counts2expr(counts_df: pd.DataFrame) -> pd.DataFrame:
    counts_df = counts_df.set_index(["ProbeID", "Name"])
    adata = sc.AnnData(counts_df.T)
    # adata.var_names_make_unique()
    # 1e6 = counts per million (cpm) normalization
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata, base=2)
    return adata.to_df().T


def tar2rcc(tar_file: str) -> str:
    """Extracts all files from a tarfile into a new directory with '_RCC' appended

    Args:
        tar_file (str): must be a tarfile

    Returns:
        list: Path to new '_RCC' directory
    """
    if not tarfile.is_tarfile(tar_file):
        raise ValueError(f"{tar_file} is not a tarfile")

    filename = os.path.basename(os.path.abspath(tar_file))
    filename = filename.split(".")[0]
    rcc_dir = filename + "_RCC"

    with tarfile.open(tar_file) as tar:
        def is_within_directory(directory, target):
            
            abs_directory = os.path.abspath(directory)
            abs_target = os.path.abspath(target)
        
            prefix = os.path.commonprefix([abs_directory, abs_target])
            
            return prefix == abs_directory
        
        def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
        
            for member in tar.getmembers():
                member_path = os.path.join(path, member.name)
                if not is_within_directory(path, member_path):
                    raise Exception("Attempted Path Traversal in Tar File")
        
            tar.extractall(path, members, numeric_owner=numeric_owner) 
            
        
        safe_extract(tar, rcc_dir)

    gunzip_dir(rcc_dir)
    return rcc_dir


def gunzip_dir(my_dir: str) -> None:
    """Unzips any gzipped files in given directory

    Args:
        my_dir (str): path to directory for files to be gunzipped
    """
    for gzip_file in glob.glob(os.path.join(my_dir, "*.gz")):
        gunzip_file = pathlib.Path(gzip_file).stem
        gunzip_file = os.path.join(my_dir, gunzip_file)
        with gzip.open(gzip_file, "rb") as file_in, open(gunzip_file, "wb") as file_out:
            shutil.copyfileobj(file_in, file_out)
            os.remove(gzip_file)
