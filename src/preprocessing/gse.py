from dataclasses import dataclass
import GEOparse
import pandas as pd
import numpy as np
import re
import os
import glob

from .preprocess import getProbeID


@dataclass
class GSE:
    accessionID: str

    def __post_init__(self):
        # Add GEOparse GSE object as attribute
        self.gse = GEOparse.get_GEO(geo=self.accessionID, silent=True)
        self.gpls = list(self.gse.gpls.keys())

        # remove downloaded soft file
        os.remove(glob.glob(f"{self.accessionID}*.soft*")[0])

    def expr(self, gpl_name: str, takeLog: bool=True) -> pd.DataFrame:
        # confirm that gsm correlates to called gpl
        for name, gsm in self.gse.gsms.items():
            if gsm.metadata["platform_id"][0] != gpl_name:
                continue
            
            # collect count data
            gsm_df = gsm.table.rename(columns={"ID_REF":"Name"})
            gsm_df = gsm_df.set_index("Name")
            gsm_df.columns = [name]

            if takeLog:
                gsm_df[name] = np.where(gsm_df[name] > 0, np.log2(gsm_df[name]), 0)

            # merge each column into one dataframe
            if "expr_df" not in locals():
                expr_df = gsm_df
            else:
                expr_df = expr_df.merge(gsm_df, left_index=True, right_index=True)

        expr_df = expr_df.reset_index()
        # add ProbeID (ENSG)
        if "ProbeID" not in expr_df.columns:
            expr_df = getProbeID(expr_df)
        return expr_df

    def survival(self, gpl_name: str) -> pd.DataFrame:
        """Creates metadata information for each GSM (sample)

        Args:
            gpl (string): gpl name associated

        Returns:
            pd.DataFrame: survival dataframe including all samples and metadata
        """
        to_drop = [
            "geo_accession",
            "status",
            "date",
            "protocol",
            "proccesing",
            "data_processing",
            "contact",
            "supplementary",
            "platform_id",
            "series_id",
            "relation",
        ]

        all_metadata = {}
        for name, gsm in self.gse.gsms.items():
            # confirm gsm is associated with input gpl
            gsm_gpl = gsm.metadata["platform_id"][0]
            if gsm_gpl != gpl_name:
                continue

            # remove keys from metadata that aren't desired in survival
            metadata = gsm.metadata.copy()
            for key in gsm.metadata:
                for drop in to_drop:
                    if re.search(drop, key):
                        metadata.pop(key, None)
            all_metadata[name] = metadata

        all_metadata = pd.DataFrame(all_metadata).T
        # convert all list values into string
        all_metadata = all_metadata.applymap(lambda x: "\t".join(x), na_action="ignore")

        for column in all_metadata.columns:
            # split columns with multiple values into seperate columns
            to_merge = all_metadata[column].str.split("\t", expand=True)
            if len(to_merge.columns) > 1:
                # create column names for additional columns
                col_names = [str(i + 1) for i in range(len(to_merge.columns))]
                col_names = [column + "_" + name for name in col_names]
                to_merge.columns = col_names
            else:
                to_merge.columns = [column]

            if "df" not in locals():
                df = to_merge
            else:
                df = df.merge(to_merge, left_index=True, right_index=True)

        df.index.name = "ArrayID"

        for column in df.columns:
            # rename columns with cell value label
            if df[column].str.contains(":").all():
                value = df[column].str.extract(r"(.*:)").iloc[0, 0][:-1]
                df = df.rename(columns={column: value})
                df[value] = df[value].str.extract(r".*: (.*)")

        # add 'c ' label for use in hegemon
        df = df.rename(columns={col: "c " + col for col in df.columns})
        df = df.reset_index()
        return df

if __name__ == "__main__":
    gse = "GSE130284"
    gpl = "GPL18573"
    gse130284 = GSE(gse)
    surv = gse130284.survival(gpl)
        