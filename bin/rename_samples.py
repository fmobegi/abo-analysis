#!/usr/bin/env python3

import sys
import re
import pandas as pd
from datetime import datetime


# Constants
DATE_FORMAT = "%Y_%m_%d"


def read_renaming_file(deobfuscation):
    """
    Read the Excel file from soft PCR HLA exported file
    """
    return pd.read_excel(deobfuscation, index_col=None, na_values=["NA"], usecols="A,C")


def preprocess_renaming_file(renaming_file):
    """
    Preprocess the renaming file by renaming columns, extracting the first value from the comma-separated Grid_number column,
    converting Grid_number column to string type, and filtering out rows where Grid_number starts with a letter
    """
    renaming_file = renaming_file.rename(
        columns={"Acc#": "Sample ID", "Patient Name": "Grid_number"}
    )
    renaming_file["Grid_number"] = renaming_file["Grid_number"].str.split(",").str[0]
    renaming_file["Grid_number"] = renaming_file["Grid_number"].astype(str)
    return renaming_file[~renaming_file["Grid_number"].str.contains("^[a-zA-Z]")]


def apply_regex_pattern(df, pattern):
    """
    Apply the regex pattern to the Sample ID column
    """
    return df.apply(
        lambda x: re.match(pattern, x).group(1) if re.match(pattern, x) else x
    )


def read_final_export_file(final_export_file):
    """
    Read the final export file into a DataFrame
    """
    return pd.read_csv(final_export_file, sep=",")


def merge_dataframes(final_export, renaming_file_filtered):
    """
    Left join final_export with samples using "Sample ID" as the key
    """
    return pd.merge(final_export, renaming_file_filtered, how="left")


def reorder_columns(final_export_grid):
    """
    Reorder the columns in the final_export_grid DataFrame
    """
    col_order = ["Grid_number"] + [
        col for col in final_export_grid.columns if col != "Grid_number"
    ]
    return final_export_grid[col_order]


def convert_grid_number_to_string(final_export_grid):
    """
    Convert Grid_number column to string type
    """
    final_export_grid["Grid_number"] = final_export_grid["Grid_number"].astype(str)
    return final_export_grid


def rename_columns(final_export_grid):
    """
    Rename 'Sample ID' to 'SequencingAcc#' and 'Grid_number' to 'Sample ID'
    """
    final_export_grid = final_export_grid.rename(
        columns={"Sample ID": "SequencingAcc#"}
    )
    final_export_grid = final_export_grid.rename(columns={"Grid_number": "Sample ID"})
    return final_export_grid


def create_copy_without_sequencing_acc(final_export_grid):
    """
    Create a copy with both sequencing Acc# and Patient ID #
    """
    final_export_grid_no_accession = final_export_grid.copy()
    final_export_grid_no_accession.drop("SequencingAcc#", axis=1, inplace=True)
    return final_export_grid_no_accession


# def write_to_file(final_export_grid, final_export_grid_no_accession):
#     """
#     Write to file
#     """
#     final_export_grid.to_csv(
#         "MatchPointExport_with_sequencingAcc_"
#         + datetime.now().strftime("%Y_%m_%d")
#         + ".txt",
#         index=False,
#         encoding="utf-8",
#         lineterminator="\n",
#     )


#     final_export_grid_no_accession.to_csv(
#         "MatchPointExport_" + datetime.now().strftime("%Y_%m_%d") + ".txt",
#         index=False,
#         encoding="utf-8",
#         lineterminator="\n",
#     )
def write_to_file(final_export_grid, final_export_grid_no_accession, directory=""):
    """
    Write data-frames to files with current date as suffix.
    """
    date_suffix = datetime.now().strftime(DATE_FORMAT)
    path_with_accession = (
        f"{directory}MatchPointExport_with_sequencingAcc_{date_suffix}.txt"
    )
    path_without_accession = f"{directory}MatchPointExport_{date_suffix}.txt"

    try:
        final_export_grid.to_csv(
            path_with_accession, index=False, encoding="utf-8", lineterminator="\n"
        )
        final_export_grid_no_accession.to_csv(
            path_without_accession, index=False, encoding="utf-8", lineterminator="\n"
        )
    except Exception as e:
        print(f"Error writing files: {e}")


def main(final_export_file, deobfuscation):
    """
    This function reads in a file exported from Soft and the ABO pipeline final_export file and
    replaces SequencingAcc# with Patient# to allow export into MatchPoint
    """
    pattern = r"(.+)_barcode\d+$"

    renaming_file = read_renaming_file(deobfuscation)
    renaming_file_filtered = preprocess_renaming_file(renaming_file)
    renaming_file_filtered.loc[:, "Sample ID"] = apply_regex_pattern(
        renaming_file_filtered["Sample ID"], pattern
    )

    final_export = read_final_export_file(final_export_file)
    final_export.loc[:, "Sample ID"] = apply_regex_pattern(
        final_export["Sample ID"], pattern
    )

    final_export_grid = merge_dataframes(final_export, renaming_file_filtered)
    final_export_grid = reorder_columns(final_export_grid)
    final_export_grid = convert_grid_number_to_string(final_export_grid)
    final_export_grid = rename_columns(final_export_grid)

    final_export_grid_no_accession = create_copy_without_sequencing_acc(
        final_export_grid
    )

    write_to_file(final_export_grid, final_export_grid_no_accession)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python rename_samples.py <final_export_file> <deobfuscation>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
