#!/home/andrew/anaconda3/envs/research/bin/python
from typing import Final
import os
import shutil
from glob import glob
import subprocess
from subprocess import CompletedProcess
import pandas as pd
from pandas import DataFrame

#===== Constants =====#
TOLERANCE: Final[float] = 1e-6
VALIDATION_DIR: Final[str] = "BaselineResults"
RESULTS_FILENAME: Final[str] = "KonzaSurface_results.txt"
NUMEXP_FILENAME: Final[str] = "KonzaSurface_results_numexp.txt"
RESULTS_PATH: Final[str] = os.path.join(VALIDATION_DIR, RESULTS_FILENAME)
NUMEXP_PATH: Final[str] = os.path.join(VALIDATION_DIR, NUMEXP_FILENAME)
BIORT_BUILD_PATH: Final[str] = "build/src/biort"
#=====================#

#===== Functions =====#
def get_result_filepaths() -> list[str]:
    all_files: list[str] = glob("output/*.txt")
    files = list(
        filter(
            lambda x: ("results" in x) and ("numexp" not in x),
            all_files
        )
    )
    # print(f"Results output files: {files}")
    return files

def get_numexp_filepaths() -> list[str]:
    all_files: list[str] = glob("output/*.txt")
    files = list(
        filter(
            lambda x: ("results" in x) and ("numexp" in x),
            all_files
        )
    )
    # print(f"Numexp output files: {files}")
    return files

def get_newest_file(lst: list[str]) -> str:
    return max(lst, key=os.path.getctime)

def get_err(df1: DataFrame, df2: DataFrame) -> float:
    err: float = 0.0

    col: str
    for col in df1.columns:
        err += df1[col].subtract(df2[col]).abs().sum()

    return err
#=====================#

#===== Main function =====#
def main() -> None:
    # Copy the build executable
    shutil.copy(BIORT_BUILD_PATH, ".")

    # Run the model
    subprocess.run(
        "./biort KonzaSurface",
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )
    
    # Get the filepaths
    newest_result_path: Final[str] = get_newest_file(get_result_filepaths())
    newest_numexp_path: Final[str] = get_newest_file(get_numexp_filepaths())

    # print(f"Results validation path: {RESULTS_PATH}")
    # print(f"Results test path: {newest_result_path}")
    # print(f"Numexp validation path: {NUMEXP_PATH}")
    # print(f"Numexp test path: {newest_numexp_path}")

    # Load the dataframes
    res_val_df: DataFrame = pd.read_csv(RESULTS_PATH, skiprows=[1], delim_whitespace=True)
    res_tst_df: DataFrame = pd.read_csv(newest_result_path, skiprows=[1], delim_whitespace=True)

    num_val_df: DataFrame = pd.read_csv(NUMEXP_PATH, skiprows=[1], delim_whitespace=True)
    num_tst_df: DataFrame = pd.read_csv(newest_numexp_path, skiprows=[1], delim_whitespace=True)

    diff_res: float = get_err(res_val_df, res_tst_df)
    diff_num: float = get_err(num_val_df, num_tst_df)

    if diff_res > TOLERANCE or diff_num > TOLERANCE:
        print(f"Error in results is too great: {diff_res}")
        diff_res_df: DataFrame = res_val_df.subtract(res_tst_df)
        diff_num_df: DataFrame = num_val_df.subtract(num_tst_df)
        res_ax_err = res_val_df.subtract(res_tst_df).abs().sum(axis=0)
        num_ax_err = (num_val_df - num_tst_df).abs().sum(axis=0)
        print("Results error:")
        print(res_ax_err)
        print("Numexp error: ")
        print(num_ax_err)
        diff_res_df.to_csv("diff_res.csv")
        diff_num_df.to_csv("diff_num.csv")
    else:
        print("Build was successful")
        # Remove outputs
        outputs: Final[list[str]] = glob("output/*.txt")
        for f in outputs:
            os.remove(f)


#=========================#

if __name__ == "__main__":
    main()