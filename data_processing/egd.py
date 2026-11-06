import pandas as pd
import numpy as np
from pathlib import Path
import re
import logging

from data_processing.data_subject import DataSubject
from data_processing.utils import count_cohort_split

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%H:%M')
SEQUENCES = {
    "T1" : 1,
    "T1GD": 2,
    "T2": 3,
    "FLAIR": 4,
    "MASK": 5
}
ROOT = Path("../data/EGD/egd_data/images/")
METADATA_PATH = "../data/EGD/egd_data/metadata/subjects.csv"
SEX_MAP = {"male":"M", "female":"F"}


def find_sequence(patient_id, suffix):
    p = ROOT / patient_id / f"MR_{patient_id}" / f"scan{SEQUENCES[suffix]}_{suffix}_{suffix}" / "NIFTI" / f"{suffix}.nii.gz"
    return p if p.exists() else None

def idh_status_map(idh):
    if idh == 'wildtype':
        return 0
    else:
        return 1


def load_egd_data():
    """
    Data pipeline to extract the relevant patient data for the EGD patient cohort.
    Extracts relevant clinical data and MRI file paths, and normalises patient sex to M/F.
    Returns a list of DataSubject objects, one for each patient.
    """

    df = pd.read_csv(METADATA_PATH)
    df = df[["subject_label", "fields/field[0]/field", "demographics/gender", "demographics/age"]]
    df = df.rename(columns={"fields/field[0]/field":"IDH", "demographics/gender":"sex", "demographics/age":"age"})
    df = df.replace([-1, "-1"], np.nan)
    df = df.replace(r"(?i)^\s*unknown\s*$", np.nan, regex=True)
    df = df.dropna().reset_index(drop=True)
    df = df.replace({"sex": SEX_MAP})
    egd_df = df.copy()

    for seq in SEQUENCES.keys():
        egd_df[seq] = egd_df["subject_label"].apply(lambda pid: find_sequence(pid, seq))

    complete = egd_df.notna().all(axis=1)
    if not complete.all():
        raise ValueError("UTSW Patient Data is incomplete - some row of DataFrame is missing columnn values")

    egd_subjects = []
    for index, row in egd_df.iterrows():
        egd_subjects.append(
            DataSubject(
                cohort = "EGD",
                patient_id = row["subject_label"],
                t1 = row["T1"],
                t1c = row["T1GD"],
                t2 = row["T2"],
                flair = row["FLAIR"],
                mask = row["MASK"],
                age = row["age"],
                sex = row["sex"],
                label = row["IDH"]
            )
        )

    wt, m = count_cohort_split(egd_subjects)
    logging.info(f"EGD Cohort - {len(egd_subjects)} subjects loaded ({wt} wildtype vs {m} mutant)")
    return egd_subjects