import pandas as pd
from pathlib import Path
import re
import logging

from data_processing.data_subject import DataSubject
from data_processing.utils import count_cohort_split

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%H:%M')

SEQUENCES = [
    "t1",
    "t1ce",
    "t2",
    "flair"
]
ROOT = Path("../../data/UCSF-PDGM/PKG - UCSF-PDGM Version 5/UCSF-PDGM-v5/")
METADATA_PATH = "../../data/UTSW-Glioma/UTSW_Glioma_Metadata.tsv"


def find_sequence(patient_id, suffix):
    p = ROOT / patient_id / f"brain_{suffix}.nii.gz"
    return p if p.exists() else None

def find_mask(patient_id):
    p = ROOT / patient_id / f"tumorseg_FeTS.nii.gz"
    return p if p.exists() else None

def idh_status_map(idh):
    if idh == 'wild type':
        return 0
    else:
        return 1


def load_utsw_data():
    """
    Data pipeline to extract the relevant patient data for the UTSW-Glioma patient cohort.
    Extracts relevant clinical data and MRI file paths, and normalises IDH mutation status to 0/1.
    Returns a list of DataSubject objects, one for each patient.
    """

    df = pd.read_csv(METADATA_PATH, sep="\t")
    df = df[["Subject ID", "Sex at birth", "Age at Imaging", "IDH", "T1", "T1GD", "T2", "T2FLAIR"]]
    df = df.dropna(subset=["IDH"])
    df = df[["Subject ID", "Sex at birth", "Age at Imaging", "IDH"]]
    df["IDH"] = df["IDH"].apply(lambda idh: idh_status_map(idh))
    utsw_df = df.copy()

    for seq in SEQUENCES:
        utsw_df[seq] = utsw_df["Subject ID"].apply(lambda pid: find_sequence(pid, seq))

    utsw_df['mask'] = utsw_df["Subject ID"].apply(lambda pid: find_mask(pid))
    complete = utsw_df.notna().all(axis=1)
    if not complete.all():
        raise ValueError("UTSW Patient Data is incomplete - some row of DataFrame is missing columnn values")

    utsw_subjects = []

    for _, row in utsw_df.iterrows():
        utsw_subjects.append(
            DataSubject(
                cohort = "UTSW-Glioma",
                patient_id = row["Subject ID"],
                t1 = row["t1"],
                t1c = row["t1ce"],
                t2 = row["t2"],
                flair = row["flair"],
                mask = row["mask"],
                age = row["Age at Imaging"],
                sex = row["Sex at birth"],
                label = row["IDH"]
            )
        )

    wt, m = count_cohort_split(utsw_subjects)
    logging.info(f"UTSW-Glioma Cohort - {len(utsw_subjects)} subjects loaded ({wt} wildtype vs {m} mutant)")
    return utsw_subjects