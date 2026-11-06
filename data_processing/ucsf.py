import pandas as pd
from pathlib import Path
import re
import logging

from data_processing.data_subject import DataSubject
from data_processing.utils import count_cohort_split

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%H:%M')
SEQUENCES = [
    "T1",
    "T1c",
    "T2",
    "FLAIR"
]
ROOT = Path("../../data/UCSF-PDGM/PKG - UCSF-PDGM Version 5/UCSF-PDGM-v5/")
METADATA_PATH = "../../data/UCSF-PDGM/UCSF-PDGM-metadata_v5.csv"


def folder_id(csv_id, width=4):
    num = int(re.search(r"(\d+)$", csv_id).group(1))
    return f"UCSF-PDGM-{num:0{width}d}"

def find_sequence(csv_id, suffix):
    fid = folder_id(csv_id)
    p = ROOT / f"{fid}_nifti" / f"{fid}_{suffix}_bias.nii.gz"
    return p if p.exists() else None

def find_mask(csv_id):
    fid = folder_id(csv_id)
    p = ROOT / f"{fid}_nifti" / f"{fid}_tumor_segmentation.nii.gz"
    return p if p.exists() else None

def idh_status_map(idh):
    if idh == 'wildtype':
        return 0
    else:
        return 1


def load_ucsf_data():
    """
    Data pipeline to extract the relevant patient data for the UCSF-PDGM patient cohort.
    Extracts relevant clinical data and MRI file paths, and normalises IDH mutation status to 0/1.
    Returns a list of DataSubject objects, one for each patient.
    """

    df = pd.read_csv(METADATA_PATH)
    df = df[["ID", "Sex", "Age at MRI", "IDH"]]

    # Remove follow up imaging of pre-existing patients
    # More information in 'Detailed Description' @ https://www.cancerimagingarchive.net/collection/ucsf-pdgm/
    ucsf_df = df[~df['ID'].str.contains(r'FU\d{3}d', regex=True, na=False)] 

    ucsf_df = ucsf_df.copy()
    ucsf_df["IDH"] = ucsf_df["IDH"].apply(lambda idh: idh_status_map(idh))
    ucsf_df["patient_id"] = ucsf_df["ID"].apply(folder_id)

    for seq in SEQUENCES:
        ucsf_df[seq] = ucsf_df["ID"].apply(lambda pid: find_sequence(pid, seq))

    ucsf_df["mask"] = ucsf_df["ID"].apply(lambda pid: find_mask(pid))
    complete = ucsf_df.notna().all(axis=1)
    if not complete.all():
        raise ValueError("UCSF Patient Data is incomplete - some row of DataFrame is missing columnn values")

    ucsf_subjects = []
    for _, row in ucsf_df.iterrows():
        ucsf_subjects.append(
            DataSubject(
                cohort = "UCSF-PDGM",
                patient_id = row["ID"],
                t1 = row["T1"],
                t1c = row["T1c"],
                t2 = row["T2"],
                flair = row["FLAIR"],
                mask = row["mask"],
                age = row["Age at MRI"],
                sex = row["Sex"],
                label = row["IDH"]
            )
        )

    wt, m = count_cohort_split(ucsf_subjects)
    logging.info(f"UCSF-PDGM Cohort - {len(ucsf_subjects)} subjects loaded ({wt} wildtype vs {m} mutant)")
    return ucsf_subjects

    