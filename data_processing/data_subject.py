import nibabel as nib
import numpy as np

class DataSubject:
    """
    Data representation of adult diffuse glioma patients for all cohorts used in this study.
    Includes:
        - Cohort name and patient id metatdata to allow for identification of patients from cohorts
        - Clinical data for age and sex (M/F) of the patient
        - Label for idh mutation status (0 represents wildtype and 1 represents mutant)
        - .nii.gz file paths for 4 MRI sequences (T1, T1c, T2 and FLAIR)
        - .nii.gz file path for the tumour segmentation mask of the patient
    """
    def __init__(self, cohort, patient_id, t1, t1c, t2, flair, mask, age, sex, label):
        self.cohort = cohort
        self.patient_id = patient_id

        self.t1 = t1
        self.t1c = t1c
        self.t2 = t2
        self.flair = flair
        self.mask = mask

        self.age = age
        self.sex = sex
        self.label = label

    def load(self, sequence):
        """Load one sequence as a numpy array, e.g. subject.load("flair")."""
        path = getattr(self, sequence)
        return np.asarray(nib.load(path).dataobj, dtype=np.float32)

    def load_nifti(self, sequence):
        "Load one sequence as a nibabel image."
        path = getattr(self, sequence)
        return nib.load(path)