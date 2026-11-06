from data_processing.data_subject import DataSubject

def count_cohort_split(cohort):
    """Helper function to calculate the data split between IDH-wildtype and IDH-mutant cases in a cohort."""
    wt = 0
    m = 0

    for subj in cohort:
        if subj.label == 0:
            wt += 1
        else:
            m += 1

    return wt, m