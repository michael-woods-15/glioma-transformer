import logging

from data_processing.ucsf import load_ucsf_data
from data_processing.utsw import load_utsw_data
from data_processing.egd import load_egd_data

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%H:%M')

def run_data_pipeline():
    ucsf_subjects = load_ucsf_data()
    utsw_subjects = load_utsw_data()
    egd_subjects = load_egd_data()


if __name__ == '__main__':
    run_data_pipeline()