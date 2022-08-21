#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# v.0.1
# PaDEL-robust: robust version of PaDEL-descriptor Python interface
# Developed in 2022 by Mikhail Markovsky <m.markovsky@gmail.com>

import os
from subprocess import Popen, DEVNULL, TimeoutExpired
from shutil import which
import uuid
import psutil
from csv import DictReader
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

PROCESS_DESCRIPTION = 'Calculating descriptors'
CLEANING_DESCRIPTION = 'Cleaning temporary data'

# Path to the main Java app
_PADEL_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'PaDEL-Descriptor',
    'PaDEL-Descriptor.jar'
)


def _kill(proc_pid):
    """
    Kills the defined process with all children if exist
    :param proc_pid: PID of the target process
    :return: None
    """
    try:
        process = psutil.Process(proc_pid)
        for proc in process.children(recursive=True):
            proc.kill()
        process.kill()
    except psutil.NoSuchProcess:
        pass


def _test_java():
    """
    Tests if Java is installed
    :return:
    """
    if which('java') is None:
        raise ReferenceError('Java JRE 6+ not found (required for PaDEL-Descriptor)')


class PadelDescriptor:
    def __init__(self,
                 headless: bool = True,
                 max_runtime: int = -1,
                 waiting_jobs: int = -1,
                 threads: int = -1,
                 d_2d: bool = False,
                 d_3d: bool = False,
                 config: str = None,
                 convert_3d: bool = False,
                 descriptor_types: str = None,
                 detect_aromaticity: bool = False,
                 fingerprints: bool = False,
                 log: bool = False,
                 max_cpd_per_file: int = 0,
                 remove_salt: bool = False,
                 retain_3d: bool = False,
                 retain_order: bool = False,
                 standardize_nitro: bool = False,
                 standardize_tautomers: bool = False,
                 tautomer_list: str = None,
                 use_filename_as_molname: bool = False,
                 timeout: int = None,
                 pool_size: int = None,
                 use_tqdm: bool = True,
                 temp_dir: str = 'padel_temp'):
        """
        Initializes the PaDEL descriptor instance
        :param headless: Use Java app without notifications?
        :param max_runtime: Maximum running time per molecule (in mS); defaults to -1 (unlimited) (not recommended)
        :param waiting_jobs: Maximum number of jobs to store in queue for worker
            threads to process; defaults to -1 (50 * max threads) (not recommended)
        :param threads: maximum number of threads to use; defaults to -1 (equal to number of CPU cores)
            (not recommended)
        :param d_2d: If `True`, calculates 2-D descriptors
        :param d_3d: If `True`, calculates 3-D descriptors
        :param config: Path to configuration file (optional)
        :param convert_3d: If `True`, converts molecule to 3-D
        :param descriptor_types: Path to descriptor types file (optional)
        :param detect_aromaticity: If `True`, removes existing aromaticity
            information and automatically detect aromaticity in the molecule
            before calculation of descriptors
        :param fingerprints: If `True`, calculates fingerprints
        :param log: if `True`, Creates a log file (same as descriptors file, with .log extension)
        :param max_cpd_per_file: Maximum number of compounds to be stored in each
            descriptor file; defaults to 0 (unlimited)
        :param remove_salt: If `True`, removes salt from the molecule
        :param retain_3d: If `True`, retains 3-D coordinates when standardizing structure
        :param retain_order: If `True`, retains order of molecules in
            structural files for descriptor file
        :param standardize_nitro: If `True`, standardizes nitro groups to N(:O):O
        :param standardize_tautomers: If `True`, standardizes tautomers
        :param tautomer_list: Path to SMIRKS tautomers file (optional)
        :param use_filename_as_molname: If `True`, uses filename (minus the
            extension) as the molecule name
        :param timeout: Timeout for each thread in seconds
        :param pool_size: Size of multiprocessing pool (number of CPU cores by default)
        :param use_tqdm: Use TQDM progress bars for a time estimation
        :param temp_dir: Temporary directory for SMILES and CSV files
        """
        self._headless = headless
        self._max_runtime = max_runtime
        self._waiting_jobs = waiting_jobs
        self._threads = threads
        self._d_2d = d_2d
        self._d_3d = d_3d
        self._config = config
        self._convert_3d = convert_3d
        self._descriptor_types = descriptor_types
        self._detect_aromaticity = detect_aromaticity
        self._fingerprints = fingerprints
        self._log = log
        self._max_cpd_per_file = max_cpd_per_file
        self._remove_salt = remove_salt
        self._retain_3d = retain_3d
        self._retain_order = retain_order
        self._standardize_nitro = standardize_nitro
        self._standardize_tautomers = standardize_tautomers
        self._tautomer_list = tautomer_list
        self._use_filename_as_molname = use_filename_as_molname
        self._timeout = timeout
        if pool_size is None:
            self._pool_size = cpu_count()
        else:
            self._pool_size = pool_size
        self._use_tqdm = use_tqdm
        self._temp_dir = temp_dir

        _test_java()

        if not os.path.exists(self._temp_dir):
            os.mkdir(self._temp_dir)

    def make_descriptor(self, smiles):
        """
        Makes a prediction for a single molecule
        :param smiles: SMILES string
        :return: Dictionary with the descriptors
        """
        if self._headless:
            command = 'java -Xms1G -Xmx1G -Djava.awt.headless=true -jar {}'.format(_PADEL_PATH)
        else:
            command = 'java -jar {}'.format(_PADEL_PATH)
        command += ' -maxruntime {}'.format(self._max_runtime)
        command += ' -waitingjobs {}'.format(self._waiting_jobs)
        command += ' -threads {}'.format(self._threads)
        command += ' -maxcpdperfile {}'.format(self._max_cpd_per_file)
        if self._d_2d is True:
            command += ' -2d'
        if self._d_3d is True:
            command += ' -3d'
        if self._config is not None:
            command += ' -config {}'.format(self._config)
        if self._convert_3d is True:
            command += ' -convert3d'
        if self._descriptor_types is not None:
            command += ' -descriptortypes {}'.format(self._descriptor_types)
        if self._detect_aromaticity is True:
            command += ' -detectaromaticity'
        if self._fingerprints is True:
            command += ' -fingerprints'
        if self._log is True:
            command += ' -log'
        if self._remove_salt is True:
            command += ' -removesalt'
        if self._retain_3d is True:
            command += ' -retain3d'
        if self._retain_order is True:
            command += ' -retainorder'
        if self._standardize_nitro is True:
            command += ' -standardizenitro'
        if self._standardize_tautomers is True:
            command += ' -standardizetautomers'
        if self._tautomer_list is not None:
            command += ' -tautomerlist {}'.format(self._tautomer_list)
        if self._use_filename_as_molname is True:
            command += ' -usefilenameasmolname'

        base_name = os.path.join(self._temp_dir, str(uuid.uuid4()))
        smi_name = base_name + '.smi'
        csv_name = base_name + '.csv'
        command += ' -dir {}'.format(smi_name)
        command += ' -file {}'.format(csv_name)
        process = None
        try:
            with open(smi_name, 'w') as smi_file:
                smi_file.write(smiles)

            process = Popen(command, stdout=DEVNULL, shell=True)
            process.wait(self._timeout)

        except TimeoutExpired:
            pass

        finally:
            _kill(process.pid)

            with open(csv_name, 'r', encoding='utf-8') as desc_file:
                reader = DictReader(desc_file)
                rows = [row for row in reader]

            os.remove(smi_name)
            os.remove(csv_name)

        if len(rows) > 0:
            rows[0]['Name'] = smiles
            rows = rows[0]
        else:
            rows = {}

        return rows

    def make_descriptors_batch(self, smiles_list):
        """
        Prepares a batch of predictors
        :param smiles_list: List of SMILES
        :return: List of dicts with descriptors
        """
        tasks = len(smiles_list)
        with Pool(self._pool_size) as p:
            if self._use_tqdm:
                desc_array = list(
                    tqdm(
                        p.imap(
                            self.make_descriptor,
                            smiles_list
                        ),
                        total=tasks,
                        desc=PROCESS_DESCRIPTION
                    )
                )
            else:
                desc_array = list(p.imap(self.make_descriptor, smiles_list))
        if self._use_tqdm:
            desc_array = [
                i for i in tqdm(
                    desc_array,
                    desc=CLEANING_DESCRIPTION
                ) if isinstance(i, dict)
            ]
        else:
            desc_array = [i for i in desc_array if isinstance(i, dict)]
        return desc_array
