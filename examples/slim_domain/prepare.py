import os
import pandas as pd
from calvados.cfg import Config, Job, Components
import subprocess
import numpy as np
from argparse import ArgumentParser

cwd = os.getcwd()
sysname = 'test'

# set the side length of the cubic box
L = 70

# set the saving interval (number of integration steps)
N_save = 10000

# PATH
path = f'{cwd}/{sysname:s}'
subprocess.run(f'mkdir -p {path}',shell=True)
subprocess.run(f'mkdir -p data',shell=True)
residues_file = f'{cwd}/input/residues_CALVADOS3.csv'

config = Config(
  # GENERAL
  sysname = sysname, # name of simulation system
  box = [L, L, L], # nm
  temp = 298, # K
  ionic = 0.15, # molar
  pH = 7.0,
  topol = 'random',

  # RUNTIME SETTINGS
  wfreq = N_save, # dcd writing interval, 1 = 10 fs
  steps = 3.1e8, # number of simulation steps
  platform = 'CPU', # or CUDA
  threads = 1,
  restart = 'checkpoint',
  frestart = 'restart.chk',
  verbose = True,
)

config.write(path,name='config.yaml')


components = Components(
  # Defaults
  molecule_type = 'complex',
  nmol = 1, # number of molecules
  restraint = True, # apply restraints
  charge_termini = 'both', # charge N or C or both

  # INPUT
  fresidues = residues_file, # residue definitions
  ffasta = f'{cwd}/input/idr.fasta', # residue definitions
  fdomains = f'{cwd}/input/domains.yaml', # domain definitions (harmonic restraints)
  pdb_folder = f'{cwd}/input', # directory for pdb and PAE files

  # RESTRAINTS
  restraint_type = 'harmonic', # harmonic or go
  use_com = True, # apply on centers of mass instead of CA
  k_harmonic = 700., # Restraint force constant

  # NONBONDED GO CONTACTS BETWEEN MOLECULES
  #if eps_go = 0 no go potential will be added.
  eps_go = 3.6667, # default well depth (kJ/mol) for molecule pairs not in nb_go_forces.yaml

)
# one copy of every protein chain in input/test.cif (chains are read as
# components test_A, test_B, ... with domains keyed the same way in domains.yaml),
# plus 10 free copies of chain B , extra_chains={'B': 10}
components.add(name='test', nmol = 0, extra_chains={'A': 2,'B':1})
components.write(path,name='components.yaml')
