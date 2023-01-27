#!python
#
# Simple RMSD and RMSF analysis of simulations
# Code 2022 by Peter Kasson

import glob
import json
import numpy
import os

gmxbin = 'srun gmx'

def analyze_one(basename):
  """Runs analysis for one trajectory.
  Args: basename: basename for trajectory
  Rets: RMSD aligned by chain
  """
  print(basename)
  rmsd = []
  for chain_idx in range(3):  
    tmpdat = numpy.loadtxt('%s_rmsd_%d.xvg' % (basename, chain_idx),
                           comments=['#', '@'])
    rmsd.append(tmpdat[:, 1]**2)
  # re-calculate RMSD
  return numpy.sqrt(numpy.mean(numpy.vstack(rmsd), axis=0))

def analyze_jointalign(basename):
  """Runs analysis for one trajectory.
  Args: basename: basename for trajectory
  Rets: RMSD aligned as whole unit
  """
  tmpdat = numpy.loadtxt('%s_rmsd.xvg' % basename,
                         comments=['#', '@'])
  return tmpdat[:, 1]

if __name__ == '__main__':
  rmsd_dict = {}
  for xtcfile in glob.glob('*prod.xtc'):
    keyname = xtcfile[:-4]
    rmsd_dict[keyname] = analyze_one(keyname).tolist()
  outfile = open('rmsd_vals.json', 'w')
  json.dump(rmsd_dict, outfile)
  outfile.close()

  rmsd_dict = {}
  for xtcfile in glob.glob('*prod.xtc'):
    keyname = xtcfile[:-4]
    rmsd_dict[keyname] = analyze_jointalign(keyname).tolist()
  outfile = open('rmsd_jointvals.json', 'w')
  json.dump(rmsd_dict, outfile)
  outfile.close()
