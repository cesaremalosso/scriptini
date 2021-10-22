#!/usr/bin/python

# not supported: CON, HRS, NCG, POL, SPR, STR, THE, WFC

###############################################
import numpy as np
from os import stat as fstat
import argparse
import re

def check_extra_steps(steps):
  # returns the indexes of the steps to be removed
  skip = []
  for t in range(len(steps)-1):
    if (steps[t] >= steps[t+1]):
      first = np.argmax(steps[:] == steps[t+1])
      print("  RESTART found at #{:d} step={:d} --> #{:d} step={:d} -- First occurance: #{:d} step={:d}  /-\-/-\-/-\  {:d} steps will be removed".format(t, int(steps[t]), t+1, int(steps[t+1]), first, int(steps[first]), t+1-first))
      skip.extend(range(first,t+1))
  return skip

def read_scalar_timeseries(filename):
  """Read a file of scalar/vector type:
  # step time Temp ...
  10  1.03  300.4 ...
  20  1.04  305.2 ...
  ..."""
  with open(filename, 'r') as f:
    line = f.readline()
    if (line.split()[0] == '#'):
       headline = line
    else:
       headline = []
  data = np.loadtxt(filename, dtype=str)
  steps = np.array(list(map(int, data[:,0])))
  return data, steps, headline

def read_matrix_timeseries(filename, nrows):
  """Read a file of matrix/atomic type:
  step time
  10.4356  -463.1021  67.10023
  -7.6432  -98.34529  343.3345
  ...(nrows lines per step)...
  step2 time2
  10.4356  -463.1021  67.10023
  -7.6432  -98.34529  343.3345
  ...(nrows lines per step)...
"""
  steps = []
  times = []
  data = []
  with open(filename, 'r') as f:
    while True:
      line = f.readline().split()
      if (len(line) == 0):  # EOF
        print("  END OF FILE")
        break
      if (len(line) == 2):    # new step line
        steps.append(int(line[0]))
        times.append(line[1])
      else:
        raise RuntimeError('ERROR. Wrong number of atoms?')
      data_t = ['']*nrows
      for i in range(nrows):
        line = f.readline().split()
        data_t[i] = line
      data.append(np.array(data_t))
  return np.array(data), np.array(steps), np.array(times)

def read_matrix_key_timeseries(filename, step_key='STEP:', ncomment_lines=1):
  """Read a file of matrix/atomic type with a step-KEY and comment lines.
  KEY  10  0.10
  comment comment comment ...
  -20.39  -20.11  -19.97  -19.89  -19.87
  -19.87  -19.82  -19.78  -19.72  -19.68
  ...
  KEY  20  0.20
  comment comment comment ...
  -20.39  -20.11  -19.97  -19.89  -19.87
  -19.87  -19.82  -19.78  -19.72  -19.68
  ...
"""
  with open(filename, 'r') as f:
    i = -1
    while True:
      line = f.readline().split()
      if (len(line) == 0):  # EOF
        raise RuntimeError('ERROR. step_key not found.')
      if (line[0] == step_key):
        if (i == -1):
          i = 0
        else:
          break
      elif (i >= 0):
        i += 1
  nrows = i - ncomment_lines
  print("  nrows = ", nrows)
  
  steps = []
  times = []
  comment = []
  data = []
  with open(filename, 'r') as f:
    while True:
      line = f.readline().split()
      if (len(line) == 0):  # EOF
        print("  END OF FILE")
        break
      if (line[0] == step_key):   # new step line
        steps.append(int(line[1]))
        times.append(line[2])
      else:
        raise RuntimeError('ERROR. Wrong number of comment lines?\n STEP={}\n line={}'.format(steps[-1],line))
      comm_t = ['']*ncomment_lines
      for i in range(ncomment_lines):
        line = f.readline()
        comm_t[i] = line
      comment.append(np.array(comm_t))
      data_t = ['']*nrows
      for i in range(nrows):
        line = f.readline()
        data_t[i] = line
      data.append(data_t)
  return data, np.array(steps), np.array(times), np.array(comment)

#################################################

def main ():
  """This program reads the outputs of a CP-QE simulation and fixes the timeseries from restarts that may have happened.
  Example:
    python fix_cp_traj.py silica silicaok

  Supported files:     .cel .eig .evp .for .nos .pos .str .vel
  Not (yet) supported: .con .hrs .ncg .pol .spr .the .wfc
  """
  _epilog = """---
  Code written by Loris Ercole.
  SISSA, Via Bonomea, 265 - 34136 Trieste ITALY
  """
  parser = argparse.ArgumentParser()
  parser = argparse.ArgumentParser(description=main.__doc__, epilog=_epilog, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('prefix', type=str, help='prefix of the files in the QE output directory.')
  parser.add_argument('outprefix', type=str, help='new prefix of the files to be written.')
  parser.add_argument('-n', '--natoms', type=int, required=True, help='number of atoms')
  args = parser.parse_args()
  prefix = args.prefix
  outprefix = args.outprefix
  if (prefix == outprefix):
    raise ValueError('For your data safety, input and output prefix should differ!')
    return 1
  natoms = args.natoms

  term = ['.cel', '.eig', '.evp', '.for', '.nos', '.pos', '.str', '.vel'] 
  read_file = {}

  for tt in term:
    try:
      filein = open(prefix + tt, 'r')
      read_file[tt] = True
      print('The ' + tt +' file is present: it will be fixed')
    except IOError as err:
      err = re.sub(r'\[.*\]', '', str(err))
      print('Warning!' + err + '. The ' + tt +' file is NOT present: it will be ignored')
      read_file[tt] = False
  
  if read_file['.cel'] == True:
    # CEL - box cell file
    filename = prefix + '.cel'
    print("* Fixing {} file...".format(filename))
    outfilename = outprefix + '.cel'
    if fstat(filename).st_size:  # if file is not empty
      try:
       data, steps, times = read_matrix_timeseries(filename, 3)
       skip = check_extra_steps(steps)
       data = np.delete(data, skip, 0)
       steps = np.delete(steps, skip, 0)
       times = np.delete(times, skip, 0)
       with open(outfilename, 'w') as f:
         for s, t, d in zip(steps, times, data):
           f.write(str(s) + ' ' + t + '\n')
           np.savetxt(f, d, fmt=' %s')
       print("  --> {} written".format(outfilename))
      except RuntimeError as e:
        print('Error reading file.\n{} {}'.format(e.errno, e.strerror))
    else:
      print(" {} is empty.".format(filename))

  if read_file['.evp'] == True:
    # EVP - thermo file
    filename = prefix + '.evp'
    print("* Fixing {} file...".format(filename))
    outfilename = outprefix + '.evp'
    if fstat(filename).st_size:  # if file is not empty
      try:
        data, steps, headline = read_scalar_timeseries(filename)
        skip = check_extra_steps(steps)
        data = np.delete(data, skip, 0)
        with open(outfilename, 'w') as f:
          if len(headline):
            f.write(headline)
          np.savetxt(f, data, fmt=' %s')
        print("  --> {} written".format(outfilename))
      except RuntimeError as e:
        print('Error reading file.\n{} {}'.format(e.errno, e.strerror))
    else:
      print(" {} is empty.".format(filename))
  
  if read_file['.for'] == True:
    # FOR - forces file
    filename = prefix + '.for'
    print("* Fixing {} file...".format(filename))
    outfilename = outprefix + '.for'
    if fstat(filename).st_size:  # if file is not empty
      try:
        data, steps, times = read_matrix_timeseries(filename, natoms)
        skip = check_extra_steps(steps)
        data = np.delete(data, skip, 0)
        steps = np.delete(steps, skip, 0)
        times = np.delete(times, skip, 0)
        with open(outfilename, 'w') as f:
          for s, t, d in zip(steps, times, data):
            f.write(str(s) + ' ' + t + '\n')
            np.savetxt(f, d, fmt=' %s')
        print("  --> {} written".format(outfilename))
      except RuntimeError as e:
        print('Error reading file.\n{} {}'.format(e.errno, e.strerror))
    else:
      print(" {} is empty.".format(filename))
    
  if read_file['.nos'] == True:
    # NOS - nos file
    filename = prefix + '.nos'
    print("* Fixing {} file...".format(filename))
    outfilename = outprefix + '.nos'
    if fstat(filename).st_size:  # if file is not empty
      try:
        data, steps, headline = read_scalar_timeseries(filename)
        skip = check_extra_steps(steps)
        data = np.delete(data, skip, 0)
        with open(outfilename, 'w') as f:
          if len(headline):
            f.write(headline)
          np.savetxt(f, data, fmt=' %s')
        print("  --> {} written".format(outfilename))
      except RuntimeError as e:
        print('Error reading file.\n{} {}'.format(e.errno, e.strerror))
    else:
      print(" {} is empty.".format(filename))
    
  if read_file['.pos'] == True:
    # POS - positions file
    filename = prefix + '.pos'
    print("* Fixing {} file...".format(filename))
    outfilename = outprefix + '.pos'
    if fstat(filename).st_size:  # if file is not empty
      try:
        data, steps, times = read_matrix_timeseries(filename, natoms)
        skip = check_extra_steps(steps)
        data = np.delete(data, skip, 0)
        steps = np.delete(steps, skip, 0)
        times = np.delete(times, skip, 0)
        with open(outfilename, 'w') as f:
          for s, t, d in zip(steps, times, data):
            f.write(str(s) + ' ' + t + '\n')
            np.savetxt(f, d, fmt=' %s')
        print("  --> {} written".format(outfilename))
      except RuntimeError as e:
        print('Error reading file.\n{} {}'.format(e.errno, e.strerror))
    else:
      print(" {} is empty.".format(filename))
    
  if read_file['.str'] == True:
    # STR - box stress tensor file
    filename = prefix + '.str'
    print("* Fixing {} file...".format(filename))
    outfilename = outprefix + '.str'
    if fstat(filename).st_size:  # if file is not empty
      try:
        data, steps, times = read_matrix_timeseries(filename, 3)
        skip = check_extra_steps(steps)
        data = np.delete(data, skip, 0)
        steps = np.delete(steps, skip, 0)
        times = np.delete(times, skip, 0)
        with open(outfilename, 'w') as f:
          for s, t, d in zip(steps, times, data):
            f.write(str(s) + ' ' + t + '\n')
            np.savetxt(f, d, fmt=' %s')
        print("  --> {} written".format(outfilename))
      except RuntimeError as e:
        print('Error reading file.\n{} {}'.format(e.errno, e.strerror))
    else:
      print("  ! {} is empty.".format(filename))
  
  if read_file['.vel'] == True:
    # VEL - velocities file
    filename = prefix + '.vel'
    print("* Fixing {} file...".format(filename))
    outfilename = outprefix + '.vel'
    if fstat(filename).st_size:  # if file is not empty
      try:
        data, steps, times = read_matrix_timeseries(filename, natoms)
        skip = check_extra_steps(steps)
        data = np.delete(data, skip, 0)
        steps = np.delete(steps, skip, 0)
        times = np.delete(times, skip, 0)
        with open(outfilename, 'w') as f:
          for s, t, d in zip(steps, times, data):
            f.write(str(s) + ' ' + t + '\n')
            np.savetxt(f, d, fmt=' %s')
        print("  --> {} written".format(outfilename))
      except RuntimeError as e:
        print('Error reading file.\n{} {}'.format(e.errno, e.strerror))
    else:
      print(" {} is empty.".format(filename))

  if read_file['.eig'] == True:
    # EIG - eigenvalues file
    filename = prefix + '.eig'
    print("* Fixing {} file...".format(filename))
    outfilename = outprefix + '.eig'
    if fstat(filename).st_size:  # if file is not empty
      try:
        data, steps, times, comment = read_matrix_key_timeseries(filename, 'STEP:', 1)
        skip = check_extra_steps(steps)
        data = np.delete(data, skip, 0)
        steps = np.delete(steps, skip, 0)
        times = np.delete(times, skip, 0)
        comment = np.delete(comment, skip, 0)
        with open(outfilename, 'w') as f:
          for s, t, c, d in zip(steps, times, comment, data):
            f.write('  STEP:  ' + str(s) + ' ' + t + '\n')
            f.write(c)
            for dd in d:
              f.write(dd)
        print("  --> {} written".format(outfilename))
      except RuntimeError as e:
        print('Error reading file.\n{} {}'.format(e.errno, e.strerror))
    else:
      print(" {} is empty.".format(filename))
  
  return 0

if __name__ == "__main__":
  main()
