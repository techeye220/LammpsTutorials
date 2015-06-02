
import fileinput

import numpy as np

def parse_stdin() :
  data = []
  for line in fileinput.input() :
    if line[0] == "#" : continue
    data.append(line.strip().split()[1])

  return np.array(data,dtype=float)

if __name__ == '__main__':

  data = parse_stdin()
  print "%.3f %.3f"%(data.mean(),data.std())
