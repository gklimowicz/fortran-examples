import re, sys
from math import sqrt

re_w = re.compile("E .* (-*\d*\.\d+e\S*)")
re_nums = re.compile("C (\d+\.\d+e\S*) (\d*\.\d+e\S*)")
def xsection(line):
  try:
    return float(re_nums.search(line).group(1)), float(re_nums.search(line).group(2))
  except AttributeError:
    print 'Could not find xsection in this line:'
    print line
    return 0.0

def weight(line):
  try:
    weight = float(re_w.search(line).group(1))
    return weight
  except AttributeError:
    print 'Could not find number in this line:'
    print line
    return 0.0

test = sys.argv[1]
filename = test + '.hepmc'
valid = True
xsec = 0.0
error = 0.0
sum_weights = [0.0, 0.0, 0.0]
NN = [0, 0, 0]
line_no = 0
region = 0
regions = range(3)

print 'Start file ' + filename
with open(filename, 'r') as infile:
  for line in infile:
    if 'C ' in line:
      xsec, error = xsection(line)
      break
  print 'Cross section is', xsec, '+-', error
with open(filename, 'r') as infile:
  for line in infile:
    line_no += 1
    if 'E ' in line:
      this_weight = weight(line)
      #print 'line', line_no, this_weight
      sum_weights[region] += this_weight
      NN[region] += 1
      region += 1
      if region == 3:
        region = 0
  print 'HepMC'
  for region in regions:
    mean = sum_weights[region] / NN[region]
    print 'region', region
    print 'NN:', NN[region]
    print 'sum_weights:', sum_weights[region]
    print 'mean:', mean
    print 50 * '-'
  print 'Overall'
  mean = sum(sum_weights) / sum(NN)
  print 'mean:', mean
  print 'abs(mean - xsec):', abs(mean-xsec)
  print 'error:', error
  pull = abs(mean-xsec) / error
  print 'pull:', pull
  valid = pull < 3
  print 'valid:', valid
  print 'End file ' + filename

returncode = 0 if valid else 1
print 'returncode:', returncode
sys.exit(returncode)
