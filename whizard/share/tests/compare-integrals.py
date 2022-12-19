import re, sys
from math import sqrt

re_num = re.compile("\w+\(\w+\) =  (\d*\.\d+E.*)")
re_nums = re.compile("\w+ *(\d+\.\d+E*\S*) *(\d*\.\d+E*\S*)")
number = lambda line: float(re_num.search(line).group(1))
numbers = lambda line: (float(re_nums.search(line).group(1)),
                        float(re_nums.search(line).group(2)))
process = sys.argv[1]
filename = process + '.log'
reference_file = sys.argv[2]

with open(filename, 'r') as infile:
  for line in infile:
    if 'integral(' in line:
      integral = number(line)
    if 'error(' in line:
      error = number(line)

with open(reference_file, 'r') as infile:
  for line in infile:
    if process + ' ' in line:
      ref_integral, ref_error = numbers(line)
print 'Reference:', ref_integral, '+-', ref_error

try:
  print process, integral, '+-', error
  error_sum = sqrt (error**2 + ref_error**2)
  pull = abs (integral - ref_integral) / error_sum
  print 'pull:', pull
  valid = pull < 3
except NameError:
  print 'Number not found. Run failed.'
  valid = False

returncode = 0 if valid else 1
sys.exit(returncode)
