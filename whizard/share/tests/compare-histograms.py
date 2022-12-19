import re, sys
from math import sqrt

delta = 1E-12
re_nums = re.compile(" +\S+ +(\d+\.\d+E\S*)  (\d*\.\d+E\S*)")
re_bin = re.compile(" +(\d*\.\d+E\S*)")
bin_midpoint = lambda line: float(re_bin.search(line).group(1))
numbers = lambda line: (float(re_nums.search(line).group(1)),
                        float(re_nums.search(line).group(2)))

def pull_and_chi2(line, ref_line):
  value, error = numbers(line)
  ref_value, ref_error = numbers(ref_line)
  error_sum = sqrt (error**2 + ref_error**2)
  diff = abs (value - ref_value)
  chi2_result = diff ** 2 / ref_value
  if abs (error_sum) > delta:
    pull_result = diff / error_sum
  else:
    pull_result = 0.0 if diff < delta else 10.0
  return pull_result, chi2_result

process = sys.argv[1]
filename = process + '_hist.dat'
reference = 'ref-output/' + process + '.ref'
pull_max = 0.0
dof = -1
sum_chi2 = 0.0

print (3 * '{:<30s}').format("Bin Midpoint", "Pull", "Chi2")
# In Python 2.7 we could write this in one line
with open(filename, 'r') as infile:
  with open(reference, 'r') as ref_file:
    for line, ref_line in zip(infile, ref_file):
      if '#' not in line and len(line) > 1:
        dof += 1
        this_pull, this_chi2 = pull_and_chi2(line, ref_line)
        sum_chi2 += this_chi2
        pull_max = max(pull_max, this_pull)
        print (3 * '{:<30.10f}').format(bin_midpoint(line), this_pull, this_chi2)
      elif 'Underflow' in line:
        break
chi2_per_dof = sum_chi2 / dof
print (4 * '{:<30s}').format("Max Pull", "Chi2", "dof", "Chi2/dof")
print (4 * '{:<30.10f}').format(pull_max, sum_chi2, dof, chi2_per_dof)
returncode = 0 if pull_max < 3 and chi2_per_dof < 3 else 1
print returncode
sys.exit(returncode)
