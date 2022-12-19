import re, sys
from math import sqrt

process = sys.argv[1]
filename = process + '.log'

integrals = []
methods = []
ref_methods = []
ref_methods.append (["omega","omega","openloops","omega"])
ref_methods.append (["openloops","omega","openloops","omega"])
ref_methods.append (["omega","openloops","openloops","omega"])
ref_methods.append (["openloops","openloops","openloops","omega"])
ref_methods.append (["omega","omega","openloops","openloops"])
ref_methods.append (["openloops","omega","openloops","openloops"])
ref_methods.append (["omega","openloops","openloops","openloops"])
ref_methods.append (["openloops","openloops","openloops","openloops"])

with open(filename, 'r') as infile:
  read = False
  for line in infile:
    if read:
        if (k == 4):
           methods.append (tmp)
           read = False 
        m = re.search("\[([A-Za-z0-9_]+)\]", line)
        tmp.append (m.group(1))
        k += 1

    if 'Process components' in line:
        read = True
        tmp = []
        k = 1
try:
    i = 0
    valid = True
    for integral in integrals:
        error = errors[i]
        error_sum = sqrt (error**2 + ref_error**2)
        pull = abs (integral - ref_integral) / error_sum
        valid = valid and pull < 3
        i += 1
    i = 0
    for method in methods:
        ref_method = ref_methods[i]
        valid = valid and (ref_method == method)
        i += 1

except NameError:
  print 'Number not found. Run failed.'
  valid = False

returncode = 0 if valid else 1
sys.exit(returncode)
