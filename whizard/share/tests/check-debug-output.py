import re, sys

re_e = re.compile("p\(0:3\) = *([0-9]+\.[0-9]+)")
def energy(line):
  try:
    return float(re_e.search(line).group(1))
  except AttributeError:
    print 'Could not find number in this line:'
    print line
    return 0.0
test = sys.argv[1]
filename = test + '.debug'
inc_energy = 0.0
out_energy = 0.0
line_no = 0
rel_delta = 10.0**-5
abs_delta = 0.002
valid = True
def nearly_equal(a, b):
  try:
    diff = abs(a - b)
    return diff / max([a, b]) < rel_delta or diff < abs_delta
  except TypeError:
    return False

print 'Start file ' + filename
print (3 * '{:<20s}').format('Line-Nr', 'Incoming-Energy', 'Outgoing-Energy')
with open(filename, 'r') as infile:
  for line in infile:
    line_no += 1
    if 'incoming momenta' in line:
      inc_energy = energy(line)
    elif 'outgoing momenta' in line:
      out_energy = energy(line)
      if not nearly_equal(out_energy, inc_energy):
        print ("{:<20d}" + 2 *"{:<20.10f}").format(line_no, inc_energy, out_energy)
        valid = False
print 'End file ' + filename
returncode = 0 if valid else 1
print 'returncode:', returncode
sys.exit(returncode)
