#! /usr/bin/env python
# Copyright (c) 2011, 2019 Harald Klimach <harald@klimachs.de>
# Copyright (c) 2014 Peter Vitt <peter.vitt2@uni-siegen.de>
#
# Parts of this file were written by Harald Klimach and Peter Vitt
# for University of Siegen.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
# OR OTHER DEALINGS IN THE SOFTWARE.
# **************************************************************************** #
import re
import ast

TASK_REGEX = "^\s*\{task\s*(\d*):\s*fc\s*(\w*)\.f90\s*->[^\[]*(\[.*)$"

class FortranCircular(Exception):
    """ An exception to identify circular Fortran dependencies. """
    def __init__(self, depstring):
        self.depstring = depstring

    def __str__(self):
        return('Circular dependency in modules encountered: {0}'.format(self.depstring))

def find_circle(start, nodename, nodelist):
    """ Find a circular dependency recursively for the given nodename.
        The "use" chain of modules leading this module has to be provided
        in start, and the list of all modules with their direct dependencies
        to consider in nodelist (dictionary with node name as key and list
        of dependencies as value).
    """
    # First we have to check whether the module is one we found during the file
    # parse. If not, then we cannot detect any circle. If we know the module,
    # proceed with the search.
    if nodename in nodelist:
        # We know the node, so we can check all the used modules
        for child in nodelist[nodename]:
            # Is the node already in the start list?
            if child in start:
                # If so, we have a circular dependency
                childidx = start.index(child)
                depstring = '{0} => {1}'.format(
                      ' -> '.join(map(str,start[childidx:])), child)
                return(depstring)

            # We haven't passed the current node yet, so we have to check it's
            # decendants
            new_start = list(start)
            new_start.append(child)
            childdeps = find_circle(new_start, child, nodelist)
            if childdeps:
                return(childdeps)

    return False

def find_circular_dependency(msg):
    """ Find the circular dependencies for the multiline error message produced
        by waf's exception when encountering a dependency deadlock.
    """
    nodes = dict()
    tasklist = {}
    re_task = re.compile(TASK_REGEX, re.I)

    for line in msg.splitlines():

        m = re_task.search(line)
        if m:
            taskid = ast.literal_eval(m.group(1))
            taskmod = m.group(2)
            taskdeps = ast.literal_eval(m.group(3))
            tasklist[taskid] = {'module': taskmod, 'deps': taskdeps}

    for tid, tsk in tasklist.items():
        depnames = []
        for dep in tsk['deps']:
            if dep in tasklist.keys():
                 depnames.append(tasklist[dep]['module'])
        nodes[tsk['module']] = depnames

    for node in nodes.keys():
        circdep = find_circle([node], node, nodes)
        if circdep:
            raise(FortranCircular(circdep))

    return(False)
