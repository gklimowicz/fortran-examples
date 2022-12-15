#! /usr/bin/env pyhton
# encoding: utf-8
# Copyright (c) 2016-2017, 2018 Harald Klimach <harald@klimachs.de>
# Copyright (c) 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
# Copyright (c) 2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
# Copyright (c) 2018 Raphael Haupt <Raphael.Haupt@student.uni-siegen.de>
#
# Parts of this file were written by Harald Klimach, Peter Vitt, Nikhil Anand
# and Raphael Haupt for University of Siegen.
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

### Definition of some default fortran flags for various compilers.
### Easily select and combine any of the desired features.
### The first key-element has to specify the compiler name as set
### by the waflib/extras/fc_* tools.
### * BGXLF (XLF on BlueGene)
### * CRAY (Cray Compiler)
### * GFORTRAN (Gnu Fortran Compiler from the GCC)
### * IFORT (Intel Fortran Compiler)
### * NAG (NAG Fortran Compiler)
### * NFORT (NEC Fortran Compiler)
### * FUJITSU (Fujitsu Fortran Compiler for ARM64FX)
### * OPEN64 (Open64 Compiler)
### * PGFC (PGI Fortran Compiler)
### * SOL (Oracle Solaris Studio Compiler)
### * XLF (IBM Fortran Compiler)

fcopts = {}

# The following sets of flags are provided for all compilers:
# * warn:     activate compile-time warnings
# * w2e:      turn warnings into errors
# * standard: check for standard compliance (Fortran 2003)
# * double:   turn default reals into double precision
# * debug:    activate debug symbols and runtime checks
# * optimize: activate aggressive optimization
# * openmp:   activate openmp support
# * noomp:    deactivate openmp support
# * profile:  activate profiling (with gprof)

fcopts['BENCH', 'warn'] = ['-Wall', '-Wconversion', '-Wimplicit-interface', '-Wunderflow', '-W', '-frange-check']
fcopts['BENCH', 'w2e'] = ['-Werror']
fcopts['BENCH', 'standard'] = ['-std=f2008']
fcopts['BENCH', 'double'] = ['-fdefault-real-8']
fcopts['BENCH', 'debug'] = ['-fbacktrace', '-fcheck=bounds,do,mem,pointer,recursion', '-finit-real=nan', '-ffpe-trap=invalid,zero,overflow', '-g']
fcopts['BENCH', 'optimize'] = ['-O3', '-march=native']
fcopts['BENCH', 'openmp'] = ['-fopenmp']
fcopts['BENCH', 'noomp'] = []
fcopts['BENCH', 'pre'] = []
fcopts['BENCH', 'profile'] = ['-pg']
fcopts['BENCH', 'fixform'] = ['-ffixed-form']
fcopts['BENCH', 'freeform'] = ['-ffree-form']


fcopts['GFORTRAN', 'warn'] = ['-Wall', '-Wno-maybe-uninitialized' ]
fcopts['GFORTRAN', 'w2e'] = ['-Werror']
fcopts['GFORTRAN', 'standard'] = ['-std=f2008', '-pedantic-errors']
fcopts['GFORTRAN', 'double'] = ['-fdefault-real-8']
fcopts['GFORTRAN', 'debug'] = ['-Og', '-fbacktrace', '-fcheck=bounds,do,mem,pointer,recursion', '-finit-real=nan', '-ffpe-trap=invalid,zero,overflow', '-g']
fcopts['GFORTRAN', 'optimize'] = ['-O3', '-march=native']
fcopts['GFORTRAN', 'openmp'] = ['-fopenmp']
fcopts['GFORTRAN', 'noomp'] = []
fcopts['GFORTRAN', 'pre'] = ['-cpp']
fcopts['GFORTRAN', 'profile'] = ['-pg']
fcopts['GFORTRAN', 'fixform'] = ['-ffixed-form']
fcopts['GFORTRAN', 'freeform'] = ['-ffree-form']

fcopts['IFORT', 'warn'] = '-warn all'.split()
fcopts['IFORT', 'w2e'] = '-warn stderrors'.split()
fcopts['IFORT', 'standard'] = ['-stand']
fcopts['IFORT', 'double'] = '-real-size 64'.split()
fcopts['IFORT', 'debug'] = '-check all -check noarg_temp_created'.split() + ['-traceback', '-g']
fcopts['IFORT', 'optimize'] = '-xHost -O3 -ipo -no-prec-div'.split()
fcopts['IFORT', 'openmp'] = ['-qopenmp']
fcopts['IFORT', 'noomp'] = []
fcopts['IFORT', 'pre'] = ['-fpp']
fcopts['IFORT', 'profile'] = ['-qopt-report=5', '-pg']
fcopts['IFORT', 'fixform'] = []
fcopts['IFORT', 'freeform'] = []

fcopts['IFORTwin', 'warn'] = '/warn:all'.split()
fcopts['IFORTwin', 'w2e'] = '/warn:stderrors'.split()
fcopts['IFORTwin', 'standard'] = ['/stand']
fcopts['IFORTwin', 'double'] = '/real-size:64'
fcopts['IFORTwin', 'debug'] = '/check:all, noarg_temp_created'.split() + ['/traceback', '/g']
fcopts['IFORTwin', 'optimize'] = '/QxHost /O3 /Qprec-div-'.split()
fcopts['IFORTwin', 'openmp'] = ['/Qopenmp']
fcopts['IFORTwin', 'noomp'] = []
fcopts['IFORTwin', 'pre'] = ['/fpp']
fcopts['IFORTwin', 'profile'] = ['/Qopt-report:5']
fcopts['IFORTwin', 'fixform'] = []
fcopts['IFORTwin', 'freeform'] = []

fcopts['SOL', 'warn'] = ['-w4']
fcopts['SOL', 'w2e'] = ['-errwarn=%all']
fcopts['SOL', 'standard'] = []
fcopts['SOL', 'double'] = ['-xtypemap=real:64']
fcopts['SOL', 'debug'] = ['-C', '-xcheck', '-traceback', '-g']
fcopts['SOL', 'optimize'] = ['-fast', '-xipo']
fcopts['SOL', 'openmp'] = []
fcopts['SOL', 'noomp'] = []
fcopts['SOL', 'pre'] = ['-cpp']
fcopts['SOL', 'profile'] = ['-pg']
fcopts['SOL', 'fixform'] = []
fcopts['SOL', 'freeform'] = []

fcopts['PGFC', 'warn'] = ['-Minform=inform', '-Minfo=all']
fcopts['PGFC', 'w2e'] = []
fcopts['PGFC', 'standard'] = ['-Mstandard']
fcopts['PGFC', 'double'] = ['-Mr8']
fcopts['PGFC', 'debug'] = ['-Mbounds', '-Mchkptr', '-Mlist', '-traceback', '-g']
fcopts['PGFC', 'optimize'] = ['-O4']
fcopts['PGFC', 'openmp'] = []
fcopts['PGFC', 'noomp'] = []
fcopts['PGFC', 'pre'] = ['-cpp']
fcopts['PGFC', 'profile'] = ['-pg']
fcopts['PGFC', 'fixform'] = []
fcopts['PGFC', 'freeform'] = []

fcopts['BGXLF', 'warn'] = []
fcopts['BGXLF', 'w2e'] = ['-qhalt=w']
fcopts['BGXLF', 'standard'] = ['-qlanglvl=2008pure']
fcopts['BGXLF', 'double'] = ['-qautodbl=dbl4']
fcopts['BGXLF', 'debug'] = ['-C', '-g', '-qflttrap', '-qfullpath']
fcopts['BGXLF', 'optimize'] = ['-O5']
fcopts['BGXLF', 'openmp'] = ['-qsmp']
fcopts['BGXLF', 'noomp'] = []
fcopts['BGXLF', 'pre'] = ['-cpp']
fcopts['BGXLF', 'profile'] = []
fcopts['BGXLF', 'fixform'] = ['-qfixed=72']
fcopts['BGXLF', 'freeform'] = ['-qfree']

fcopts['CRAY', 'warn'] = ['-m0']
fcopts['CRAY', 'w2e'] = []
fcopts['CRAY', 'standard'] = ['-e', 'n']
fcopts['CRAY', 'double'] = ['-s', 'real64']
fcopts['CRAY', 'debug'] = ['-e', 'DcI', '-R', 'bcps']
fcopts['CRAY', 'optimize'] = ['-O3']
fcopts['CRAY', 'openmp'] = ['-h', 'omp']
fcopts['CRAY', 'noomp'] = ['-h', 'noomp']
fcopts['CRAY', 'pre'] = ['-cpp']
fcopts['CRAY', 'profile'] = ['-h', 'profile_generate', '-h', 'func_trace', '-h', 'keepfiles']
fcopts['CRAY', 'fixform'] = []
fcopts['CRAY', 'freeform'] = []

fcopts['NAG', 'warn'] = []
fcopts['NAG', 'w2e'] = []
fcopts['NAG', 'standard'] = ['-f2008']
fcopts['NAG', 'double'] = ['-r8']
fcopts['NAG', 'debug'] = ['-C=all', '-mtrace=all', '-nan', '-gline', '-g', '-g90']
fcopts['NAG', 'optimize'] = ['-O4']
fcopts['NAG', 'openmp'] = ['-openmp']
fcopts['NAG', 'noomp'] = []
fcopts['NAG', 'pre'] = ['-cpp']
fcopts['NAG', 'profile'] = ['-pg']
fcopts['NAG', 'fixform'] = []
fcopts['NAG', 'freeform'] = []

fcopts['NEC', 'warn'] = ['-w,all']
fcopts['NEC', 'w2e'] = []
fcopts['NEC', 'standard'] = ['-f2008']
fcopts['NEC', 'double'] = []
fcopts['NEC', 'debug'] = ['-check,all','-traceback','-init,stack=nan','-mtrace,full','-R,diaglist']
fcopts['NEC', 'optimize'] = ['-C,hopt']
fcopts['NEC', 'openmp'] = ['-openmp']
fcopts['NEC', 'noomp'] = []
fcopts['NEC', 'pre'] = ['-cpp']
fcopts['NEC', 'profile'] = ['-ftrace','-O,fullmsg','-pvctl,fullmsg','-R,fmtlist,diaglist,summary']
fcopts['NEC', 'fixform'] = []
fcopts['NEC', 'freeform'] = []

fcopts['NFORT', 'warn'] = ['-Wall', '-Woverflow-errors', '-Wextension', '-Wobsolescent']
fcopts['NFORT', 'w2e'] = ['-Werror']
fcopts['NFORT', 'standard'] = ['-std=f2008']
fcopts['NFORT', 'double'] = []
fcopts['NFORT', 'debug'] = ['-g', '-fcheck=all', '-traceback', '-minit-stack=nan', '-mmemory-trace-full']
fcopts['NFORT', 'optimize'] = ['-O4']
fcopts['NFORT', 'openmp'] = ['-fopenmp']
fcopts['NFORT', 'noomp'] = []
fcopts['NFORT', 'pre'] = ['-fpp']
fcopts['NFORT', 'profile'] = ['-ftrace','-fdiag-vector=2','-report-all']
fcopts['NFORT', 'fixform'] = ['-ffixed-form']
fcopts['NFORT', 'freeform'] = ['-ffree-form']

fcopts['FUJITSU', 'warn'] = ['-fw', '-N maxserious=1']
fcopts['FUJITSU', 'w2e'] = []
fcopts['FUJITSU', 'standard'] = ['-N check_std=08s']
fcopts['FUJITSU', 'double'] = []
fcopts['FUJITSU', 'debug'] = ['-g', '-NRtrap', '-Haefosux']
fcopts['FUJITSU', 'optimize'] = ['-Kfast']
fcopts['FUJITSU', 'openmp'] = ['-Kopenmp']
fcopts['FUJITSU', 'noomp'] = []
fcopts['FUJITSU', 'pre'] = []
fcopts['FUJITSU', 'profile'] = ['-Nlst=t']
fcopts['FUJITSU', 'fixform'] = ['-Fixed']
fcopts['FUJITSU', 'freeform'] = ['-Free']

### End of set of Fortran flags
#########################################################################

# Helpers to manage the flags:
def options(opt):
    opt.load('compiler_fc')
    fcopts = opt.add_option_group('Fortran Compiler specific settings')
    fcopts.add_option('--fc_addflags', action='store',
                      help='Additional Fortran compiler flags to use.')

    fcopts.add_option('--fc_delflags', action='store',
                      help='Fortran compiler flags to remove.')

    fcopts.add_option('--bench', action='store_true', default=False,
                      help='Set this flag to use the explicit compiler given in bench_compiler.py')


def configure(conf):
    from waflib import Logs
    if not hasattr(conf.options, 'bench'):
       Logs.error("Error: use opt.load('fortran_compiler') in your options!")
    if not conf.options.bench:
      conf.load('compiler_fc')

    else:
      print('loading bench compiler')
      from bench_compiler import set_bench_f_compiler
      set_bench_f_compiler(conf)
      conf.load('fc')


def set_fc_flags(conf, flagset, osflags=None):
    ''' Define a set of flags to use for Fortran compilations.
        The flagset has to be an array of strings indicating the
        set of options to pick from the fcopts table.
    '''
    myflags = list(osflags) or []
    if getattr(conf.env, 'IFORT_WIN32', False):
      fcname = 'IFORTwin'
    else:
      fcname = conf.env.FC_NAME

    # Flag sets to add from the fcopts table.
    for fs in flagset:
        myflags += fcopts[fcname, fs]
        # Add the sanitize option to the debug flags for gfortran >= 4.8
        if (fs == 'debug') and (fcname == 'GFORTRAN'):
            if ( (int(conf.env.FC_VERSION[0]) == 4)
                 and (int(conf.env.FC_VERSION[1]) >= 8)
                 or int(conf.env.FC_VERSION[0]) > 4):
                myflags += ['-fsanitize=address']


    # Additional flags to set due to the fc_addflags option.
    if conf.options.fc_addflags:
        myflags += conf.options.fc_addflags.split(' ')

    # Flags to filter out according to the fc_delflags option.
    if conf.options.fc_delflags:
        delflags = conf.options.fc_delflags.split(' ')
    else:
        delflags = []

    # Now set the FCFLAGS in the environment accordingly.
    conf.env.FCFLAGS = []
    for flag in myflags:
        if not flag in delflags:
            conf.env.append_value('FCFLAGS', [flag])
    if getattr(conf.env, 'LINK_FC', conf.env.FC) == conf.env.FC:
        conf.env.LINKFLAGS = conf.env.FCFLAGS

from waflib import TaskGen

@TaskGen.feature('fc')
@TaskGen.after_method('apply_link')
def addremove_fc_flags(self):
  if hasattr(self, 'delflags'):
    delflags = self.delflags
  else:
    delflags = []
  self.env.detach()
  if hasattr(self, 'bld'):
    if hasattr(self.bld, 'options'):
      if self.bld.options.fc_addflags:
        for flag in self.bld.options.fc_addflags.split(' '):
          self.env.append_unique('FCFLAGS', [flag])
      if self.bld.options.fc_delflags:
        delflags += self.bld.options.fc_delflags.split(' ')
  for flag in delflags:
    try:
      self.env.FCFLAGS.remove(flag)
    except ValueError:
      pass

@TaskGen.feature('fcprogram', 'fcshlib', 'fcstlib')
@TaskGen.after_method('apply_link')
def addremove_fclink_flags(self):
  if hasattr(self, 'bld'):
    if hasattr(self.bld, 'options'):
      self.env.detach()
      if self.bld.options.fc_addflags:
        for flag in self.bld.options.fc_addflags.split(' '):
          self.env.append_unique('LINKFLAGS', [flag])
      if self.bld.options.fc_delflags:
        for flag in self.bld.options.fc_delflags.split(' '):
          try:
            self.env.LINKFLAGS.remove(flag)
          except ValueError:
            pass
