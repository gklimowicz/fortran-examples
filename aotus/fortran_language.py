# Copyright (c) 2016, 2018-2019 Harald Klimach <harald@klimachs.de>
#
# Parts of this file were written by Harald Klimach for University of Siegen.
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

# Checks for supported Fortran language features
#
# Results are stored in the context environment with a fortsupp_ prefix.
#
def options(opt):
    flopts = opt.add_option_group('Fortran Language specific settings')
    flopts.add_option('--quad_kind', action='store',
                      help='Kind for quadruple precision reals')
    flopts.add_option('--xdble_kind', action='store',
                      help='Kind for extended double precision reals')
    flopts.add_option('--fortran_isNaN', action='store',
                      help='Function name to use for vendor specific isNaN checking.')


def supports_iso_c(conf, mandatory=True):
  '''
     Check whether the Fotran compiler from the context in conf supports the
     Fortran 2003 ISO_C_Binding.
     Set mandatory to False if the config should not abort upon missing
     ISO_C_Binding support.
     Result stored in conf.env.fortsupp_iso_c
  '''

  fcenv = conf.env.derive()
  fcenv.detach()

  conf.check_fc( fragment = '''
program check_iso_c
  use, intrinsic :: iso_c_binding
  implicit none
  write(*,*) c_int
end program check_iso_c''',
                 msg = "Checking for ISO_C_Binding support",
                 define_name = 'iso_c_binding',
                 mandatory = mandatory)

  fcenv['fortsupp_iso_c'] = conf.is_defined('iso_c_binding')
  conf.env = fcenv


################
# isNaN support:
################
isNaN_vendor_stub = '''
program check_isnan
  implicit none
  real :: a_real
  logical :: noNumber
  noNumber = {0}(a_real)
  write(*,*) noNumber
end program check_isnan
'''

def supports_ieee_is_NaN(conf, mandatory=True):
  '''
     Check for IEEE isNaN support.
  '''
  fcenv = conf.env.derive()
  fcenv.detach()

  conf.check_fc( fragment = '''
program check_isnan
  use, intrinsic :: ieee_arithmetic
  implicit none
  real :: a_real
  logical :: noNumber
  noNumber = ieee_is_NaN(a_real)
  write(*,*) noNumber
end program check_isnan''',
                 msg = 'Checking for IEEE_is_NaN',
                 mandatory = mandatory, define_name='isNaN')
  fcenv['fortsupp_ieee_is_NaN'] = conf.is_defined('isNaN')

  conf.env = fcenv

def supports_vendor_is_NaN(conf, mandatory=True):
  '''
     Check for vendor specific isNaN support.
  '''
  fcenv = conf.env.derive()
  fcenv.detach()

  isNaN_name = 'isNaN'
  if conf.options.fortran_isNaN:
    isNaN_name = conf.options.fortran_isNaN
  conf.check_fc(fragment = isNaN_vendor_stub.format(isNaN_name),
                msg = 'Checking for vendor specific isNaN',
                mandatory = mandatory, define_name='isNaN')
  if conf.is_defined('isNaN'):
    fcenv['fortsupp_vendor_is_NaN'] = isNaN_name
  else:
    fcenv['fortsupp_vendor_is_NaN'] = None

  conf.env = fcenv


def supports_c_sizeof(conf, mandatory=True):
  '''
     Check for F2008 c_sizeof support.
  '''
  fcenv = conf.env.derive()
  fcenv.detach()

  conf.check_fc( fragment = '''
program check_c_sizeof
  use, intrinsic :: iso_c_binding, only: c_sizeof, c_float
  implicit none
  real(kind=c_float) :: a_real
  integer :: realsize
  realsize = c_sizeof(a_real)
  write(*,*) realsize
end program check_c_sizeof
''',
                 msg = 'Checking for F2008 c_sizeof support',
                 mandatory = mandatory, define_name='c_sizeof')
  fcenv['fortsupp_c_sizeof'] = conf.is_defined('c_sizeof')

  conf.env = fcenv

def supports_vendor_sizeof(conf, mandatory=True):
  '''
     Check for vendor sizeof support.
  '''
  fcenv = conf.env.derive()
  fcenv.detach()

  conf.check_fc( fragment = '''
program check_vendor_sizeof
  implicit none
  real :: a_real
  integer :: realsize
  realsize = sizeof(a_real)
  write(*,*) realsize
end program check_vendor_sizeof
''',
                 msg = 'Checking for vendor sizeof support',
                 mandatory = mandatory, define_name='vendor_sizeof')
  fcenv['fortsupp_vendor_sizeof'] = conf.is_defined('vendor_sizeof')

  conf.env = fcenv

##########################
# F2008 special functions:
##########################

def supports_f2008_gamma(conf, mandatory=True):
  '''
     Check for F2008 gamma function support.
  '''

  fcenv = conf.env.derive()
  fcenv.detach()

  conf.check_fc(fragment = '''
program check_gamma
  implicit none
  real :: a_real
  real :: res
  read(*,*) a_real
  res = gamma(a_real)
  write(*,*) res
end program check_gamma''',
                  msg = 'Checking for F2008 Gamma function',
                  mandatory = mandatory, define_name='has_f2008_gamma')
  fcenv['fortsupp_f2008_gamma'] = conf.is_defined('has_f2008_gamma')
  conf.env = fcenv

def supports_f2008_bessel(conf, mandatory=True):
  '''
     Check for F2008 bessel function support.
  '''

  fcenv = conf.env.derive()
  fcenv.detach()

  conf.check_fc(fragment = '''
program check_bessel
  implicit none
  real :: a_real
  real :: res
  integer :: an_int
  read(*,*) an_int
  read(*,*) a_real
  res = bessel_jn(an_int, a_real)
  write(*,*) res
end program check_bessel''',
                msg = 'Checking for F2008 Bessel_jn function',
                mandatory = mandatory, define_name = 'has_f2008_bessel')
  fcenv['fortsupp_f2008_bessel'] = conf.is_defined('has_f2008_bessel')

  conf.env = fcenv


################
# Integer kinds:
################

int_prec_stub = '''
program check_i_kind
  implicit none
  integer, parameter :: int_k = selected_int_kind({0})
  integer(kind=int_k) :: an_int
  write(*,*) int_k
end program check_i_kind
'''

default_int_stub = '''
program check_default_ikind
  implicit none
  write(*,*) kind(1)
end program check_default_ikind
'''

def supported_int_kinds(conf, precisions = [ ( 1, True),
                                             ( 2, True),
                                             ( 4, True),
                                             ( 8, True),
                                             ( 9, True),
                                             (16, False),
                                             (18, False),
                                             (32, False),
                                             (38, False) ]):
  ''' Check the available integer kinds.
      Results in a dictionary that maps the tested precisions
      to their respective integer kinds if supported.
      The precisions needs to be a list of tuples, where each
      tuple indicates a precision to check (as first index) and
      whether it is mandatory or not (as second index).
      By default, powers of two precisions are checked up to 32
      and the usual maximal decimal places that fit into 4 byte (9),
      8 byte (18) and 16 byte (38).
      9 supported decimal places (4 byte integers) are assumed
      mandatory.
  '''
  conf.env['fortsupp_intkind'] = {}
  for prec in precisions:
    fcenv = conf.env.derive()
    fcenv.detach()

    conf.check_fc( fragment = int_prec_stub.format(prec[0]),
                   mandatory = prec[1],
                   msg = 'Checking integer kind with {0} decimal places.'.format(prec[0]),
                   define_name = 'intkind',
                   execute = True, define_ret = True)

    if conf.is_defined('intkind'):
      fcenv['fortsupp_intkind'][prec[0]] = int( conf.get_define('intkind')
                                             .replace('"', '').strip() )
    conf.env = fcenv

def default_int_kind(conf):
  fcenv = conf.env.derive()
  fcenv.detach()

  conf.check_fc( fragment = default_int_stub,
                 mandatory = False,
                 msg = 'Getting default integer kind.',
                 define_name = 'defintkind',
                 execute = True, define_ret = True)

  return( int( conf.get_define('defintkind').replace('"', '').strip() ) )


#############
# Real kinds:
#############

real_prec_stub = '''
program check_r_kind
  implicit none
  integer, parameter :: real_k = selected_real_kind({0})
  real(kind=real_k) :: a_real
  write(*,*) real_k
end program check_r_kind
'''

real_kind_stub = '''
program check_r_kind
  implicit none
  real(kind={0}) :: a_real
  write(*,*) {0}
end program check_r_kind
'''

def supports_quad_kind(conf, mandatory=True):
  '''
     Check whether the Fotran compiler from the context in conf supports the
     quadruple precision real numbers.
     Set mandatory to False if the config should not abort upon missing
     quadruple support.
     Resulting kind value for quadruple reals is stored in
     conf.env.fortsupp_quad_kind (negative if not available)
  '''

  conf.start_msg('Checking for Quadruple precision')

  fcenv = conf.env.derive()
  fcenv.detach()

  if conf.options.quad_kind:
    conf.check_fc( fragment = real_kind_stub.format(conf.options.quad_kind),
                   mandatory = mandatory,
                   define_name = 'quadruple')
    if conf.is_defined('quadruple'):
      fcenv['fortsupp_quad_kind'] = conf.options.quad_kind
    else:
      fcenv['fortsupp_quad_kind'] = -1
  else:
    conf.check_fc( fragment = real_prec_stub.format(33),
                   mandatory = mandatory,
                   define_name = 'quadruple',
                   execute = True, define_ret = True)

    if conf.is_defined('quadruple'):
      fcenv['fortsupp_quad_kind'] = int( conf.get_define('quadruple')
                                             .replace('"', '').strip() )
    else:
      fcenv['fortsupp_quad_kind'] = -1

  conf.env = fcenv

  if conf.env['fortsupp_quad_kind'] > 0:
     conf.end_msg('yes', color='GREEN')
  else:
     conf.end_msg('NO', color='RED')


def supports_xdble_kind(conf, mandatory=True):
  '''
     Check whether the Fotran compiler from the context in conf supports the
     extended double precision real numbers.
     Set mandatory to False if the config should not abort upon missing
     extended double support.
     Resulting kind value for extended double precision reals is stored in
     conf.env.fortsupp_xdble_kind (negative if not available)
  '''

  # Make sure to look for quadruple precision first:
  if not hasattr(conf.env, 'fortsupp_quad_kind'):
    supports_quadruple_kind(conf)

  conf.start_msg('Checking for Extended double precision')

  fcenv = conf.env.derive()
  fcenv.detach()

  if conf.options.xdble_kind:
    conf.check_fc( fragment = real_kind_stub.format(conf.options.xdble_kind),
                   mandatory = mandatory,
                   define_name = 'xdble')
    if conf.is_defined('xdble'):
      fcenv['fortsupp_xdble_kind'] = conf.options.xdble_kind
    else:
      fcenv['fortsupp_xdble_kind'] = -1
  else:

    conf.check_fc( fragment = real_prec_stub.format(18),
                   mandatory = mandatory,
                   define_name = 'xdble',
                   execute = True, define_ret = True)

    if conf.is_defined('xdble'):
      fcenv['fortsupp_xdble_kind'] = int( conf.get_define('xdble')
                                              .replace('"', '').strip() )
    else:
      fcenv['fortsupp_xdble_kind'] = -1

  # If the found kind is actually the quadruple precision, do not set the
  # extended double precision kind.
  if fcenv['fortsupp_xdble_kind'] == fcenv['fortsupp_quad_kind']:
    fcenv['fortsupp_xdble_kind'] = -1
  conf.env = fcenv

  if conf.env['fortsupp_xdble_kind'] > 0:
     conf.end_msg('yes', color='GREEN')
  else:
     conf.end_msg('NO', color='RED')

