# This file was automatically created by FeynRules 1.7.214
# Mathematica version: 9.0 for Mac OS X x86 (64-bit) (January 24, 2013)
# Date: Mon 26 Aug 2013 03:14:07



from object_library import all_parameters, Parameter


from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot

# This is a default parameter object representing 0.
ZERO = Parameter(name = 'ZERO',
                 nature = 'internal',
                 type = 'real',
                 value = '0.0',
                 texname = '0')

# User-defined parameters.
RRd1x1 = Parameter(name = 'RRd1x1',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRd1x1}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 1, 1 ])

RRd2x2 = Parameter(name = 'RRd2x2',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRd2x2}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 2, 2 ])

RRd3x3 = Parameter(name = 'RRd3x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.938737896,
                   texname = '\\text{RRd3x3}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 3, 3 ])

RRd3x6 = Parameter(name = 'RRd3x6',
                   nature = 'external',
                   type = 'real',
                   value = 0.344631925,
                   texname = '\\text{RRd3x6}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 3, 6 ])

RRd4x4 = Parameter(name = 'RRd4x4',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRd4x4}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 4, 4 ])

RRd5x5 = Parameter(name = 'RRd5x5',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRd5x5}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 5, 5 ])

RRd6x3 = Parameter(name = 'RRd6x3',
                   nature = 'external',
                   type = 'real',
                   value = -0.344631925,
                   texname = '\\text{RRd6x3}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 6, 3 ])

RRd6x6 = Parameter(name = 'RRd6x6',
                   nature = 'external',
                   type = 'real',
                   value = 0.938737896,
                   texname = '\\text{RRd6x6}',
                   lhablock = 'DSQMIX',
                   lhacode = [ 6, 6 ])

alp = Parameter(name = 'alp',
                nature = 'external',
                type = 'real',
                value = -0.11382521,
                texname = '\\alpha',
                lhablock = 'FRALPHA',
                lhacode = [ 1 ])

RMUH = Parameter(name = 'RMUH',
                 nature = 'external',
                 type = 'real',
                 value = 357.680977,
                 texname = '\\text{RMUH}',
                 lhablock = 'HMIX',
                 lhacode = [ 1 ])

tb = Parameter(name = 'tb',
               nature = 'external',
               type = 'real',
               value = 9.74862403,
               texname = 't_b',
               lhablock = 'HMIX',
               lhacode = [ 2 ])

RmD21x1 = Parameter(name = 'RmD21x1',
                    nature = 'external',
                    type = 'real',
                    value = 273684.674,
                    texname = '\\text{RmD21x1}',
                    lhablock = 'MSD2',
                    lhacode = [ 1, 1 ])

RmD22x2 = Parameter(name = 'RmD22x2',
                    nature = 'external',
                    type = 'real',
                    value = 273684.674,
                    texname = '\\text{RmD22x2}',
                    lhablock = 'MSD2',
                    lhacode = [ 2, 2 ])

RmD23x3 = Parameter(name = 'RmD23x3',
                    nature = 'external',
                    type = 'real',
                    value = 270261.969,
                    texname = '\\text{RmD23x3}',
                    lhablock = 'MSD2',
                    lhacode = [ 3, 3 ])

RmE21x1 = Parameter(name = 'RmE21x1',
                    nature = 'external',
                    type = 'real',
                    value = 18630.6287,
                    texname = '\\text{RmE21x1}',
                    lhablock = 'MSE2',
                    lhacode = [ 1, 1 ])

RmE22x2 = Parameter(name = 'RmE22x2',
                    nature = 'external',
                    type = 'real',
                    value = 18630.6287,
                    texname = '\\text{RmE22x2}',
                    lhablock = 'MSE2',
                    lhacode = [ 2, 2 ])

RmE23x3 = Parameter(name = 'RmE23x3',
                    nature = 'external',
                    type = 'real',
                    value = 17967.6406,
                    texname = '\\text{RmE23x3}',
                    lhablock = 'MSE2',
                    lhacode = [ 3, 3 ])

RmL21x1 = Parameter(name = 'RmL21x1',
                    nature = 'external',
                    type = 'real',
                    value = 38155.67,
                    texname = '\\text{RmL21x1}',
                    lhablock = 'MSL2',
                    lhacode = [ 1, 1 ])

RmL22x2 = Parameter(name = 'RmL22x2',
                    nature = 'external',
                    type = 'real',
                    value = 38155.67,
                    texname = '\\text{RmL22x2}',
                    lhablock = 'MSL2',
                    lhacode = [ 2, 2 ])

RmL23x3 = Parameter(name = 'RmL23x3',
                    nature = 'external',
                    type = 'real',
                    value = 37828.6769,
                    texname = '\\text{RmL23x3}',
                    lhablock = 'MSL2',
                    lhacode = [ 3, 3 ])

RMx1 = Parameter(name = 'RMx1',
                 nature = 'external',
                 type = 'real',
                 value = 101.396534,
                 texname = '\\text{RMx1}',
                 lhablock = 'MSOFT',
                 lhacode = [ 1 ])

RMx2 = Parameter(name = 'RMx2',
                 nature = 'external',
                 type = 'real',
                 value = 191.504241,
                 texname = '\\text{RMx2}',
                 lhablock = 'MSOFT',
                 lhacode = [ 2 ])

RMx3 = Parameter(name = 'RMx3',
                 nature = 'external',
                 type = 'real',
                 value = 588.263031,
                 texname = '\\text{RMx3}',
                 lhablock = 'MSOFT',
                 lhacode = [ 3 ])

mHd2 = Parameter(name = 'mHd2',
                 nature = 'external',
                 type = 'real',
                 value = 32337.4943,
                 texname = '\\text{Subsuperscript}\\left[m,H_d,2\\right]',
                 lhablock = 'MSOFT',
                 lhacode = [ 21 ])

mHu2 = Parameter(name = 'mHu2',
                 nature = 'external',
                 type = 'real',
                 value = -128800.134,
                 texname = '\\text{Subsuperscript}\\left[m,H_u,2\\right]',
                 lhablock = 'MSOFT',
                 lhacode = [ 22 ])

RmQ21x1 = Parameter(name = 'RmQ21x1',
                    nature = 'external',
                    type = 'real',
                    value = 299836.701,
                    texname = '\\text{RmQ21x1}',
                    lhablock = 'MSQ2',
                    lhacode = [ 1, 1 ])

RmQ22x2 = Parameter(name = 'RmQ22x2',
                    nature = 'external',
                    type = 'real',
                    value = 299836.701,
                    texname = '\\text{RmQ22x2}',
                    lhablock = 'MSQ2',
                    lhacode = [ 2, 2 ])

RmQ23x3 = Parameter(name = 'RmQ23x3',
                    nature = 'external',
                    type = 'real',
                    value = 248765.367,
                    texname = '\\text{RmQ23x3}',
                    lhablock = 'MSQ2',
                    lhacode = [ 3, 3 ])

RmU21x1 = Parameter(name = 'RmU21x1',
                    nature = 'external',
                    type = 'real',
                    value = 280382.106,
                    texname = '\\text{RmU21x1}',
                    lhablock = 'MSU2',
                    lhacode = [ 1, 1 ])

RmU22x2 = Parameter(name = 'RmU22x2',
                    nature = 'external',
                    type = 'real',
                    value = 280382.106,
                    texname = '\\text{RmU22x2}',
                    lhablock = 'MSU2',
                    lhacode = [ 2, 2 ])

RmU23x3 = Parameter(name = 'RmU23x3',
                    nature = 'external',
                    type = 'real',
                    value = 179137.072,
                    texname = '\\text{RmU23x3}',
                    lhablock = 'MSU2',
                    lhacode = [ 3, 3 ])

RNN1x1 = Parameter(name = 'RNN1x1',
                   nature = 'external',
                   type = 'real',
                   value = 0.98636443,
                   texname = '\\text{RNN1x1}',
                   lhablock = 'NMIX',
                   lhacode = [ 1, 1 ])

RNN1x2 = Parameter(name = 'RNN1x2',
                   nature = 'external',
                   type = 'real',
                   value = -0.0531103553,
                   texname = '\\text{RNN1x2}',
                   lhablock = 'NMIX',
                   lhacode = [ 1, 2 ])

RNN1x3 = Parameter(name = 'RNN1x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.146433995,
                   texname = '\\text{RNN1x3}',
                   lhablock = 'NMIX',
                   lhacode = [ 1, 3 ])

RNN1x4 = Parameter(name = 'RNN1x4',
                   nature = 'external',
                   type = 'real',
                   value = -0.0531186117,
                   texname = '\\text{RNN1x4}',
                   lhablock = 'NMIX',
                   lhacode = [ 1, 4 ])

RNN2x1 = Parameter(name = 'RNN2x1',
                   nature = 'external',
                   type = 'real',
                   value = 0.0993505358,
                   texname = '\\text{RNN2x1}',
                   lhablock = 'NMIX',
                   lhacode = [ 2, 1 ])

RNN2x2 = Parameter(name = 'RNN2x2',
                   nature = 'external',
                   type = 'real',
                   value = 0.944949299,
                   texname = '\\text{RNN2x2}',
                   lhablock = 'NMIX',
                   lhacode = [ 2, 2 ])

RNN2x3 = Parameter(name = 'RNN2x3',
                   nature = 'external',
                   type = 'real',
                   value = -0.26984672,
                   texname = '\\text{RNN2x3}',
                   lhablock = 'NMIX',
                   lhacode = [ 2, 3 ])

RNN2x4 = Parameter(name = 'RNN2x4',
                   nature = 'external',
                   type = 'real',
                   value = 0.156150698,
                   texname = '\\text{RNN2x4}',
                   lhablock = 'NMIX',
                   lhacode = [ 2, 4 ])

RNN3x1 = Parameter(name = 'RNN3x1',
                   nature = 'external',
                   type = 'real',
                   value = -0.0603388002,
                   texname = '\\text{RNN3x1}',
                   lhablock = 'NMIX',
                   lhacode = [ 3, 1 ])

RNN3x2 = Parameter(name = 'RNN3x2',
                   nature = 'external',
                   type = 'real',
                   value = 0.0877004854,
                   texname = '\\text{RNN3x2}',
                   lhablock = 'NMIX',
                   lhacode = [ 3, 2 ])

RNN3x3 = Parameter(name = 'RNN3x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.695877493,
                   texname = '\\text{RNN3x3}',
                   lhablock = 'NMIX',
                   lhacode = [ 3, 3 ])

RNN3x4 = Parameter(name = 'RNN3x4',
                   nature = 'external',
                   type = 'real',
                   value = 0.710226984,
                   texname = '\\text{RNN3x4}',
                   lhablock = 'NMIX',
                   lhacode = [ 3, 4 ])

RNN4x1 = Parameter(name = 'RNN4x1',
                   nature = 'external',
                   type = 'real',
                   value = -0.116507132,
                   texname = '\\text{RNN4x1}',
                   lhablock = 'NMIX',
                   lhacode = [ 4, 1 ])

RNN4x2 = Parameter(name = 'RNN4x2',
                   nature = 'external',
                   type = 'real',
                   value = 0.310739017,
                   texname = '\\text{RNN4x2}',
                   lhablock = 'NMIX',
                   lhacode = [ 4, 2 ])

RNN4x3 = Parameter(name = 'RNN4x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.64922596,
                   texname = '\\text{RNN4x3}',
                   lhablock = 'NMIX',
                   lhacode = [ 4, 3 ])

RNN4x4 = Parameter(name = 'RNN4x4',
                   nature = 'external',
                   type = 'real',
                   value = -0.684377823,
                   texname = '\\text{RNN4x4}',
                   lhablock = 'NMIX',
                   lhacode = [ 4, 4 ])

RRl1x1 = Parameter(name = 'RRl1x1',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRl1x1}',
                   lhablock = 'SELMIX',
                   lhacode = [ 1, 1 ])

RRl2x2 = Parameter(name = 'RRl2x2',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRl2x2}',
                   lhablock = 'SELMIX',
                   lhacode = [ 2, 2 ])

RRl3x3 = Parameter(name = 'RRl3x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.28248719,
                   texname = '\\text{RRl3x3}',
                   lhablock = 'SELMIX',
                   lhacode = [ 3, 3 ])

RRl3x6 = Parameter(name = 'RRl3x6',
                   nature = 'external',
                   type = 'real',
                   value = 0.959271071,
                   texname = '\\text{RRl3x6}',
                   lhablock = 'SELMIX',
                   lhacode = [ 3, 6 ])

RRl4x4 = Parameter(name = 'RRl4x4',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRl4x4}',
                   lhablock = 'SELMIX',
                   lhacode = [ 4, 4 ])

RRl5x5 = Parameter(name = 'RRl5x5',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRl5x5}',
                   lhablock = 'SELMIX',
                   lhacode = [ 5, 5 ])

RRl6x3 = Parameter(name = 'RRl6x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.959271071,
                   texname = '\\text{RRl6x3}',
                   lhablock = 'SELMIX',
                   lhacode = [ 6, 3 ])

RRl6x6 = Parameter(name = 'RRl6x6',
                   nature = 'external',
                   type = 'real',
                   value = -0.28248719,
                   texname = '\\text{RRl6x6}',
                   lhablock = 'SELMIX',
                   lhacode = [ 6, 6 ])

aEWM1 = Parameter(name = 'aEWM1',
                  nature = 'external',
                  type = 'real',
                  value = 127.934,
                  texname = '\\text{Subsuperscript}[\\alpha ,w,-1]',
                  lhablock = 'SMINPUTS',
                  lhacode = [ 1 ])

aS = Parameter(name = 'aS',
               nature = 'external',
               type = 'real',
               value = 0.118,
               texname = '\\alpha _s',
               lhablock = 'SMINPUTS',
               lhacode = [ 3 ])

RRn1x1 = Parameter(name = 'RRn1x1',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRn1x1}',
                   lhablock = 'SNUMIX',
                   lhacode = [ 1, 1 ])

RRn2x2 = Parameter(name = 'RRn2x2',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRn2x2}',
                   lhablock = 'SNUMIX',
                   lhacode = [ 2, 2 ])

RRn3x3 = Parameter(name = 'RRn3x3',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRn3x3}',
                   lhablock = 'SNUMIX',
                   lhacode = [ 3, 3 ])

Rtd3x3 = Parameter(name = 'Rtd3x3',
                   nature = 'external',
                   type = 'real',
                   value = -110.693742,
                   texname = '\\text{Rtd3x3}',
                   lhablock = 'TD',
                   lhacode = [ 3, 3 ])

Rte3x3 = Parameter(name = 'Rte3x3',
                   nature = 'external',
                   type = 'real',
                   value = -25.4019727,
                   texname = '\\text{Rte3x3}',
                   lhablock = 'TE',
                   lhacode = [ 3, 3 ])

Rtu3x3 = Parameter(name = 'Rtu3x3',
                   nature = 'external',
                   type = 'real',
                   value = -444.752457,
                   texname = '\\text{Rtu3x3}',
                   lhablock = 'TU',
                   lhacode = [ 3, 3 ])

RUU1x1 = Parameter(name = 'RUU1x1',
                   nature = 'external',
                   type = 'real',
                   value = 0.916834859,
                   texname = '\\text{RUU1x1}',
                   lhablock = 'UMIX',
                   lhacode = [ 1, 1 ])

RUU1x2 = Parameter(name = 'RUU1x2',
                   nature = 'external',
                   type = 'real',
                   value = -0.399266629,
                   texname = '\\text{RUU1x2}',
                   lhablock = 'UMIX',
                   lhacode = [ 1, 2 ])

RUU2x1 = Parameter(name = 'RUU2x1',
                   nature = 'external',
                   type = 'real',
                   value = 0.399266629,
                   texname = '\\text{RUU2x1}',
                   lhablock = 'UMIX',
                   lhacode = [ 2, 1 ])

RUU2x2 = Parameter(name = 'RUU2x2',
                   nature = 'external',
                   type = 'real',
                   value = 0.916834859,
                   texname = '\\text{RUU2x2}',
                   lhablock = 'UMIX',
                   lhacode = [ 2, 2 ])

RMNS1x1 = Parameter(name = 'RMNS1x1',
                    nature = 'external',
                    type = 'real',
                    value = 1.,
                    texname = '\\text{RMNS1x1}',
                    lhablock = 'UPMNS',
                    lhacode = [ 1, 1 ])

RMNS2x2 = Parameter(name = 'RMNS2x2',
                    nature = 'external',
                    type = 'real',
                    value = 1.,
                    texname = '\\text{RMNS2x2}',
                    lhablock = 'UPMNS',
                    lhacode = [ 2, 2 ])

RMNS3x3 = Parameter(name = 'RMNS3x3',
                    nature = 'external',
                    type = 'real',
                    value = 1.,
                    texname = '\\text{RMNS3x3}',
                    lhablock = 'UPMNS',
                    lhacode = [ 3, 3 ])

RRu1x1 = Parameter(name = 'RRu1x1',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRu1x1}',
                   lhablock = 'USQMIX',
                   lhacode = [ 1, 1 ])

RRu2x2 = Parameter(name = 'RRu2x2',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRu2x2}',
                   lhablock = 'USQMIX',
                   lhacode = [ 2, 2 ])

RRu3x3 = Parameter(name = 'RRu3x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.55364496,
                   texname = '\\text{RRu3x3}',
                   lhablock = 'USQMIX',
                   lhacode = [ 3, 3 ])

RRu3x6 = Parameter(name = 'RRu3x6',
                   nature = 'external',
                   type = 'real',
                   value = 0.83275282,
                   texname = '\\text{RRu3x6}',
                   lhablock = 'USQMIX',
                   lhacode = [ 3, 6 ])

RRu4x4 = Parameter(name = 'RRu4x4',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRu4x4}',
                   lhablock = 'USQMIX',
                   lhacode = [ 4, 4 ])

RRu5x5 = Parameter(name = 'RRu5x5',
                   nature = 'external',
                   type = 'real',
                   value = 1.,
                   texname = '\\text{RRu5x5}',
                   lhablock = 'USQMIX',
                   lhacode = [ 5, 5 ])

RRu6x3 = Parameter(name = 'RRu6x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.83275282,
                   texname = '\\text{RRu6x3}',
                   lhablock = 'USQMIX',
                   lhacode = [ 6, 3 ])

RRu6x6 = Parameter(name = 'RRu6x6',
                   nature = 'external',
                   type = 'real',
                   value = -0.55364496,
                   texname = '\\text{RRu6x6}',
                   lhablock = 'USQMIX',
                   lhacode = [ 6, 6 ])

RCKM1x1 = Parameter(name = 'RCKM1x1',
                    nature = 'external',
                    type = 'real',
                    value = 1.,
                    texname = '\\text{RCKM1x1}',
                    lhablock = 'VCKM',
                    lhacode = [ 1, 1 ])

RCKM2x2 = Parameter(name = 'RCKM2x2',
                    nature = 'external',
                    type = 'real',
                    value = 1.,
                    texname = '\\text{RCKM2x2}',
                    lhablock = 'VCKM',
                    lhacode = [ 2, 2 ])

RCKM3x3 = Parameter(name = 'RCKM3x3',
                    nature = 'external',
                    type = 'real',
                    value = 1.,
                    texname = '\\text{RCKM3x3}',
                    lhablock = 'VCKM',
                    lhacode = [ 3, 3 ])

RVV1x1 = Parameter(name = 'RVV1x1',
                   nature = 'external',
                   type = 'real',
                   value = 0.972557835,
                   texname = '\\text{RVV1x1}',
                   lhablock = 'VMIX',
                   lhacode = [ 1, 1 ])

RVV1x2 = Parameter(name = 'RVV1x2',
                   nature = 'external',
                   type = 'real',
                   value = -0.232661249,
                   texname = '\\text{RVV1x2}',
                   lhablock = 'VMIX',
                   lhacode = [ 1, 2 ])

RVV2x1 = Parameter(name = 'RVV2x1',
                   nature = 'external',
                   type = 'real',
                   value = 0.232661249,
                   texname = '\\text{RVV2x1}',
                   lhablock = 'VMIX',
                   lhacode = [ 2, 1 ])

RVV2x2 = Parameter(name = 'RVV2x2',
                   nature = 'external',
                   type = 'real',
                   value = 0.972557835,
                   texname = '\\text{RVV2x2}',
                   lhablock = 'VMIX',
                   lhacode = [ 2, 2 ])

Ryd3x3 = Parameter(name = 'Ryd3x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.138840206,
                   texname = '\\text{Ryd3x3}',
                   lhablock = 'YD',
                   lhacode = [ 3, 3 ])

Rye3x3 = Parameter(name = 'Rye3x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.10089081,
                   texname = '\\text{Rye3x3}',
                   lhablock = 'YE',
                   lhacode = [ 3, 3 ])

Ryu3x3 = Parameter(name = 'Ryu3x3',
                   nature = 'external',
                   type = 'real',
                   value = 0.89284455,
                   texname = '\\text{Ryu3x3}',
                   lhablock = 'YU',
                   lhacode = [ 3, 3 ])

MZ = Parameter(name = 'MZ',
               nature = 'external',
               type = 'real',
               value = 91.1876,
               texname = '\\text{MZ}',
               lhablock = 'MASS',
               lhacode = [ 23 ])

MW = Parameter(name = 'MW',
               nature = 'external',
               type = 'real',
               value = 79.8290131,
               texname = '\\text{MW}',
               lhablock = 'MASS',
               lhacode = [ 24 ])

Mneu1 = Parameter(name = 'Mneu1',
                  nature = 'external',
                  type = 'real',
                  value = 96.6880686,
                  texname = '\\text{Mneu1}',
                  lhablock = 'MASS',
                  lhacode = [ 1000022 ])

Mneu2 = Parameter(name = 'Mneu2',
                  nature = 'external',
                  type = 'real',
                  value = 181.088157,
                  texname = '\\text{Mneu2}',
                  lhablock = 'MASS',
                  lhacode = [ 1000023 ])

Mneu3 = Parameter(name = 'Mneu3',
                  nature = 'external',
                  type = 'real',
                  value = -363.756027,
                  texname = '\\text{Mneu3}',
                  lhablock = 'MASS',
                  lhacode = [ 1000025 ])

Mneu4 = Parameter(name = 'Mneu4',
                  nature = 'external',
                  type = 'real',
                  value = 381.729382,
                  texname = '\\text{Mneu4}',
                  lhablock = 'MASS',
                  lhacode = [ 1000035 ])

Mch1 = Parameter(name = 'Mch1',
                 nature = 'external',
                 type = 'real',
                 value = 181.696474,
                 texname = '\\text{Mch1}',
                 lhablock = 'MASS',
                 lhacode = [ 1000024 ])

Mch2 = Parameter(name = 'Mch2',
                 nature = 'external',
                 type = 'real',
                 value = 379.93932,
                 texname = '\\text{Mch2}',
                 lhablock = 'MASS',
                 lhacode = [ 1000037 ])

Mgo = Parameter(name = 'Mgo',
                nature = 'external',
                type = 'real',
                value = 607.713704,
                texname = '\\text{Mgo}',
                lhablock = 'MASS',
                lhacode = [ 1000021 ])

MH01 = Parameter(name = 'MH01',
                 nature = 'external',
                 type = 'real',
                 value = 110.899057,
                 texname = '\\text{MH01}',
                 lhablock = 'MASS',
                 lhacode = [ 25 ])

MH02 = Parameter(name = 'MH02',
                 nature = 'external',
                 type = 'real',
                 value = 399.960116,
                 texname = '\\text{MH02}',
                 lhablock = 'MASS',
                 lhacode = [ 35 ])

MA0 = Parameter(name = 'MA0',
                nature = 'external',
                type = 'real',
                value = 399.583917,
                texname = '\\text{MA0}',
                lhablock = 'MASS',
                lhacode = [ 36 ])

MH = Parameter(name = 'MH',
               nature = 'external',
               type = 'real',
               value = 407.879012,
               texname = '\\text{MH}',
               lhablock = 'MASS',
               lhacode = [ 37 ])

Mta = Parameter(name = 'Mta',
                nature = 'external',
                type = 'real',
                value = 1.777,
                texname = '\\text{Mta}',
                lhablock = 'MASS',
                lhacode = [ 15 ])

MT = Parameter(name = 'MT',
               nature = 'external',
               type = 'real',
               value = 175.,
               texname = '\\text{MT}',
               lhablock = 'MASS',
               lhacode = [ 6 ])

MB = Parameter(name = 'MB',
               nature = 'external',
               type = 'real',
               value = 4.88991651,
               texname = '\\text{MB}',
               lhablock = 'MASS',
               lhacode = [ 5 ])

Msn1 = Parameter(name = 'Msn1',
                 nature = 'external',
                 type = 'real',
                 value = 185.258326,
                 texname = '\\text{Msn1}',
                 lhablock = 'MASS',
                 lhacode = [ 1000012 ])

Msn2 = Parameter(name = 'Msn2',
                 nature = 'external',
                 type = 'real',
                 value = 185.258326,
                 texname = '\\text{Msn2}',
                 lhablock = 'MASS',
                 lhacode = [ 1000014 ])

Msn3 = Parameter(name = 'Msn3',
                 nature = 'external',
                 type = 'real',
                 value = 184.708464,
                 texname = '\\text{Msn3}',
                 lhablock = 'MASS',
                 lhacode = [ 1000016 ])

Msl1 = Parameter(name = 'Msl1',
                 nature = 'external',
                 type = 'real',
                 value = 202.91569,
                 texname = '\\text{Msl1}',
                 lhablock = 'MASS',
                 lhacode = [ 1000011 ])

Msl2 = Parameter(name = 'Msl2',
                 nature = 'external',
                 type = 'real',
                 value = 202.91569,
                 texname = '\\text{Msl2}',
                 lhablock = 'MASS',
                 lhacode = [ 1000013 ])

Msl3 = Parameter(name = 'Msl3',
                 nature = 'external',
                 type = 'real',
                 value = 134.490864,
                 texname = '\\text{Msl3}',
                 lhablock = 'MASS',
                 lhacode = [ 1000015 ])

Msl4 = Parameter(name = 'Msl4',
                 nature = 'external',
                 type = 'real',
                 value = 144.102799,
                 texname = '\\text{Msl4}',
                 lhablock = 'MASS',
                 lhacode = [ 2000011 ])

Msl5 = Parameter(name = 'Msl5',
                 nature = 'external',
                 type = 'real',
                 value = 144.102799,
                 texname = '\\text{Msl5}',
                 lhablock = 'MASS',
                 lhacode = [ 2000013 ])

Msl6 = Parameter(name = 'Msl6',
                 nature = 'external',
                 type = 'real',
                 value = 206.867805,
                 texname = '\\text{Msl6}',
                 lhablock = 'MASS',
                 lhacode = [ 2000015 ])

Msu1 = Parameter(name = 'Msu1',
                 nature = 'external',
                 type = 'real',
                 value = 561.119014,
                 texname = '\\text{Msu1}',
                 lhablock = 'MASS',
                 lhacode = [ 1000002 ])

Msu2 = Parameter(name = 'Msu2',
                 nature = 'external',
                 type = 'real',
                 value = 561.119014,
                 texname = '\\text{Msu2}',
                 lhablock = 'MASS',
                 lhacode = [ 1000004 ])

Msu3 = Parameter(name = 'Msu3',
                 nature = 'external',
                 type = 'real',
                 value = 399.668493,
                 texname = '\\text{Msu3}',
                 lhablock = 'MASS',
                 lhacode = [ 1000006 ])

Msu4 = Parameter(name = 'Msu4',
                 nature = 'external',
                 type = 'real',
                 value = 549.259265,
                 texname = '\\text{Msu4}',
                 lhablock = 'MASS',
                 lhacode = [ 2000002 ])

Msu5 = Parameter(name = 'Msu5',
                 nature = 'external',
                 type = 'real',
                 value = 549.259265,
                 texname = '\\text{Msu5}',
                 lhablock = 'MASS',
                 lhacode = [ 2000004 ])

Msu6 = Parameter(name = 'Msu6',
                 nature = 'external',
                 type = 'real',
                 value = 585.785818,
                 texname = '\\text{Msu6}',
                 lhablock = 'MASS',
                 lhacode = [ 2000006 ])

Msd1 = Parameter(name = 'Msd1',
                 nature = 'external',
                 type = 'real',
                 value = 568.441109,
                 texname = '\\text{Msd1}',
                 lhablock = 'MASS',
                 lhacode = [ 1000001 ])

Msd2 = Parameter(name = 'Msd2',
                 nature = 'external',
                 type = 'real',
                 value = 568.441109,
                 texname = '\\text{Msd2}',
                 lhablock = 'MASS',
                 lhacode = [ 1000003 ])

Msd3 = Parameter(name = 'Msd3',
                 nature = 'external',
                 type = 'real',
                 value = 513.065179,
                 texname = '\\text{Msd3}',
                 lhablock = 'MASS',
                 lhacode = [ 1000005 ])

Msd4 = Parameter(name = 'Msd4',
                 nature = 'external',
                 type = 'real',
                 value = 545.228462,
                 texname = '\\text{Msd4}',
                 lhablock = 'MASS',
                 lhacode = [ 2000001 ])

Msd5 = Parameter(name = 'Msd5',
                 nature = 'external',
                 type = 'real',
                 value = 545.228462,
                 texname = '\\text{Msd5}',
                 lhablock = 'MASS',
                 lhacode = [ 2000003 ])

Msd6 = Parameter(name = 'Msd6',
                 nature = 'external',
                 type = 'real',
                 value = 543.726676,
                 texname = '\\text{Msd6}',
                 lhablock = 'MASS',
                 lhacode = [ 2000005 ])

WZ = Parameter(name = 'WZ',
               nature = 'external',
               type = 'real',
               value = 2.41143316,
               texname = '\\text{WZ}',
               lhablock = 'DECAY',
               lhacode = [ 23 ])

WW = Parameter(name = 'WW',
               nature = 'external',
               type = 'real',
               value = 2.00282196,
               texname = '\\text{WW}',
               lhablock = 'DECAY',
               lhacode = [ 24 ])

Wneu2 = Parameter(name = 'Wneu2',
                  nature = 'external',
                  type = 'real',
                  value = 0.0207770048,
                  texname = '\\text{Wneu2}',
                  lhablock = 'DECAY',
                  lhacode = [ 1000023 ])

Wneu3 = Parameter(name = 'Wneu3',
                  nature = 'external',
                  type = 'real',
                  value = 1.91598495,
                  texname = '\\text{Wneu3}',
                  lhablock = 'DECAY',
                  lhacode = [ 1000025 ])

Wneu4 = Parameter(name = 'Wneu4',
                  nature = 'external',
                  type = 'real',
                  value = 2.58585079,
                  texname = '\\text{Wneu4}',
                  lhablock = 'DECAY',
                  lhacode = [ 1000035 ])

Wch1 = Parameter(name = 'Wch1',
                 nature = 'external',
                 type = 'real',
                 value = 0.0170414503,
                 texname = '\\text{Wch1}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000024 ])

Wch2 = Parameter(name = 'Wch2',
                 nature = 'external',
                 type = 'real',
                 value = 2.4868951,
                 texname = '\\text{Wch2}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000037 ])

Wgo = Parameter(name = 'Wgo',
                nature = 'external',
                type = 'real',
                value = 5.50675438,
                texname = '\\text{Wgo}',
                lhablock = 'DECAY',
                lhacode = [ 1000021 ])

WH01 = Parameter(name = 'WH01',
                 nature = 'external',
                 type = 'real',
                 value = 0.00198610799,
                 texname = '\\text{WH01}',
                 lhablock = 'DECAY',
                 lhacode = [ 25 ])

WH02 = Parameter(name = 'WH02',
                 nature = 'external',
                 type = 'real',
                 value = 0.574801389,
                 texname = '\\text{WH02}',
                 lhablock = 'DECAY',
                 lhacode = [ 35 ])

WA0 = Parameter(name = 'WA0',
                nature = 'external',
                type = 'real',
                value = 0.632178488,
                texname = '\\text{WA0}',
                lhablock = 'DECAY',
                lhacode = [ 36 ])

WH = Parameter(name = 'WH',
               nature = 'external',
               type = 'real',
               value = 0.546962813,
               texname = '\\text{WH}',
               lhablock = 'DECAY',
               lhacode = [ 37 ])

WT = Parameter(name = 'WT',
               nature = 'external',
               type = 'real',
               value = 1.56194983,
               texname = '\\text{WT}',
               lhablock = 'DECAY',
               lhacode = [ 6 ])

Wsn1 = Parameter(name = 'Wsn1',
                 nature = 'external',
                 type = 'real',
                 value = 0.149881634,
                 texname = '\\text{Wsn1}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000012 ])

Wsn2 = Parameter(name = 'Wsn2',
                 nature = 'external',
                 type = 'real',
                 value = 0.149881634,
                 texname = '\\text{Wsn2}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000014 ])

Wsn3 = Parameter(name = 'Wsn3',
                 nature = 'external',
                 type = 'real',
                 value = 0.147518977,
                 texname = '\\text{Wsn3}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000016 ])

Wsl1 = Parameter(name = 'Wsl1',
                 nature = 'external',
                 type = 'real',
                 value = 0.213682161,
                 texname = '\\text{Wsl1}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000011 ])

Wsl2 = Parameter(name = 'Wsl2',
                 nature = 'external',
                 type = 'real',
                 value = 0.213682161,
                 texname = '\\text{Wsl2}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000013 ])

Wsl3 = Parameter(name = 'Wsl3',
                 nature = 'external',
                 type = 'real',
                 value = 0.148327268,
                 texname = '\\text{Wsl3}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000015 ])

Wsl4 = Parameter(name = 'Wsl4',
                 nature = 'external',
                 type = 'real',
                 value = 0.216121626,
                 texname = '\\text{Wsl4}',
                 lhablock = 'DECAY',
                 lhacode = [ 2000011 ])

Wsl5 = Parameter(name = 'Wsl5',
                 nature = 'external',
                 type = 'real',
                 value = 0.216121626,
                 texname = '\\text{Wsl5}',
                 lhablock = 'DECAY',
                 lhacode = [ 2000013 ])

Wsl6 = Parameter(name = 'Wsl6',
                 nature = 'external',
                 type = 'real',
                 value = 0.269906096,
                 texname = '\\text{Wsl6}',
                 lhablock = 'DECAY',
                 lhacode = [ 2000015 ])

Wsu1 = Parameter(name = 'Wsu1',
                 nature = 'external',
                 type = 'real',
                 value = 5.47719539,
                 texname = '\\text{Wsu1}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000002 ])

Wsu2 = Parameter(name = 'Wsu2',
                 nature = 'external',
                 type = 'real',
                 value = 5.47719539,
                 texname = '\\text{Wsu2}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000004 ])

Wsu3 = Parameter(name = 'Wsu3',
                 nature = 'external',
                 type = 'real',
                 value = 2.02159578,
                 texname = '\\text{Wsu3}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000006 ])

Wsu4 = Parameter(name = 'Wsu4',
                 nature = 'external',
                 type = 'real',
                 value = 1.15297292,
                 texname = '\\text{Wsu4}',
                 lhablock = 'DECAY',
                 lhacode = [ 2000002 ])

Wsu5 = Parameter(name = 'Wsu5',
                 nature = 'external',
                 type = 'real',
                 value = 1.15297292,
                 texname = '\\text{Wsu5}',
                 lhablock = 'DECAY',
                 lhacode = [ 2000004 ])

Wsu6 = Parameter(name = 'Wsu6',
                 nature = 'external',
                 type = 'real',
                 value = 7.37313275,
                 texname = '\\text{Wsu6}',
                 lhablock = 'DECAY',
                 lhacode = [ 2000006 ])

Wsd1 = Parameter(name = 'Wsd1',
                 nature = 'external',
                 type = 'real',
                 value = 5.31278772,
                 texname = '\\text{Wsd1}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000001 ])

Wsd2 = Parameter(name = 'Wsd2',
                 nature = 'external',
                 type = 'real',
                 value = 5.31278772,
                 texname = '\\text{Wsd2}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000003 ])

Wsd3 = Parameter(name = 'Wsd3',
                 nature = 'external',
                 type = 'real',
                 value = 3.73627601,
                 texname = '\\text{Wsd3}',
                 lhablock = 'DECAY',
                 lhacode = [ 1000005 ])

Wsd4 = Parameter(name = 'Wsd4',
                 nature = 'external',
                 type = 'real',
                 value = 0.285812308,
                 texname = '\\text{Wsd4}',
                 lhablock = 'DECAY',
                 lhacode = [ 2000001 ])

Wsd5 = Parameter(name = 'Wsd5',
                 nature = 'external',
                 type = 'real',
                 value = 0.285812308,
                 texname = '\\text{Wsd5}',
                 lhablock = 'DECAY',
                 lhacode = [ 2000003 ])

Wsd6 = Parameter(name = 'Wsd6',
                 nature = 'external',
                 type = 'real',
                 value = 0.801566294,
                 texname = '\\text{Wsd6}',
                 lhablock = 'DECAY',
                 lhacode = [ 2000005 ])

beta = Parameter(name = 'beta',
                 nature = 'internal',
                 type = 'real',
                 value = 'cmath.atan(tb)',
                 texname = '\\beta')

CKM1x1 = Parameter(name = 'CKM1x1',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RCKM1x1',
                   texname = '\\text{CKM1x1}')

CKM2x2 = Parameter(name = 'CKM2x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RCKM2x2',
                   texname = '\\text{CKM2x2}')

CKM3x3 = Parameter(name = 'CKM3x3',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RCKM3x3',
                   texname = '\\text{CKM3x3}')

cw = Parameter(name = 'cw',
               nature = 'internal',
               type = 'real',
               value = 'MW/MZ',
               texname = 'c_w')

mD21x1 = Parameter(name = 'mD21x1',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmD21x1',
                   texname = '\\text{mD21x1}')

mD22x2 = Parameter(name = 'mD22x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmD22x2',
                   texname = '\\text{mD22x2}')

mD23x3 = Parameter(name = 'mD23x3',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmD23x3',
                   texname = '\\text{mD23x3}')

mE21x1 = Parameter(name = 'mE21x1',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmE21x1',
                   texname = '\\text{mE21x1}')

mE22x2 = Parameter(name = 'mE22x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmE22x2',
                   texname = '\\text{mE22x2}')

mE23x3 = Parameter(name = 'mE23x3',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmE23x3',
                   texname = '\\text{mE23x3}')

mL21x1 = Parameter(name = 'mL21x1',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmL21x1',
                   texname = '\\text{mL21x1}')

mL22x2 = Parameter(name = 'mL22x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmL22x2',
                   texname = '\\text{mL22x2}')

mL23x3 = Parameter(name = 'mL23x3',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmL23x3',
                   texname = '\\text{mL23x3}')

mQ21x1 = Parameter(name = 'mQ21x1',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmQ21x1',
                   texname = '\\text{mQ21x1}')

mQ22x2 = Parameter(name = 'mQ22x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmQ22x2',
                   texname = '\\text{mQ22x2}')

mQ23x3 = Parameter(name = 'mQ23x3',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmQ23x3',
                   texname = '\\text{mQ23x3}')

mU21x1 = Parameter(name = 'mU21x1',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmU21x1',
                   texname = '\\text{mU21x1}')

mU22x2 = Parameter(name = 'mU22x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmU22x2',
                   texname = '\\text{mU22x2}')

mU23x3 = Parameter(name = 'mU23x3',
                   nature = 'internal',
                   type = 'complex',
                   value = 'RmU23x3',
                   texname = '\\text{mU23x3}')

MUH = Parameter(name = 'MUH',
                nature = 'internal',
                type = 'complex',
                value = 'RMUH',
                texname = '\\mu')

Mx1 = Parameter(name = 'Mx1',
                nature = 'internal',
                type = 'complex',
                value = 'RMx1',
                texname = 'M_1')

Mx2 = Parameter(name = 'Mx2',
                nature = 'internal',
                type = 'complex',
                value = 'RMx2',
                texname = 'M_2')

Mx3 = Parameter(name = 'Mx3',
                nature = 'internal',
                type = 'complex',
                value = 'RMx3',
                texname = 'M_3')

NN1x1 = Parameter(name = 'NN1x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN1x1',
                  texname = '\\text{NN1x1}')

NN1x2 = Parameter(name = 'NN1x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN1x2',
                  texname = '\\text{NN1x2}')

NN1x3 = Parameter(name = 'NN1x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN1x3',
                  texname = '\\text{NN1x3}')

NN1x4 = Parameter(name = 'NN1x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN1x4',
                  texname = '\\text{NN1x4}')

NN2x1 = Parameter(name = 'NN2x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN2x1',
                  texname = '\\text{NN2x1}')

NN2x2 = Parameter(name = 'NN2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN2x2',
                  texname = '\\text{NN2x2}')

NN2x3 = Parameter(name = 'NN2x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN2x3',
                  texname = '\\text{NN2x3}')

NN2x4 = Parameter(name = 'NN2x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN2x4',
                  texname = '\\text{NN2x4}')

NN3x1 = Parameter(name = 'NN3x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN3x1',
                  texname = '\\text{NN3x1}')

NN3x2 = Parameter(name = 'NN3x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN3x2',
                  texname = '\\text{NN3x2}')

NN3x3 = Parameter(name = 'NN3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN3x3',
                  texname = '\\text{NN3x3}')

NN3x4 = Parameter(name = 'NN3x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN3x4',
                  texname = '\\text{NN3x4}')

NN4x1 = Parameter(name = 'NN4x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN4x1',
                  texname = '\\text{NN4x1}')

NN4x2 = Parameter(name = 'NN4x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN4x2',
                  texname = '\\text{NN4x2}')

NN4x3 = Parameter(name = 'NN4x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN4x3',
                  texname = '\\text{NN4x3}')

NN4x4 = Parameter(name = 'NN4x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RNN4x4',
                  texname = '\\text{NN4x4}')

Rd1x1 = Parameter(name = 'Rd1x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd1x1',
                  texname = '\\text{Rd1x1}')

Rd2x2 = Parameter(name = 'Rd2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd2x2',
                  texname = '\\text{Rd2x2}')

Rd3x3 = Parameter(name = 'Rd3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd3x3',
                  texname = '\\text{Rd3x3}')

Rd3x6 = Parameter(name = 'Rd3x6',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd3x6',
                  texname = '\\text{Rd3x6}')

Rd4x4 = Parameter(name = 'Rd4x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd4x4',
                  texname = '\\text{Rd4x4}')

Rd5x5 = Parameter(name = 'Rd5x5',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd5x5',
                  texname = '\\text{Rd5x5}')

Rd6x3 = Parameter(name = 'Rd6x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd6x3',
                  texname = '\\text{Rd6x3}')

Rd6x6 = Parameter(name = 'Rd6x6',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRd6x6',
                  texname = '\\text{Rd6x6}')

Rl1x1 = Parameter(name = 'Rl1x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl1x1',
                  texname = '\\text{Rl1x1}')

Rl2x2 = Parameter(name = 'Rl2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl2x2',
                  texname = '\\text{Rl2x2}')

Rl3x3 = Parameter(name = 'Rl3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl3x3',
                  texname = '\\text{Rl3x3}')

Rl3x6 = Parameter(name = 'Rl3x6',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl3x6',
                  texname = '\\text{Rl3x6}')

Rl4x4 = Parameter(name = 'Rl4x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl4x4',
                  texname = '\\text{Rl4x4}')

Rl5x5 = Parameter(name = 'Rl5x5',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl5x5',
                  texname = '\\text{Rl5x5}')

Rl6x3 = Parameter(name = 'Rl6x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl6x3',
                  texname = '\\text{Rl6x3}')

Rl6x6 = Parameter(name = 'Rl6x6',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRl6x6',
                  texname = '\\text{Rl6x6}')

Rn1x1 = Parameter(name = 'Rn1x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRn1x1',
                  texname = '\\text{Rn1x1}')

Rn2x2 = Parameter(name = 'Rn2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRn2x2',
                  texname = '\\text{Rn2x2}')

Rn3x3 = Parameter(name = 'Rn3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRn3x3',
                  texname = '\\text{Rn3x3}')

Ru1x1 = Parameter(name = 'Ru1x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu1x1',
                  texname = '\\text{Ru1x1}')

Ru2x2 = Parameter(name = 'Ru2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu2x2',
                  texname = '\\text{Ru2x2}')

Ru3x3 = Parameter(name = 'Ru3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu3x3',
                  texname = '\\text{Ru3x3}')

Ru3x6 = Parameter(name = 'Ru3x6',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu3x6',
                  texname = '\\text{Ru3x6}')

Ru4x4 = Parameter(name = 'Ru4x4',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu4x4',
                  texname = '\\text{Ru4x4}')

Ru5x5 = Parameter(name = 'Ru5x5',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu5x5',
                  texname = '\\text{Ru5x5}')

Ru6x3 = Parameter(name = 'Ru6x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu6x3',
                  texname = '\\text{Ru6x3}')

Ru6x6 = Parameter(name = 'Ru6x6',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RRu6x6',
                  texname = '\\text{Ru6x6}')

UU1x1 = Parameter(name = 'UU1x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RUU1x1',
                  texname = '\\text{UU1x1}')

UU1x2 = Parameter(name = 'UU1x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RUU1x2',
                  texname = '\\text{UU1x2}')

UU2x1 = Parameter(name = 'UU2x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RUU2x1',
                  texname = '\\text{UU2x1}')

UU2x2 = Parameter(name = 'UU2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RUU2x2',
                  texname = '\\text{UU2x2}')

VV1x1 = Parameter(name = 'VV1x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RVV1x1',
                  texname = '\\text{VV1x1}')

VV1x2 = Parameter(name = 'VV1x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RVV1x2',
                  texname = '\\text{VV1x2}')

VV2x1 = Parameter(name = 'VV2x1',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RVV2x1',
                  texname = '\\text{VV2x1}')

VV2x2 = Parameter(name = 'VV2x2',
                  nature = 'internal',
                  type = 'complex',
                  value = 'RVV2x2',
                  texname = '\\text{VV2x2}')

ee = Parameter(name = 'ee',
               nature = 'internal',
               type = 'real',
               value = '2*cmath.sqrt(1/aEWM1)*cmath.sqrt(cmath.pi)',
               texname = 'e')

G = Parameter(name = 'G',
              nature = 'internal',
              type = 'real',
              value = '2*cmath.sqrt(aS)*cmath.sqrt(cmath.pi)',
              texname = 'G')

td3x3 = Parameter(name = 'td3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rtd3x3',
                  texname = '\\text{td3x3}')

te3x3 = Parameter(name = 'te3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rte3x3',
                  texname = '\\text{te3x3}')

tu3x3 = Parameter(name = 'tu3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rtu3x3',
                  texname = '\\text{tu3x3}')

yd3x3 = Parameter(name = 'yd3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Ryd3x3',
                  texname = '\\text{yd3x3}')

ye3x3 = Parameter(name = 'ye3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rye3x3',
                  texname = '\\text{ye3x3}')

yu3x3 = Parameter(name = 'yu3x3',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Ryu3x3',
                  texname = '\\text{yu3x3}')

bb = Parameter(name = 'bb',
               nature = 'internal',
               type = 'complex',
               value = '((-mHd2 + mHu2)*cmath.tan(2*alp))/2. - MZ**2*(cmath.sin(2*beta)/2. + cmath.cos(2*beta)*cmath.tan(2*alp))',
               texname = 'b')

sw = Parameter(name = 'sw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(1 - cw**2)',
               texname = 's_w')

gp = Parameter(name = 'gp',
               nature = 'internal',
               type = 'real',
               value = 'ee/cw',
               texname = 'g\'')

gw = Parameter(name = 'gw',
               nature = 'internal',
               type = 'real',
               value = 'ee/sw',
               texname = 'g_w')

vev = Parameter(name = 'vev',
                nature = 'internal',
                type = 'real',
                value = '(2*cw*MZ*sw)/ee',
                texname = 'v')

vd = Parameter(name = 'vd',
               nature = 'internal',
               type = 'real',
               value = 'vev*cmath.cos(beta)',
               texname = 'v_d')

vu = Parameter(name = 'vu',
               nature = 'internal',
               type = 'real',
               value = 'vev*cmath.sin(beta)',
               texname = 'v_u')

I1a33 = Parameter(name = 'I1a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complexconjugate(CKM3x3)*complexconjugate(yu3x3)',
                  texname = '\\text{I1a33}')

I10a11 = Parameter(name = 'I10a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(CKM1x1)',
                   texname = '\\text{I10a11}')

I10a22 = Parameter(name = 'I10a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(CKM2x2)',
                   texname = '\\text{I10a22}')

I10a33 = Parameter(name = 'I10a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(CKM3x3)',
                   texname = '\\text{I10a33}')

I10a36 = Parameter(name = 'I10a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(CKM3x3)',
                   texname = '\\text{I10a36}')

I100a33 = Parameter(name = 'I100a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x6*complexconjugate(Rd3x6)',
                    texname = '\\text{I100a33}')

I100a36 = Parameter(name = 'I100a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x6*complexconjugate(Rd3x6)',
                    texname = '\\text{I100a36}')

I100a44 = Parameter(name = 'I100a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd4x4*complexconjugate(Rd4x4)',
                    texname = '\\text{I100a44}')

I100a55 = Parameter(name = 'I100a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd5x5*complexconjugate(Rd5x5)',
                    texname = '\\text{I100a55}')

I100a63 = Parameter(name = 'I100a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd3x6*complexconjugate(Rd6x6)',
                    texname = '\\text{I100a63}')

I100a66 = Parameter(name = 'I100a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rd6x6*complexconjugate(Rd6x6)',
                    texname = '\\text{I100a66}')

I101a33 = Parameter(name = 'I101a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rl3x6*complexconjugate(Rl3x6)',
                    texname = '\\text{I101a33}')

I101a36 = Parameter(name = 'I101a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rl6x6*complexconjugate(Rl3x6)',
                    texname = '\\text{I101a36}')

I101a44 = Parameter(name = 'I101a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rl4x4*complexconjugate(Rl4x4)',
                    texname = '\\text{I101a44}')

I101a55 = Parameter(name = 'I101a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rl5x5*complexconjugate(Rl5x5)',
                    texname = '\\text{I101a55}')

I101a63 = Parameter(name = 'I101a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rl3x6*complexconjugate(Rl6x6)',
                    texname = '\\text{I101a63}')

I101a66 = Parameter(name = 'I101a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Rl6x6*complexconjugate(Rl6x6)',
                    texname = '\\text{I101a66}')

I102a33 = Parameter(name = 'I102a33',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I102a33}')

I102a36 = Parameter(name = 'I102a36',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru3x6)',
                    texname = '\\text{I102a36}')

I102a44 = Parameter(name = 'I102a44',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru4x4*complexconjugate(Ru4x4)',
                    texname = '\\text{I102a44}')

I102a55 = Parameter(name = 'I102a55',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru5x5*complexconjugate(Ru5x5)',
                    texname = '\\text{I102a55}')

I102a63 = Parameter(name = 'I102a63',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru3x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I102a63}')

I102a66 = Parameter(name = 'I102a66',
                    nature = 'internal',
                    type = 'complex',
                    value = 'Ru6x6*complexconjugate(Ru6x6)',
                    texname = '\\text{I102a66}')

I11a33 = Parameter(name = 'I11a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(yu3x3)',
                   texname = '\\text{I11a33}')

I11a36 = Parameter(name = 'I11a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(yu3x3)',
                   texname = '\\text{I11a36}')

I12a33 = Parameter(name = 'I12a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(CKM3x3)',
                   texname = '\\text{I12a33}')

I12a36 = Parameter(name = 'I12a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(CKM3x3)',
                   texname = '\\text{I12a36}')

I13a33 = Parameter(name = 'I13a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(yd3x3)',
                   texname = '\\text{I13a33}')

I13a36 = Parameter(name = 'I13a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(yd3x3)',
                   texname = '\\text{I13a36}')

I14a33 = Parameter(name = 'I14a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3',
                   texname = '\\text{I14a33}')

I14a36 = Parameter(name = 'I14a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3',
                   texname = '\\text{I14a36}')

I15a11 = Parameter(name = 'I15a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I15a11}')

I15a22 = Parameter(name = 'I15a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I15a22}')

I15a33 = Parameter(name = 'I15a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I15a33}')

I15a36 = Parameter(name = 'I15a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I15a36}')

I15a63 = Parameter(name = 'I15a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I15a63}')

I15a66 = Parameter(name = 'I15a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I15a66}')

I16a33 = Parameter(name = 'I16a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*complexconjugate(Rd3x6)',
                   texname = '\\text{I16a33}')

I16a36 = Parameter(name = 'I16a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*complexconjugate(Rd3x6)',
                   texname = '\\text{I16a36}')

I16a44 = Parameter(name = 'I16a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd4x4*complexconjugate(Rd4x4)',
                   texname = '\\text{I16a44}')

I16a55 = Parameter(name = 'I16a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd5x5*complexconjugate(Rd5x5)',
                   texname = '\\text{I16a55}')

I16a63 = Parameter(name = 'I16a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*complexconjugate(Rd6x6)',
                   texname = '\\text{I16a63}')

I16a66 = Parameter(name = 'I16a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*complexconjugate(Rd6x6)',
                   texname = '\\text{I16a66}')

I17a33 = Parameter(name = 'I17a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd3x6)*complexconjugate(td3x3)',
                   texname = '\\text{I17a33}')

I17a36 = Parameter(name = 'I17a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd3x6)*complexconjugate(td3x3)',
                   texname = '\\text{I17a36}')

I17a63 = Parameter(name = 'I17a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd6x6)*complexconjugate(td3x3)',
                   texname = '\\text{I17a63}')

I17a66 = Parameter(name = 'I17a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd6x6)*complexconjugate(td3x3)',
                   texname = '\\text{I17a66}')

I18a33 = Parameter(name = 'I18a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I18a33}')

I18a36 = Parameter(name = 'I18a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I18a36}')

I18a63 = Parameter(name = 'I18a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I18a63}')

I18a66 = Parameter(name = 'I18a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I18a66}')

I19a33 = Parameter(name = 'I19a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*td3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I19a33}')

I19a36 = Parameter(name = 'I19a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*td3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I19a36}')

I19a63 = Parameter(name = 'I19a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*td3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I19a63}')

I19a66 = Parameter(name = 'I19a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*td3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I19a66}')

I2a33 = Parameter(name = 'I2a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yd3x3*complexconjugate(CKM3x3)',
                  texname = '\\text{I2a33}')

I20a33 = Parameter(name = 'I20a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*yd3x3*complexconjugate(Rd3x3)*complexconjugate(yd3x3)',
                   texname = '\\text{I20a33}')

I20a36 = Parameter(name = 'I20a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*yd3x3*complexconjugate(Rd3x3)*complexconjugate(yd3x3)',
                   texname = '\\text{I20a36}')

I20a63 = Parameter(name = 'I20a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*yd3x3*complexconjugate(Rd6x3)*complexconjugate(yd3x3)',
                   texname = '\\text{I20a63}')

I20a66 = Parameter(name = 'I20a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*yd3x3*complexconjugate(Rd6x3)*complexconjugate(yd3x3)',
                   texname = '\\text{I20a66}')

I21a33 = Parameter(name = 'I21a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I21a33}')

I21a36 = Parameter(name = 'I21a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I21a36}')

I21a63 = Parameter(name = 'I21a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I21a63}')

I21a66 = Parameter(name = 'I21a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I21a66}')

I22a33 = Parameter(name = 'I22a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I22a33}')

I22a36 = Parameter(name = 'I22a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I22a36}')

I22a63 = Parameter(name = 'I22a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I22a63}')

I22a66 = Parameter(name = 'I22a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I22a66}')

I23a33 = Parameter(name = 'I23a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I23a33}')

I23a36 = Parameter(name = 'I23a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I23a36}')

I24a33 = Parameter(name = 'I24a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'ye3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I24a33}')

I24a36 = Parameter(name = 'I24a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'ye3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I24a36}')

I25a11 = Parameter(name = 'I25a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I25a11}')

I25a22 = Parameter(name = 'I25a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I25a22}')

I25a33 = Parameter(name = 'I25a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I25a33}')

I25a36 = Parameter(name = 'I25a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I25a36}')

I25a63 = Parameter(name = 'I25a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I25a63}')

I25a66 = Parameter(name = 'I25a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I25a66}')

I26a33 = Parameter(name = 'I26a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I26a33}')

I26a36 = Parameter(name = 'I26a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I26a36}')

I26a44 = Parameter(name = 'I26a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I26a44}')

I26a55 = Parameter(name = 'I26a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I26a55}')

I26a63 = Parameter(name = 'I26a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I26a63}')

I26a66 = Parameter(name = 'I26a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I26a66}')

I27a33 = Parameter(name = 'I27a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(ye3x3)',
                   texname = '\\text{I27a33}')

I27a36 = Parameter(name = 'I27a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(ye3x3)',
                   texname = '\\text{I27a36}')

I28a33 = Parameter(name = 'I28a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3',
                   texname = '\\text{I28a33}')

I28a36 = Parameter(name = 'I28a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3',
                   texname = '\\text{I28a36}')

I29a11 = Parameter(name = 'I29a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1',
                   texname = '\\text{I29a11}')

I29a22 = Parameter(name = 'I29a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2',
                   texname = '\\text{I29a22}')

I29a33 = Parameter(name = 'I29a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3',
                   texname = '\\text{I29a33}')

I29a36 = Parameter(name = 'I29a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3',
                   texname = '\\text{I29a36}')

I3a33 = Parameter(name = 'I3a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM3x3*complexconjugate(yd3x3)',
                  texname = '\\text{I3a33}')

I30a33 = Parameter(name = 'I30a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3',
                   texname = '\\text{I30a33}')

I30a36 = Parameter(name = 'I30a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3',
                   texname = '\\text{I30a36}')

I31a11 = Parameter(name = 'I31a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I31a11}')

I31a22 = Parameter(name = 'I31a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I31a22}')

I31a33 = Parameter(name = 'I31a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I31a33}')

I31a36 = Parameter(name = 'I31a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I31a36}')

I31a63 = Parameter(name = 'I31a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I31a63}')

I31a66 = Parameter(name = 'I31a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I31a66}')

I32a33 = Parameter(name = 'I32a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I32a33}')

I32a36 = Parameter(name = 'I32a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl3x6)',
                   texname = '\\text{I32a36}')

I32a44 = Parameter(name = 'I32a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl4x4*complexconjugate(Rl4x4)',
                   texname = '\\text{I32a44}')

I32a55 = Parameter(name = 'I32a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl5x5*complexconjugate(Rl5x5)',
                   texname = '\\text{I32a55}')

I32a63 = Parameter(name = 'I32a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I32a63}')

I32a66 = Parameter(name = 'I32a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*complexconjugate(Rl6x6)',
                   texname = '\\text{I32a66}')

I33a33 = Parameter(name = 'I33a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl3x6)*complexconjugate(te3x3)',
                   texname = '\\text{I33a33}')

I33a36 = Parameter(name = 'I33a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl3x6)*complexconjugate(te3x3)',
                   texname = '\\text{I33a36}')

I33a63 = Parameter(name = 'I33a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl6x6)*complexconjugate(te3x3)',
                   texname = '\\text{I33a63}')

I33a66 = Parameter(name = 'I33a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl6x6)*complexconjugate(te3x3)',
                   texname = '\\text{I33a66}')

I34a33 = Parameter(name = 'I34a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I34a33}')

I34a36 = Parameter(name = 'I34a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I34a36}')

I34a63 = Parameter(name = 'I34a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I34a63}')

I34a66 = Parameter(name = 'I34a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I34a66}')

I35a33 = Parameter(name = 'I35a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*te3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I35a33}')

I35a36 = Parameter(name = 'I35a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*te3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I35a36}')

I35a63 = Parameter(name = 'I35a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*te3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I35a63}')

I35a66 = Parameter(name = 'I35a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*te3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I35a66}')

I36a33 = Parameter(name = 'I36a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*ye3x3*complexconjugate(Rl3x3)*complexconjugate(ye3x3)',
                   texname = '\\text{I36a33}')

I36a36 = Parameter(name = 'I36a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*ye3x3*complexconjugate(Rl3x3)*complexconjugate(ye3x3)',
                   texname = '\\text{I36a36}')

I36a63 = Parameter(name = 'I36a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*ye3x3*complexconjugate(Rl6x3)*complexconjugate(ye3x3)',
                   texname = '\\text{I36a63}')

I36a66 = Parameter(name = 'I36a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*ye3x3*complexconjugate(Rl6x3)*complexconjugate(ye3x3)',
                   texname = '\\text{I36a66}')

I37a33 = Parameter(name = 'I37a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I37a33}')

I37a36 = Parameter(name = 'I37a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I37a36}')

I37a63 = Parameter(name = 'I37a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I37a63}')

I37a66 = Parameter(name = 'I37a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I37a66}')

I38a33 = Parameter(name = 'I38a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I38a33}')

I38a36 = Parameter(name = 'I38a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I38a36}')

I38a63 = Parameter(name = 'I38a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I38a63}')

I38a66 = Parameter(name = 'I38a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I38a66}')

I39a11 = Parameter(name = 'I39a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rn1x1)',
                   texname = '\\text{I39a11}')

I39a22 = Parameter(name = 'I39a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rn2x2)',
                   texname = '\\text{I39a22}')

I39a33 = Parameter(name = 'I39a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rn3x3)',
                   texname = '\\text{I39a33}')

I39a36 = Parameter(name = 'I39a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rn3x3)',
                   texname = '\\text{I39a36}')

I4a33 = Parameter(name = 'I4a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'CKM3x3*yu3x3',
                  texname = '\\text{I4a33}')

I40a33 = Parameter(name = 'I40a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*te3x3*complexconjugate(Rn3x3)',
                   texname = '\\text{I40a33}')

I40a36 = Parameter(name = 'I40a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*te3x3*complexconjugate(Rn3x3)',
                   texname = '\\text{I40a36}')

I41a33 = Parameter(name = 'I41a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*ye3x3*complexconjugate(Rn3x3)*complexconjugate(ye3x3)',
                   texname = '\\text{I41a33}')

I41a36 = Parameter(name = 'I41a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*ye3x3*complexconjugate(Rn3x3)*complexconjugate(ye3x3)',
                   texname = '\\text{I41a36}')

I42a33 = Parameter(name = 'I42a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x6*ye3x3*complexconjugate(Rn3x3)',
                   texname = '\\text{I42a33}')

I42a36 = Parameter(name = 'I42a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x6*ye3x3*complexconjugate(Rn3x3)',
                   texname = '\\text{I42a36}')

I43a11 = Parameter(name = 'I43a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rn1x1',
                   texname = '\\text{I43a11}')

I43a22 = Parameter(name = 'I43a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rn2x2',
                   texname = '\\text{I43a22}')

I43a33 = Parameter(name = 'I43a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rn3x3',
                   texname = '\\text{I43a33}')

I44a33 = Parameter(name = 'I44a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rn3x3*complexconjugate(ye3x3)',
                   texname = '\\text{I44a33}')

I45a11 = Parameter(name = 'I45a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rn1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I45a11}')

I45a22 = Parameter(name = 'I45a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rn2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I45a22}')

I45a33 = Parameter(name = 'I45a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rn3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I45a33}')

I45a36 = Parameter(name = 'I45a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rn3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I45a36}')

I46a33 = Parameter(name = 'I46a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rn3x3*complexconjugate(Rl3x6)*complexconjugate(te3x3)',
                   texname = '\\text{I46a33}')

I46a36 = Parameter(name = 'I46a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rn3x3*complexconjugate(Rl6x6)*complexconjugate(te3x3)',
                   texname = '\\text{I46a36}')

I47a33 = Parameter(name = 'I47a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rn3x3*complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I47a33}')

I47a36 = Parameter(name = 'I47a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rn3x3*complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I47a36}')

I48a33 = Parameter(name = 'I48a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rn3x3*ye3x3*complexconjugate(Rl3x3)*complexconjugate(ye3x3)',
                   texname = '\\text{I48a33}')

I48a36 = Parameter(name = 'I48a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rn3x3*ye3x3*complexconjugate(Rl6x3)*complexconjugate(ye3x3)',
                   texname = '\\text{I48a36}')

I49a33 = Parameter(name = 'I49a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                   texname = '\\text{I49a33}')

I49a36 = Parameter(name = 'I49a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                   texname = '\\text{I49a36}')

I5a33 = Parameter(name = 'I5a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complexconjugate(ye3x3)',
                  texname = '\\text{I5a33}')

I50a33 = Parameter(name = 'I50a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'yu3x3*complexconjugate(Ru3x3)',
                   texname = '\\text{I50a33}')

I50a36 = Parameter(name = 'I50a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'yu3x3*complexconjugate(Ru6x3)',
                   texname = '\\text{I50a36}')

I51a11 = Parameter(name = 'I51a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru1x1*complexconjugate(Ru1x1)',
                   texname = '\\text{I51a11}')

I51a22 = Parameter(name = 'I51a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru2x2*complexconjugate(Ru2x2)',
                   texname = '\\text{I51a22}')

I51a33 = Parameter(name = 'I51a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x3*complexconjugate(Ru3x3)',
                   texname = '\\text{I51a33}')

I51a36 = Parameter(name = 'I51a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x3*complexconjugate(Ru3x3)',
                   texname = '\\text{I51a36}')

I51a63 = Parameter(name = 'I51a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x3*complexconjugate(Ru6x3)',
                   texname = '\\text{I51a63}')

I51a66 = Parameter(name = 'I51a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x3*complexconjugate(Ru6x3)',
                   texname = '\\text{I51a66}')

I52a33 = Parameter(name = 'I52a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x6*complexconjugate(Ru3x6)',
                   texname = '\\text{I52a33}')

I52a36 = Parameter(name = 'I52a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x6*complexconjugate(Ru3x6)',
                   texname = '\\text{I52a36}')

I52a44 = Parameter(name = 'I52a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru4x4*complexconjugate(Ru4x4)',
                   texname = '\\text{I52a44}')

I52a55 = Parameter(name = 'I52a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru5x5*complexconjugate(Ru5x5)',
                   texname = '\\text{I52a55}')

I52a63 = Parameter(name = 'I52a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x6*complexconjugate(Ru6x6)',
                   texname = '\\text{I52a63}')

I52a66 = Parameter(name = 'I52a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x6*complexconjugate(Ru6x6)',
                   texname = '\\text{I52a66}')

I53a11 = Parameter(name = 'I53a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                   texname = '\\text{I53a11}')

I53a22 = Parameter(name = 'I53a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                   texname = '\\text{I53a22}')

I53a33 = Parameter(name = 'I53a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                   texname = '\\text{I53a33}')

I53a36 = Parameter(name = 'I53a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                   texname = '\\text{I53a36}')

I53a63 = Parameter(name = 'I53a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                   texname = '\\text{I53a63}')

I53a66 = Parameter(name = 'I53a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                   texname = '\\text{I53a66}')

I54a33 = Parameter(name = 'I54a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                   texname = '\\text{I54a33}')

I54a36 = Parameter(name = 'I54a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                   texname = '\\text{I54a36}')

I54a63 = Parameter(name = 'I54a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                   texname = '\\text{I54a63}')

I54a66 = Parameter(name = 'I54a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                   texname = '\\text{I54a66}')

I55a33 = Parameter(name = 'I55a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x6)*complexconjugate(tu3x3)',
                   texname = '\\text{I55a33}')

I55a36 = Parameter(name = 'I55a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x6)*complexconjugate(tu3x3)',
                   texname = '\\text{I55a36}')

I55a63 = Parameter(name = 'I55a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x6)*complexconjugate(tu3x3)',
                   texname = '\\text{I55a63}')

I55a66 = Parameter(name = 'I55a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x6)*complexconjugate(tu3x3)',
                   texname = '\\text{I55a66}')

I56a33 = Parameter(name = 'I56a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*td3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                   texname = '\\text{I56a33}')

I56a36 = Parameter(name = 'I56a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*td3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                   texname = '\\text{I56a36}')

I56a63 = Parameter(name = 'I56a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*td3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                   texname = '\\text{I56a63}')

I56a66 = Parameter(name = 'I56a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*td3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                   texname = '\\text{I56a66}')

I57a33 = Parameter(name = 'I57a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)*complexconjugate(yd3x3)',
                   texname = '\\text{I57a33}')

I57a36 = Parameter(name = 'I57a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)*complexconjugate(yd3x3)',
                   texname = '\\text{I57a36}')

I57a63 = Parameter(name = 'I57a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)*complexconjugate(yd3x3)',
                   texname = '\\text{I57a63}')

I57a66 = Parameter(name = 'I57a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)*complexconjugate(yd3x3)',
                   texname = '\\text{I57a66}')

I58a33 = Parameter(name = 'I58a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                   texname = '\\text{I58a33}')

I58a36 = Parameter(name = 'I58a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                   texname = '\\text{I58a36}')

I58a63 = Parameter(name = 'I58a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                   texname = '\\text{I58a63}')

I58a66 = Parameter(name = 'I58a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                   texname = '\\text{I58a66}')

I59a33 = Parameter(name = 'I59a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                   texname = '\\text{I59a33}')

I59a36 = Parameter(name = 'I59a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                   texname = '\\text{I59a36}')

I59a63 = Parameter(name = 'I59a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                   texname = '\\text{I59a63}')

I59a66 = Parameter(name = 'I59a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x6*yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                   texname = '\\text{I59a66}')

I6a33 = Parameter(name = 'I6a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                  texname = '\\text{I6a33}')

I6a36 = Parameter(name = 'I6a36',
                  nature = 'internal',
                  type = 'complex',
                  value = 'complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                  texname = '\\text{I6a36}')

I60a33 = Parameter(name = 'I60a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*yu3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)*complexconjugate(yu3x3)',
                   texname = '\\text{I60a33}')

I60a36 = Parameter(name = 'I60a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*yu3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)*complexconjugate(yu3x3)',
                   texname = '\\text{I60a36}')

I60a63 = Parameter(name = 'I60a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*yu3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)*complexconjugate(yu3x3)',
                   texname = '\\text{I60a63}')

I60a66 = Parameter(name = 'I60a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*yu3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)*complexconjugate(yu3x3)',
                   texname = '\\text{I60a66}')

I61a33 = Parameter(name = 'I61a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x3*complexconjugate(yu3x3)',
                   texname = '\\text{I61a33}')

I61a36 = Parameter(name = 'I61a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x3*complexconjugate(yu3x3)',
                   texname = '\\text{I61a36}')

I62a33 = Parameter(name = 'I62a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x6*yu3x3',
                   texname = '\\text{I62a33}')

I62a36 = Parameter(name = 'I62a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x6*yu3x3',
                   texname = '\\text{I62a36}')

I63a11 = Parameter(name = 'I63a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM1x1*Ru1x1',
                   texname = '\\text{I63a11}')

I63a22 = Parameter(name = 'I63a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM2x2*Ru2x2',
                   texname = '\\text{I63a22}')

I63a33 = Parameter(name = 'I63a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru3x3',
                   texname = '\\text{I63a33}')

I63a36 = Parameter(name = 'I63a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru6x3',
                   texname = '\\text{I63a36}')

I64a33 = Parameter(name = 'I64a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru3x3*complexconjugate(yd3x3)',
                   texname = '\\text{I64a33}')

I64a36 = Parameter(name = 'I64a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru6x3*complexconjugate(yd3x3)',
                   texname = '\\text{I64a36}')

I65a33 = Parameter(name = 'I65a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru3x6*yu3x3',
                   texname = '\\text{I65a33}')

I65a36 = Parameter(name = 'I65a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru6x6*yu3x3',
                   texname = '\\text{I65a36}')

I66a11 = Parameter(name = 'I66a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM1x1*Ru1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I66a11}')

I66a22 = Parameter(name = 'I66a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM2x2*Ru2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I66a22}')

I66a33 = Parameter(name = 'I66a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I66a33}')

I66a36 = Parameter(name = 'I66a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru6x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I66a36}')

I66a63 = Parameter(name = 'I66a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I66a63}')

I66a66 = Parameter(name = 'I66a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru6x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I66a66}')

I67a33 = Parameter(name = 'I67a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru3x3*complexconjugate(Rd3x6)*complexconjugate(td3x3)',
                   texname = '\\text{I67a33}')

I67a36 = Parameter(name = 'I67a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru6x3*complexconjugate(Rd3x6)*complexconjugate(td3x3)',
                   texname = '\\text{I67a36}')

I67a63 = Parameter(name = 'I67a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru3x3*complexconjugate(Rd6x6)*complexconjugate(td3x3)',
                   texname = '\\text{I67a63}')

I67a66 = Parameter(name = 'I67a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru6x3*complexconjugate(Rd6x6)*complexconjugate(td3x3)',
                   texname = '\\text{I67a66}')

I68a33 = Parameter(name = 'I68a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I68a33}')

I68a36 = Parameter(name = 'I68a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru6x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I68a36}')

I68a63 = Parameter(name = 'I68a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I68a63}')

I68a66 = Parameter(name = 'I68a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru6x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I68a66}')

I69a33 = Parameter(name = 'I69a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru3x6*tu3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I69a33}')

I69a36 = Parameter(name = 'I69a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru6x6*tu3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I69a36}')

I69a63 = Parameter(name = 'I69a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru3x6*tu3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I69a63}')

I69a66 = Parameter(name = 'I69a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru6x6*tu3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I69a66}')

I7a33 = Parameter(name = 'I7a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yd3x3*complexconjugate(Rd3x3)',
                  texname = '\\text{I7a33}')

I7a36 = Parameter(name = 'I7a36',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yd3x3*complexconjugate(Rd6x3)',
                  texname = '\\text{I7a36}')

I70a33 = Parameter(name = 'I70a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru3x3*yd3x3*complexconjugate(Rd3x3)*complexconjugate(yd3x3)',
                   texname = '\\text{I70a33}')

I70a36 = Parameter(name = 'I70a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru6x3*yd3x3*complexconjugate(Rd3x3)*complexconjugate(yd3x3)',
                   texname = '\\text{I70a36}')

I70a63 = Parameter(name = 'I70a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru3x3*yd3x3*complexconjugate(Rd6x3)*complexconjugate(yd3x3)',
                   texname = '\\text{I70a63}')

I70a66 = Parameter(name = 'I70a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru6x3*yd3x3*complexconjugate(Rd6x3)*complexconjugate(yd3x3)',
                   texname = '\\text{I70a66}')

I71a33 = Parameter(name = 'I71a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru3x6*yu3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I71a33}')

I71a36 = Parameter(name = 'I71a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru6x6*yu3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I71a36}')

I71a63 = Parameter(name = 'I71a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru3x6*yu3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I71a63}')

I71a66 = Parameter(name = 'I71a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru6x6*yu3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I71a66}')

I72a33 = Parameter(name = 'I72a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru3x3*yu3x3*complexconjugate(Rd3x3)*complexconjugate(yu3x3)',
                   texname = '\\text{I72a33}')

I72a36 = Parameter(name = 'I72a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru6x3*yu3x3*complexconjugate(Rd3x3)*complexconjugate(yu3x3)',
                   texname = '\\text{I72a36}')

I72a63 = Parameter(name = 'I72a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru3x3*yu3x3*complexconjugate(Rd6x3)*complexconjugate(yu3x3)',
                   texname = '\\text{I72a63}')

I72a66 = Parameter(name = 'I72a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru6x3*yu3x3*complexconjugate(Rd6x3)*complexconjugate(yu3x3)',
                   texname = '\\text{I72a66}')

I73a33 = Parameter(name = 'I73a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru3x6*yu3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I73a33}')

I73a36 = Parameter(name = 'I73a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru6x6*yu3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I73a36}')

I73a63 = Parameter(name = 'I73a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru3x6*yu3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I73a63}')

I73a66 = Parameter(name = 'I73a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru6x6*yu3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I73a66}')

I74a11 = Parameter(name = 'I74a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru1x1*complexconjugate(Ru1x1)',
                   texname = '\\text{I74a11}')

I74a22 = Parameter(name = 'I74a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru2x2*complexconjugate(Ru2x2)',
                   texname = '\\text{I74a22}')

I74a33 = Parameter(name = 'I74a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x3*complexconjugate(Ru3x3)',
                   texname = '\\text{I74a33}')

I74a36 = Parameter(name = 'I74a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x3*complexconjugate(Ru3x3)',
                   texname = '\\text{I74a36}')

I74a63 = Parameter(name = 'I74a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x3*complexconjugate(Ru6x3)',
                   texname = '\\text{I74a63}')

I74a66 = Parameter(name = 'I74a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x3*complexconjugate(Ru6x3)',
                   texname = '\\text{I74a66}')

I75a33 = Parameter(name = 'I75a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x6*complexconjugate(Ru3x6)',
                   texname = '\\text{I75a33}')

I75a36 = Parameter(name = 'I75a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x6*complexconjugate(Ru3x6)',
                   texname = '\\text{I75a36}')

I75a44 = Parameter(name = 'I75a44',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru4x4*complexconjugate(Ru4x4)',
                   texname = '\\text{I75a44}')

I75a55 = Parameter(name = 'I75a55',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru5x5*complexconjugate(Ru5x5)',
                   texname = '\\text{I75a55}')

I75a63 = Parameter(name = 'I75a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x6*complexconjugate(Ru6x6)',
                   texname = '\\text{I75a63}')

I75a66 = Parameter(name = 'I75a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x6*complexconjugate(Ru6x6)',
                   texname = '\\text{I75a66}')

I76a33 = Parameter(name = 'I76a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x3*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                   texname = '\\text{I76a33}')

I76a36 = Parameter(name = 'I76a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x3*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                   texname = '\\text{I76a36}')

I76a63 = Parameter(name = 'I76a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x3*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                   texname = '\\text{I76a63}')

I76a66 = Parameter(name = 'I76a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x3*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                   texname = '\\text{I76a66}')

I77a33 = Parameter(name = 'I77a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x3*complexconjugate(Ru3x6)*complexconjugate(tu3x3)',
                   texname = '\\text{I77a33}')

I77a36 = Parameter(name = 'I77a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x3*complexconjugate(Ru3x6)*complexconjugate(tu3x3)',
                   texname = '\\text{I77a36}')

I77a63 = Parameter(name = 'I77a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x3*complexconjugate(Ru6x6)*complexconjugate(tu3x3)',
                   texname = '\\text{I77a63}')

I77a66 = Parameter(name = 'I77a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x3*complexconjugate(Ru6x6)*complexconjugate(tu3x3)',
                   texname = '\\text{I77a66}')

I78a33 = Parameter(name = 'I78a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x6*tu3x3*complexconjugate(Ru3x3)',
                   texname = '\\text{I78a33}')

I78a36 = Parameter(name = 'I78a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x6*tu3x3*complexconjugate(Ru3x3)',
                   texname = '\\text{I78a36}')

I78a63 = Parameter(name = 'I78a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x6*tu3x3*complexconjugate(Ru6x3)',
                   texname = '\\text{I78a63}')

I78a66 = Parameter(name = 'I78a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x6*tu3x3*complexconjugate(Ru6x3)',
                   texname = '\\text{I78a66}')

I79a33 = Parameter(name = 'I79a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x6*yu3x3*complexconjugate(Ru3x3)',
                   texname = '\\text{I79a33}')

I79a36 = Parameter(name = 'I79a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x6*yu3x3*complexconjugate(Ru3x3)',
                   texname = '\\text{I79a36}')

I79a63 = Parameter(name = 'I79a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x6*yu3x3*complexconjugate(Ru6x3)',
                   texname = '\\text{I79a63}')

I79a66 = Parameter(name = 'I79a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x6*yu3x3*complexconjugate(Ru6x3)',
                   texname = '\\text{I79a66}')

I8a11 = Parameter(name = 'I8a11',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd1x1*complexconjugate(Rd1x1)',
                  texname = '\\text{I8a11}')

I8a22 = Parameter(name = 'I8a22',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd2x2*complexconjugate(Rd2x2)',
                  texname = '\\text{I8a22}')

I8a33 = Parameter(name = 'I8a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd3x3*complexconjugate(Rd3x3)',
                  texname = '\\text{I8a33}')

I8a36 = Parameter(name = 'I8a36',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd6x3*complexconjugate(Rd3x3)',
                  texname = '\\text{I8a36}')

I8a63 = Parameter(name = 'I8a63',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd3x3*complexconjugate(Rd6x3)',
                  texname = '\\text{I8a63}')

I8a66 = Parameter(name = 'I8a66',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd6x3*complexconjugate(Rd6x3)',
                  texname = '\\text{I8a66}')

I80a33 = Parameter(name = 'I80a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x3*yu3x3*complexconjugate(Ru3x3)*complexconjugate(yu3x3)',
                   texname = '\\text{I80a33}')

I80a36 = Parameter(name = 'I80a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x3*yu3x3*complexconjugate(Ru3x3)*complexconjugate(yu3x3)',
                   texname = '\\text{I80a36}')

I80a63 = Parameter(name = 'I80a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x3*yu3x3*complexconjugate(Ru6x3)*complexconjugate(yu3x3)',
                   texname = '\\text{I80a63}')

I80a66 = Parameter(name = 'I80a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x3*yu3x3*complexconjugate(Ru6x3)*complexconjugate(yu3x3)',
                   texname = '\\text{I80a66}')

I81a33 = Parameter(name = 'I81a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x6*yu3x3*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                   texname = '\\text{I81a33}')

I81a36 = Parameter(name = 'I81a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x6*yu3x3*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                   texname = '\\text{I81a36}')

I81a63 = Parameter(name = 'I81a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x6*yu3x3*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                   texname = '\\text{I81a63}')

I81a66 = Parameter(name = 'I81a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x6*yu3x3*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                   texname = '\\text{I81a66}')

I82a11 = Parameter(name = 'I82a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I82a11}')

I82a22 = Parameter(name = 'I82a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I82a22}')

I82a33 = Parameter(name = 'I82a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I82a33}')

I82a36 = Parameter(name = 'I82a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I82a36}')

I83a33 = Parameter(name = 'I83a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*complexconjugate(Rd3x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I83a33}')

I83a36 = Parameter(name = 'I83a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*complexconjugate(Rd6x6)*complexconjugate(yd3x3)',
                   texname = '\\text{I83a36}')

I84a33 = Parameter(name = 'I84a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*yu3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I84a33}')

I84a36 = Parameter(name = 'I84a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*yu3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I84a36}')

I85a11 = Parameter(name = 'I85a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(Rl1x1)',
                   texname = '\\text{I85a11}')

I85a22 = Parameter(name = 'I85a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(Rl2x2)',
                   texname = '\\text{I85a22}')

I85a33 = Parameter(name = 'I85a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(Rl3x3)',
                   texname = '\\text{I85a33}')

I85a36 = Parameter(name = 'I85a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(Rl6x3)',
                   texname = '\\text{I85a36}')

I86a33 = Parameter(name = 'I86a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(Rl3x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I86a33}')

I86a36 = Parameter(name = 'I86a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(Rl6x6)*complexconjugate(ye3x3)',
                   texname = '\\text{I86a36}')

I87a11 = Parameter(name = 'I87a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(Rn1x1)',
                   texname = '\\text{I87a11}')

I87a22 = Parameter(name = 'I87a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(Rn2x2)',
                   texname = '\\text{I87a22}')

I87a33 = Parameter(name = 'I87a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(Rn3x3)',
                   texname = '\\text{I87a33}')

I88a33 = Parameter(name = 'I88a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'ye3x3*complexconjugate(Rn3x3)',
                   texname = '\\text{I88a33}')

I89a11 = Parameter(name = 'I89a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                   texname = '\\text{I89a11}')

I89a22 = Parameter(name = 'I89a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                   texname = '\\text{I89a22}')

I89a33 = Parameter(name = 'I89a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                   texname = '\\text{I89a33}')

I89a36 = Parameter(name = 'I89a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                   texname = '\\text{I89a36}')

I9a33 = Parameter(name = 'I9a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd3x6*complexconjugate(Rd3x6)',
                  texname = '\\text{I9a33}')

I9a36 = Parameter(name = 'I9a36',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd6x6*complexconjugate(Rd3x6)',
                  texname = '\\text{I9a36}')

I9a44 = Parameter(name = 'I9a44',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd4x4*complexconjugate(Rd4x4)',
                  texname = '\\text{I9a44}')

I9a55 = Parameter(name = 'I9a55',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd5x5*complexconjugate(Rd5x5)',
                  texname = '\\text{I9a55}')

I9a63 = Parameter(name = 'I9a63',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd3x6*complexconjugate(Rd6x6)',
                  texname = '\\text{I9a63}')

I9a66 = Parameter(name = 'I9a66',
                  nature = 'internal',
                  type = 'complex',
                  value = 'Rd6x6*complexconjugate(Rd6x6)',
                  texname = '\\text{I9a66}')

I90a33 = Parameter(name = 'I90a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(CKM3x3)*complexconjugate(Ru3x6)*complexconjugate(yu3x3)',
                   texname = '\\text{I90a33}')

I90a36 = Parameter(name = 'I90a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'complexconjugate(CKM3x3)*complexconjugate(Ru6x6)*complexconjugate(yu3x3)',
                   texname = '\\text{I90a36}')

I91a33 = Parameter(name = 'I91a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                   texname = '\\text{I91a33}')

I91a36 = Parameter(name = 'I91a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'yd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                   texname = '\\text{I91a36}')

I92a11 = Parameter(name = 'I92a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM1x1*Ru1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I92a11}')

I92a22 = Parameter(name = 'I92a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM2x2*Ru2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I92a22}')

I92a33 = Parameter(name = 'I92a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I92a33}')

I92a36 = Parameter(name = 'I92a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru6x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I92a36}')

I92a63 = Parameter(name = 'I92a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I92a63}')

I92a66 = Parameter(name = 'I92a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'CKM3x3*Ru6x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I92a66}')

I93a11 = Parameter(name = 'I93a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rn1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I93a11}')

I93a22 = Parameter(name = 'I93a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rn2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I93a22}')

I93a33 = Parameter(name = 'I93a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rn3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I93a33}')

I93a36 = Parameter(name = 'I93a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rn3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I93a36}')

I94a11 = Parameter(name = 'I94a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(CKM1x1)*complexconjugate(Ru1x1)',
                   texname = '\\text{I94a11}')

I94a22 = Parameter(name = 'I94a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(CKM2x2)*complexconjugate(Ru2x2)',
                   texname = '\\text{I94a22}')

I94a33 = Parameter(name = 'I94a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                   texname = '\\text{I94a33}')

I94a36 = Parameter(name = 'I94a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                   texname = '\\text{I94a36}')

I94a63 = Parameter(name = 'I94a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru3x3)',
                   texname = '\\text{I94a63}')

I94a66 = Parameter(name = 'I94a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(CKM3x3)*complexconjugate(Ru6x3)',
                   texname = '\\text{I94a66}')

I95a11 = Parameter(name = 'I95a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rn1x1)',
                   texname = '\\text{I95a11}')

I95a22 = Parameter(name = 'I95a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rn2x2)',
                   texname = '\\text{I95a22}')

I95a33 = Parameter(name = 'I95a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rn3x3)',
                   texname = '\\text{I95a33}')

I95a36 = Parameter(name = 'I95a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rn3x3)',
                   texname = '\\text{I95a36}')

I96a11 = Parameter(name = 'I96a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd1x1*complexconjugate(Rd1x1)',
                   texname = '\\text{I96a11}')

I96a22 = Parameter(name = 'I96a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd2x2*complexconjugate(Rd2x2)',
                   texname = '\\text{I96a22}')

I96a33 = Parameter(name = 'I96a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I96a33}')

I96a36 = Parameter(name = 'I96a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd3x3)',
                   texname = '\\text{I96a36}')

I96a63 = Parameter(name = 'I96a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd3x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I96a63}')

I96a66 = Parameter(name = 'I96a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rd6x3*complexconjugate(Rd6x3)',
                   texname = '\\text{I96a66}')

I97a11 = Parameter(name = 'I97a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl1x1*complexconjugate(Rl1x1)',
                   texname = '\\text{I97a11}')

I97a22 = Parameter(name = 'I97a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl2x2*complexconjugate(Rl2x2)',
                   texname = '\\text{I97a22}')

I97a33 = Parameter(name = 'I97a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I97a33}')

I97a36 = Parameter(name = 'I97a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl3x3)',
                   texname = '\\text{I97a36}')

I97a63 = Parameter(name = 'I97a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl3x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I97a63}')

I97a66 = Parameter(name = 'I97a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Rl6x3*complexconjugate(Rl6x3)',
                   texname = '\\text{I97a66}')

I98a11 = Parameter(name = 'I98a11',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru1x1*complexconjugate(Ru1x1)',
                   texname = '\\text{I98a11}')

I98a22 = Parameter(name = 'I98a22',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru2x2*complexconjugate(Ru2x2)',
                   texname = '\\text{I98a22}')

I98a33 = Parameter(name = 'I98a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x3*complexconjugate(Ru3x3)',
                   texname = '\\text{I98a33}')

I98a36 = Parameter(name = 'I98a36',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x3*complexconjugate(Ru3x3)',
                   texname = '\\text{I98a36}')

I98a63 = Parameter(name = 'I98a63',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru3x3*complexconjugate(Ru6x3)',
                   texname = '\\text{I98a63}')

I98a66 = Parameter(name = 'I98a66',
                   nature = 'internal',
                   type = 'complex',
                   value = 'Ru6x3*complexconjugate(Ru6x3)',
                   texname = '\\text{I98a66}')

I99a33 = Parameter(name = 'I99a33',
                   nature = 'internal',
                   type = 'complex',
                   value = 'ye3x3',
                   texname = '\\text{I99a33}')

