#! /usr/bin/env python
# encoding: utf-8
# Copyright (c) 2011-2013, 2015-2022 Harald Klimach <harald@klimachs.de>
# Copyright (c) 2011 Konstantin Kleinheinz <k.kleinheinz@grs-sim.de>
# Copyright (c) 2011 Manuel Hasert <m.hasert@grs-sim.de>
# Copyright (c) 2011 Jens Zudrop <j.zudrop@grs-sim.de>
# Copyright (c) 2012 Kartik Jain <k.jain@grs-sim.de>
# Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
#
# Parts of this file were written by Harald Klimach, Konstantin Kleinheinz,
# Manuel Hasert, Jens Zudrop and Kartik Jain for German Research School of
# Simulation Sciences.
# Parts of this file were written by Harald Klimach and Kannan Masilamani
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

import os

APPNAME = 'aotus'
VERSION = '1'

top = '.'
out = 'build'

def append_aotmodpaths(ctx):
    ''' Add directories with Python modules for waf to sys.path. '''
    import sys
    myabspath = ctx.path.abspath()
    if not myabspath in sys.path:
        sys.path.append(myabspath)

def options(opt):
    append_aotmodpaths(opt)
    opt.load('fortran_compiler')
    opt.load('fortran_language')
    opt.load('compiler_c')
    opt.load('clang_compilation_database')
    opt.load('waf_unit_test')
    opt.load('utest_results')
    opt.load('make_fordoc')
    opt.add_option('--command_sequence', action='store_true', default=False,
                   help='Collect all executed commands into a single file.',
                   dest='cmdsequence')

def configure(conf):

    conf.load('make_fordoc')

    # Use a function for the first part to make it callable
    # from parent projects without setting the flags.
    subconf(conf)

    from fortran_compiler import set_fc_flags
    osfcflags = conf.env.FCFLAGS

    # Flags for the default (production) variant
    set_fc_flags(conf, ['optimize', 'warn'], osfcflags)

    # Set flags for the debug variant
    conf.setenv('debug',conf.env)
    set_fc_flags(conf, ['standard', 'warn', 'w2e', 'debug'],
                 osfcflags)

MKSTEMP_FRAG = '''
#include <stdlib.h>
int main(int argc, char **argv) {
  char fname[] = "aotempXXXXXX";
  return mkstemp(fname);
}
'''

POPEN_FRAG = '''
#include <stdio.h>
int main(int argc, char **argv) {
  FILE *channel;
  const char fname[] = "aotchannel";
  const char fmode[] = "r";
  channel = popen(fname, fmode);
}
'''

SRANDOM_FRAG = '''
#include <stdlib.h>
int main(int argc, char **argv) {
  int seed = 42;
  srandom(seed);
}
'''

def subconf(conf):
    """
    Configure parts, which are relevant, even when called
    from parent wscripts.
    Useful to restrict parent recursions to just this part
    of the configuration.
    """

    append_aotmodpaths(conf)
    conf.load('waf_unit_test')

    conf.env.fordurl_aotus = 'https://geb.inf.tu-dresden.de/doxy/aotus/'

    # Load the C compiler information
    conf.setenv('cenv',conf.env)
    if not conf.options.bench:
        conf.load('compiler_c')
    else:
        print('loading bench c compiler')
        from bench_compiler import set_bench_c_compiler
        set_bench_c_compiler(conf)
        conf.load('c')

    # Load the Fortran compiler information
    conf.setenv('')
    conf.env.DEST_OS = conf.all_envs['cenv'].DEST_OS
    conf.env.CC_NAME = conf.all_envs['cenv'].CC_NAME
    conf.load('fortran_compiler')
    conf.check_fortran()

    conf.setenv('cenv')

    conf.env.append_unique('DEFINES', ['LUA_ANSI'])

    # Do not change the DEFINES themselves, but use a temporary copy.
    tmpenv = conf.env.derive()
    tmpenv.detach()

    conf.start_msg('Can use POSIX features in Lua')
    conf.check_cc(header_name=['stdlib.h', 'unistd.h', 'stdio.h', 'math.h'],
                  define_name='POSIXH',
                  mandatory=False)
    if conf.is_defined('POSIXH'):
        conf.check_cc(fragment=MKSTEMP_FRAG,
                      define_name='MKSTEMP',
                      mandatory=False)
        conf.check_cc(fragment=POPEN_FRAG,
                      define_name='POPEN',
                      mandatory=False)
        conf.check_cc(fragment=SRANDOM_FRAG,
                      define_name='SRANDOM',
                      mandatory=False)

    if ( conf.is_defined('POPEN') and
         conf.is_defined('MKSTEMP') and
         conf.is_defined('SRANDOM') ):

        conf.env = tmpenv
        conf.all_envs['cenv'].DEFINES_LUA_POSIX = ['LUA_USE_POSIX']
        conf.end_msg('yes')

    else:
        conf.env = tmpenv
        conf.end_msg('NO')

    # Only required to build the Lua interpreter
    conf.check_cc(lib='m', uselib_store='MATH', mandatory=False)

    conf.setenv('cenv_debug', conf.env)
    conf.setenv('cenv_profile', conf.env)
    conf.setenv('')

    import fortran_language

    # Check for ISO_C_Binding support
    fortran_language.supports_iso_c(conf)

    # Check for higher precision real kinds:
    fortran_language.supports_quad_kind(conf = conf, mandatory = False)
    fortran_language.supports_xdble_kind(conf = conf, mandatory = False)


def build(bld):
    from waflib import Runner
    orig_refill = Runner.Parallel.refill_task_list
    def catch_circular(self):
        from waflib import Errors
        import fortran_circular
        try:
            return(orig_refill(self))
        except Errors.WafError as e:
            if not fortran_circular.find_circular_dependency(str(e)):
                raise(e)
    Runner.Parallel.refill_task_list = catch_circular

    append_aotmodpaths(bld)
    if bld.options.cmdsequence:
        import waflib.extras.command_sequence

    core_sources = ['external/lua-5.4.4/src/lapi.c',
                    'external/lua-5.4.4/src/lcode.c',
                    'external/lua-5.4.4/src/lctype.c',
                    'external/lua-5.4.4/src/ldebug.c',
                    'external/lua-5.4.4/src/ldo.c',
                    'external/lua-5.4.4/src/ldump.c',
                    'external/lua-5.4.4/src/lfunc.c',
                    'external/lua-5.4.4/src/lgc.c',
                    'external/lua-5.4.4/src/llex.c',
                    'external/lua-5.4.4/src/lmem.c',
                    'external/lua-5.4.4/src/lobject.c',
                    'external/lua-5.4.4/src/lopcodes.c',
                    'external/lua-5.4.4/src/lparser.c',
                    'external/lua-5.4.4/src/lstate.c',
                    'external/lua-5.4.4/src/lstring.c',
                    'external/lua-5.4.4/src/ltable.c',
                    'external/lua-5.4.4/src/ltm.c',
                    'external/lua-5.4.4/src/lundump.c',
                    'external/lua-5.4.4/src/lvm.c',
                    'external/lua-5.4.4/src/lzio.c']
    lib_sources = ['external/lua-5.4.4/src/lauxlib.c',
                   'external/lua-5.4.4/src/lbaselib.c',
                   'external/lua-5.4.4/src/lcorolib.c',
                   'external/lua-5.4.4/src/ldblib.c',
                   'external/lua-5.4.4/src/liolib.c',
                   'external/lua-5.4.4/src/lmathlib.c',
                   'external/lua-5.4.4/src/loslib.c',
                   'external/lua-5.4.4/src/ltablib.c',
                   'external/lua-5.4.4/src/lstrlib.c',
                   'external/lua-5.4.4/src/lutf8lib.c',
                   'external/lua-5.4.4/src/loadlib.c',
                   'external/lua-5.4.4/src/linit.c']
    lua_sources = ['external/lua-5.4.4/src/lua.c']
    luac_sources = ['external/lua-5.4.4/src/luac.c']

    wrap_sources = ['LuaFortran/wrap_lua_dump.c']

    flu_sources = ['LuaFortran/lua_fif.f90',
                   'LuaFortran/dump_lua_fif_module.f90',
                   'LuaFortran/lua_parameters.f90',
                   'LuaFortran/flu_kinds_module.f90',
                   'LuaFortran/flu_binding.f90']

    aotus_sources = ['source/aotus_module.f90',
                     'source/aot_err_module.f90',
                     'source/aot_fun_module.f90',
                     'source/aot_fun_declaration_module.f90',
                     'source/aot_table_module.f90',
                     'source/aot_table_ops_module.f90',
                     'source/aot_top_module.f90',
                     'source/aot_out_module.f90',
                     'source/aot_out_general_module.f90',
                     'source/aot_path_module.f90',
                     'source/aot_references_module.f90',
                     'source/aot_vector_module.f90']

    if bld.env['fortsupp_quad_kind'] > 0:
        aotus_sources += ['source/quadruple/aot_quadruple_fun_module.f90']
        aotus_sources += ['source/quadruple/aot_quadruple_table_module.f90']
        aotus_sources += ['source/quadruple/aot_quadruple_top_module.f90']
        aotus_sources += ['source/quadruple/aot_quadruple_out_module.f90']
        aotus_sources += ['source/quadruple/aot_quadruple_vector_module.f90']
    else:
        aotus_sources += ['source/quadruple/dummy_quadruple_fun_module.f90']
        aotus_sources += ['source/quadruple/dummy_quadruple_table_module.f90']
        aotus_sources += ['source/quadruple/dummy_quadruple_top_module.f90']
        aotus_sources += ['source/quadruple/dummy_quadruple_out_module.f90']
        aotus_sources += ['source/quadruple/dummy_quadruple_vector_module.f90']

    if bld.env['fortsupp_xdble_kind'] > 0:
        aotus_sources += ['source/extdouble/aot_extdouble_fun_module.f90']
        aotus_sources += ['source/extdouble/aot_extdouble_table_module.f90']
        aotus_sources += ['source/extdouble/aot_extdouble_top_module.f90']
        aotus_sources += ['source/extdouble/aot_extdouble_out_module.f90']
        aotus_sources += ['source/extdouble/aot_extdouble_vector_module.f90']
    else:
        aotus_sources += ['source/extdouble/dummy_extdouble_fun_module.f90']
        aotus_sources += ['source/extdouble/dummy_extdouble_table_module.f90']
        aotus_sources += ['source/extdouble/dummy_extdouble_top_module.f90']
        aotus_sources += ['source/extdouble/dummy_extdouble_out_module.f90']
        aotus_sources += ['source/extdouble/dummy_extdouble_vector_module.f90']

    if bld.cmd == 'docu':
        from waflib.extras.make_fordoc import gendoc
        import os
        source_nodes = []
        for source in aotus_sources+flu_sources:
            resource = bld.path.find_resource(source)
            if resource:
                source_nodes.append(resource)
        tgt = bld.path.get_bld().make_node('docu/modules.json')
        bld( rule = gendoc,
             source = source_nodes,
             target = tgt,
             extern = [],
             mainpage = os.path.join(bld.top_dir, 'aotus', 'aot_mainpage.md') )
        bld.env.fordext_aotus = tgt
        return None

    # C parts
    bld(
        features = 'c',
        source = core_sources + lib_sources,
        use = ['LUA_POSIX'],
        target = 'luaobjs')

    bld(
        features = 'c cstlib',
        use = 'luaobjs',
        target = 'lualib')

    bld(
        features = 'c',
        source = wrap_sources,
        use = 'luaobjs',
        includes = 'external/lua-5.4.4/src',
        target = 'wrapobjs')

    ## Building the lua interpreter (usually not needed).
    ## Only built if libm available.
    if 'LIB_MATH' in bld.all_envs['cenv']:
      bld(
          features = 'c cprogram',
          use = ['lualib', 'MATH'],
          source = lua_sources,
          target = 'lua')

    # Fortran parts
    bld(
        features = 'fc',
        source = flu_sources,
        target = 'fluobjs')

    bld(
        features = 'fc fcstlib',
        use = ['luaobjs', 'fluobjs', 'wrapobjs'],
        target = 'flu')

    bld(
        features = 'fc fcstlib',
        source = aotus_sources,
        use = ['luaobjs', 'fluobjs', 'wrapobjs'],
        target = 'aotus')

    bld(
        features = 'fc fcprogram',
        source = ['sample/aotus_sample.f90'],
        use = 'aotus',
        target = 'aotus_sample')

    bld(
        features = 'fc fcprogram',
        source = ['LuaFortran/examples/test.f90'],
        use = 'flu',
        target = 'flu_sample')

    from waflib.extras import utest_results
    utest_results.utests(bld, 'aotus')
    if bld.env['fortsupp_quad_kind'] > 0:
        utest_results.utests(bld, use = 'aotus', path = 'utests/quadruple')
    bld.add_post_fun(utest_results.summary)

    # install_files actually only done, if in install mode.
    # However, the if here avoids the ant_glob in the build directory
    # to be run if not in the install phase...
    if bld.cmd == 'install':
        bld.install_files('${PREFIX}/include',
                          bld.path.get_bld().ant_glob('*.mod'))
        bld.install_files('${PREFIX}/lib', 'libaotus.a')


from waflib.Build import BuildContext
from waflib import TaskGen

# Modify Fortran tasks to not contain SHLIB and STLIB markers if not
# explicitly requested.
@TaskGen.feature('fcprogram', 'fcshlib', 'fcstlib')
@TaskGen.before_method('process_use')
@TaskGen.after_method('apply_link')
def kill_marker_flags(self):
    if not self.env.LIB and not self.env.LIBPATH:
        self.env.FCSHLIB_MARKER = []
    if not self.env.STLIB and not self.env.STLIBPATH:
        self.env.FCSTLIB_MARKER = []

# Modify C tasks to use a dedicated C environment.
@TaskGen.feature('c', 'cstlib', 'cprogram', 'cxx')
@TaskGen.before('process_rule')
def enter_cenv(self):
    try:
        self.env = self.bld.all_envs['cenv_'+self.bld.variant].derive()
    except KeyError:
        self.env = self.bld.all_envs['cenv'].derive()


# A class to describe the debug variant
class debug(BuildContext):
    "Build a debug executable"
    cmd = 'debug'
    variant = 'debug'

from waflib import Logs
class Dumper(BuildContext):
    """Create a Makefile from the configured project."""
    fun = 'dump'
    cmd = 'dump'

def dump(bld):
    import os
    # Do not perform tests when dumping a Makefile.
    bld.options.no_tests = True
    # call the build function as if a real build were performed
    build(bld)

    from waflib import Task
    bld.commands = []
    bld.targets = []
    bld.make_ins = []

    # store the command executed
    def exec_command(self, *k, **kw):
        self.command_executed = k[0]
        self.path = kw['cwd'] or self.generator.bld.cwd
        for x in self.outputs:
            x.write('dummy', 'wb')
        return False
    Task.TaskBase.exec_command = exec_command

    # perform a fake build, and accumulate the makefile bits
    old_process = Task.TaskBase.process
    def process(self):
        old_process(self)

        lst = []
        for x in self.outputs:
            lst.append(x.path_from(self.generator.bld.bldnode))
        bld.targets.extend(lst)
        lst.append(':')
        for x in self.inputs + self.dep_nodes + self.generator.bld.node_deps.get(self.uid(), []):
            lst.append(x.path_from(self.generator.bld.bldnode))
            bld.make_ins.append(x.path_from(self.generator.bld.bldnode))
        try:
            if isinstance(self.command_executed, list):
                self.command_executed = ' '.join(self.command_executed)
            # Replace absolute paths.
            self.command_executed = self.command_executed.replace(bld.out_dir, '.')
            self.command_executed = self.command_executed.replace(bld.top_dir, os.path.relpath(bld.top_dir, bld.out_dir))
        except Exception as e:
            print(e)
        else:
            bld.commands.append(' '.join(lst))
            bld.commands.append('\t{0}\n'.format(self.command_executed).replace(bld.out_dir, '.'))

    Task.TaskBase.process = process

    # write the makefile after the build is complete
    def output_makefile(self):
        self.commands.insert(0, "all: %s\n" % " ".join([x for x in self.targets if x not in self.make_ins]))
        self.commands.insert(0, ".SUFFIXES:\n")
        self.commands.insert(0, ".POSIX:")
        self.commands.append('clean:')
        targetlist = ' '.join(self.targets)
        self.commands.append('\trm -f {0}\n'.format(targetlist))
        node = self.bldnode.make_node('Makefile')
        node.write("\n".join(self.commands))
        Logs.warn('Wrote %s' % node.abspath())
        for dummy in self.targets:
            os.remove(os.path.join(self.out_dir,dummy))

    bld.add_post_fun(output_makefile)
