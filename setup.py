#!/usr/bin/env python
# -*- coding: latin-1 -*-

# Hedge - the Hybrid'n'Easy DG Environment
# Copyright (C) 2007 Andreas Kloeckner
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import glob
import os
import os.path
import sys

def main():
    try:
        conf = {}
        execfile("siteconf.py", conf)
    except IOError:
        print "*** Please run configure first."
        sys.exit(1)

    from distutils.core import setup,Extension

    def non_matching_config():
        print "*** The version of your configuration template does not match"
        print "*** the version of the setup script. Please re-run configure."
        sys.exit(1)

    if "PYRTICLE_CONF_TEMPLATE_VERSION" not in conf:
        non_matching_config()

    if conf["PYRTICLE_CONF_TEMPLATE_VERSION"] != 1:
        non_matching_config()

    INCLUDE_DIRS = ["src/cpp"] \
            + conf["BOOST_MATH_TOOLKIT_INCLUDE_DIRS"] \
            + conf["BOOST_BINDINGS_INCLUDE_DIRS"] \
            + conf["BOOST_INCLUDE_DIRS"]

    LIBRARY_DIRS = conf["BOOST_LIBRARY_DIRS"]
    LIBRARIES = conf["BPL_LIBRARIES"]

    EXTRA_DEFINES = {}
    EXTRA_INCLUDE_DIRS = []
    EXTRA_LIBRARY_DIRS = []
    EXTRA_LIBRARIES = []

    if conf["HAVE_MPI"]:
        EXTRA_DEFINES["USE_MPI"] = 1
        EXTRA_DEFINES["OMPI_SKIP_MPICXX"] = 1
        LIBRARIES.extend(conf["BOOST_MPI_LIBRARIES"])

        cvars = sysconfig.get_config_vars()
        cvars["CC"] = conf["MPICC"]
        cvars["CXX"] = conf["MPICXX"]

    def handle_component(comp):
        if conf["USE_"+comp]:
            EXTRA_DEFINES["USE_"+comp] = 1
            EXTRA_INCLUDE_DIRS.extend(conf[comp+"_INCLUDE_DIRS"])
            EXTRA_LIBRARY_DIRS.extend(conf[comp+"_LIBRARY_DIRS"])
            EXTRA_LIBRARIES.extend(conf[comp+"_LIBRARIES"])

    setup(name="pyrticle",
          version="0.90",
          description="A high-order PIC code using Hedge",
          author=u"Andreas Kloeckner",
          author_email="inform@tiker.net",
          license = "Proprietary",
          #url="http://news.tiker.net/software/hedge",
          packages=["pyrticle"],
          package_dir={"pyrticle": "src/python"},
          ext_package="pyrticle",
          ext_modules=[
            Extension("_internal", 
                [
                    "src/cpp/meshdata.cpp",
                    "src/cpp/rec_shape.cpp",
                    "src/wrapper/wrap_tools.cpp",
                    "src/wrapper/wrap_meshdata.cpp",
                    "src/wrapper/wrap_reconstructor.cpp",
                    "src/wrapper/wrap_pusher.cpp",
                    "src/wrapper/wrap_pic.cpp", 
                    "src/wrapper/wrap_main.cpp", 
                ],
                include_dirs=INCLUDE_DIRS + EXTRA_INCLUDE_DIRS,
                library_dirs=LIBRARY_DIRS + EXTRA_LIBRARY_DIRS,
                libraries=LIBRARIES + EXTRA_LIBRARIES,
                extra_compile_args=conf["EXTRA_COMPILE_ARGS"],
                define_macros=list(EXTRA_DEFINES.iteritems()),
                )]
         )



if __name__ == '__main__':
    # hack distutils.sysconfig to eliminate debug flags
    # stolen from mpi4py
    import sys
    if not sys.platform.lower().startswith("win"):
        from distutils import sysconfig

        cvars = sysconfig.get_config_vars()
        cflags = cvars.get('OPT')
        if cflags:
            cflags = cflags.split()
            for bad_prefix in ('-g', '-O', '-Wstrict-prototypes'):
                for i, flag in enumerate(cflags):
                    if flag.startswith(bad_prefix):
                        cflags.pop(i)
                        break
                if flag in cflags:
                    cflags.remove(flag)
            cflags.append("-O3")
            cvars['OPT'] = str.join(' ', cflags)
            cvars["CFLAGS"] = cvars["BASECFLAGS"] + " " + cvars["OPT"]
    # and now call main
    main()
