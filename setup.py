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

try:
    execfile("siteconf.py")
except IOError:
    print "*** Please run configure first."
    sys.exit(1)

from distutils.core import setup,Extension

def non_matching_config():
    print "*** The version of your configuration template does not match"
    print "*** the version of the setup script. Please re-run configure."
    sys.exit(1)

try:
    PYRTICLE_CONF_TEMPLATE_VERSION
except NameError:
    non_matching_config()

if PYRTICLE_CONF_TEMPLATE_VERSION != 1:
    non_matching_config()

INCLUDE_DIRS = ["src/cpp"] \
        + BOOST_INCLUDE_DIRS \
        + BOOST_BINDINGS_INCLUDE_DIRS

LIBRARY_DIRS = BOOST_LIBRARY_DIRS
LIBRARIES = BPL_LIBRARIES

EXTRA_DEFINES = {}
EXTRA_INCLUDE_DIRS = []
EXTRA_LIBRARY_DIRS = []
EXTRA_LIBRARIES = []

def handle_component(comp):
    if globals()["USE_"+comp]:
        globals()["EXTRA_DEFINES"]["USE_"+comp] = 1
        globals()["EXTRA_INCLUDE_DIRS"] += globals()[comp+"_INCLUDE_DIRS"]
        globals()["EXTRA_LIBRARY_DIRS"] += globals()[comp+"_LIBRARY_DIRS"]
        globals()["EXTRA_LIBRARIES"] += globals()[comp+"_LIBRARIES"]

ext_modules=[
        Extension("_internal", 
            ["src/cpp/main.cpp", ],
            include_dirs=INCLUDE_DIRS + EXTRA_INCLUDE_DIRS,
            library_dirs=LIBRARY_DIRS + EXTRA_LIBRARY_DIRS,
            libraries=LIBRARIES + EXTRA_LIBRARIES,
            extra_compile_args=EXTRA_COMPILE_ARGS,
            define_macros=list(EXTRA_DEFINES.iteritems()),
            )]




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
      ext_modules=ext_modules
     )
