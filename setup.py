#!/usr/bin/env python
# Copyright (c) 2018.  Mark E. Madsen <mark@madsenlab.org>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.

"""
Description here

"""

from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages, Command
from setuptools.command.develop import develop
from setuptools.command.install import install

import subprocess
import os
import re
import shutil

# create a decorator that wraps the normal develop and
# install commands to first update the version

def versioned(command_subclass):
    """
    A decorator for ensuring that installations, whether through setup.py install
    or setup.py develop, always update the version number with the Git revision and
    most recent tag.

    :param command_subclass:
    :return:
    """
    VERSION_PY = """
# This file is updated from Git information by running 'python setup.py
# version'.
__version__ = '%s'
"""
    orig_callable = command_subclass.run

    def modified_callable(self):
        if not os.path.isdir(".git"):
            print ("This does not appear to be a Git repository.")
            return
        try:
            p = subprocess.Popen(["git", "describe",
                                  "--tags", "--always"],
                                 stdout=subprocess.PIPE)
        except EnvironmentError:
            print ("unable to run git, leaving network-drift/_version.py alone")
            return
        stdout = p.communicate()[0]
        if p.returncode != 0:
            print ("unable to run git, leaving network-drift/_version.py alone")
            return
        # our tags are like:  v2.2
        ver = stdout[len("v"):].strip()
        f = open("seriationct/_version.py", "w")
        f.write(VERSION_PY % ver)
        f.close()
        print("updated _version.py to '%s'" % ver)
        orig_callable(self)




    command_subclass.run = modified_callable
    return command_subclass


@versioned
class CustomDevelopClass(develop):
    pass

@versioned
class CustomInstallClass(install):

    def run(self):
        install.run(self)

def get_version():
    try:
        f = open("seriationct/_version.py")
    except EnvironmentError:
        return None
    for line in f.readlines():
        mo = re.match("__version__ = '([^']+)'", line)
        if mo:
            ver = mo.group(1)
            return ver
    return None


class VersionUpdate(Command):
    """setuptools Command"""
    description = "update version number"
    user_options = []


    def initialize_options(self):
        pass


    def finalize_options(self):
        pass


    def run(self):
        VERSION = """
# This file is updated from Git information by running 'python setup.py
# version'.
__version__ = '%s'
"""
        if not os.path.isdir(".git"):
            print "This does not appear to be a Git repository."
            return
        try:
            p = subprocess.Popen(["git", "describe",
                                  "--tags", "--always"],
                                 stdout=subprocess.PIPE)
        except EnvironmentError:
            print "unable to run git, leaving network-drift/_version.py alone"
            return
        stdout = p.communicate()[0]
        if p.returncode != 0:
            print "unable to run git, leaving network-drift/_version.py alone"
            return
        # our tags are like:  v2.2
        ver = stdout[len("v"):].strip()
        f = open("network-drift/_version.py", "w")
        f.write(VERSION % ver)
        f.close()
        print ("updated _version.py to '%s'" % ver)

setup(name="network-drift",
      version=get_version(),
      packages = find_packages(exclude=["simuPOP_examples","testdata","build","conf","doc","notebooks","simulations","test"]),
      scripts = [
          'models/network-structure-eval.py',
          'models/network-subpop-eval.py',
          'models/parameter-sweep-for-localization.py'
      ],
      author='Carl P. Lipo',
      author_email='clipo@binghamton.edu',
      url='https://github.com/clipo/network-drift',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: Apache Software License',
          'Operating System :: POSIX',
          'Programming Language :: Python :: 3.7',
          'Topic :: Scientific/Engineering',
      ],
      license="APACHE",
      cmdclass={ "version": VersionUpdate, "install": CustomInstallClass, "develop": CustomDevelopClass },
      )
