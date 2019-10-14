# This file is part of xtb.
#
# Copyright (C) 2019 Sebastian Ehlert
#
# xtb is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# xtb is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with xtb.  If not, see <https://www.gnu.org/licenses/>.

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name="xtb",
      version="6.2.1",
      author="Sebastian Ehlert",
      author_email="ehlert@thch.uni-bonn.de",
      description="Wrapper for the extended tight binding program",
      long_description=long_description,
      long_description_content_type="text/markdown",
      keywords="tight-binding",
      url="https://github.com/grimme-lab/xtb",
      project_urls={
          "Documentation": "https://xtb-docs.readthedocs.io/en/latest/contents.html",
      },
      packages=find_packages(),
      install_requires=['ase'],
      classifiers=[
          "Programming Language :: Python :: 3",
          "License :: OSI Approved :: LGPL3",
          "Operating System :: Linux",
      ],
      python_requires=">=3.5",
)
