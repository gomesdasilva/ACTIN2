import os
import subprocess

from setuptools import setup, find_packages

# try: from setuptools import setup
# except: from distutils.core import setup


path = os.path.dirname(os.path.realpath(__file__))
version_file = os.path.join(path, "actin2", "VERSION")

try:
    with open(version_file, 'r') as file:
        version = file.read()
except: version = "unknown"


setup(name = 'actin2',
      version = version,
      description = 'Activity Indices Calculator',
      url = 'http://github.com/gomesdasilva/ACTIN2',
      #download_url = 'https://github.com/gomesdasilva/ACTIN/archive/v1.3.2.tar.gz',
      author = 'Joao Gomes da Silva',
      author_email = 'Joao.Silva@astro.up.pt',
      license = 'MIT',
      keywords = ['astronomy', 'activity', 'fits', 'harps', 'harps-n', 'espresso', 'radial velocity', 'exoplanets', 'spirou'],
      #packages = ['actin2'],
      packages=find_packages(),
      package_data={
        '': ['spectrographs/*.py', 'test/*'],
    },
      #entry_points = {
      #  "console_scripts": ['actin = actin.actin:main']
      #  },
      include_package_data = True,
      install_requires = ['numpy', 'pandas', 'astropy', 'matplotlib', 'scipy'],
      zip_safe=False
      )


