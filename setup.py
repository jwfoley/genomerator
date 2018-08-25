from setuptools import setup, find_packages

# Read long description from readme.md.
try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst')
except:
    long_description = None

# Read version from package.
from ionize.__version__ import __version__

setup(name='ionize',
      version=__version__,
      author='Joseph W. Foley',
      author_email='jwfoley@gmail.com',
      url="https://github.com/jwfoley/genomerator"
      classifiers=[
          "Programming Language :: Python",
          "Development Status :: 3 - Alpha",
          "Intended Audience :: Science/Research",
          "Operating System :: OS Independent",
          "Topic :: Software Development :: Libraries :: Python Modules",
          ],
      use_2to3=False,
      license='LICENSE',
      description='',
      long_description=long_description,
      packages=find_packages(),
      requires=[],
      package_data=dict(),
      entry_points=dict(),
      test_suite="genomerator.tests",
      )
