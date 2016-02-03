#pylint: disable=line-too-long
from setuptools import setup
from setuptools import find_packages
import os
import ampliconsoftclipper

def read(*paths):
    """Build a file path from *paths* and return the contents."""
    with open(os.path.join(*paths), 'r') as filename:
        return filename.read()

setup(name='AmpliconSoftClipper',
      version = ampliconsoftclipper.__version__,
      description=('Command-line tool to soft-clip reads based on primer locations.'),
      long_description=(read('README.rst') + '\n\n' +
                        read('CHANGELOG.rst') + '\n\n' +
                        read('AUTHORS.rst')),
      url='https://github.com/umich-brcf-bioinf/AmpliconSoftClipper',
      author='University of Michigan Bioinformatics Core',
      author_email='bfx-ampliconsoftclipper@umich.edu',
      license='Apache',
      packages=find_packages(exclude=['test*']),
      classifiers=['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: Apache Software License',
                   'Operating System :: Unix',
                   'Programming Language :: Python :: 2.7',
                   'Topic :: Scientific/Engineering :: Bio-Informatics'],
      keywords='bioinformatic exome-seq DNA-seq BAM',
      install_requires=['pysam', 'natsort'],
      entry_points={'console_scripts': ['clipper=ampliconsoftclipper.clipper:main']},
      test_suite='nose.collector',
      tests_require=['nose', 'pysam', 'natsort','testfixtures'],
      zip_safe=False)
