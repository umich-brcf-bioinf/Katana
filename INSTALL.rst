Installing AmpliconSoftClipper
==================
AmpliconSoftClipper has been tested with Python 2.7 on OSX and \*nix.

Prerequisites
-------------
.. note:: Pip installs all required libraries; see [Installing] below.


* natsort (3.5.2)  
* nosetests, testfixtures required for running automated tests


Installing
----------
You can install from source from github:

``$ pip install git+https://github.com/umich-brcf-bioinf/AmpliconSoftClipper``

If you don't have root permissions, you can install locally:

``$ pip install git+https://github.com/umich-brcf-bioinf/AmpliconSoftClipper --user``

Following the pip install, you may need to adjust your path settings to include home/.local/bin. 


If you already have pysam installed, you can also clone from github and run directly from the source like so:

``$ git clone https://github.com/umich-brcf-bioinf/AmpliconSoftClipper``

``$ AmpliconSoftClipper/clipper-runner.py manifest.txt input.bam output.bam``
  
