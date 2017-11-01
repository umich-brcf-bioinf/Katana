Installing Katana
==================
Katana has been tested with:

* Python 2.7 and 3.4
* pysam 0.8 - 0.12
* OSX and \*nix

Prerequisites
-------------
.. note:: Pip installs all required libraries; see [Installing] below.


* natsort (3.5.2)  
* nosetests, testfixtures required for running automated tests


Installing
----------
The easiest way to install Katana is through PyPI. Get pip if it's
not available in your system:

::

   $ pip install katana

You can install from source from github:

``$ pip install git+https://github.com/umich-brcf-bioinf/Katana``

If you don't have root permissions, you can install locally:

``$ pip install git+https://github.com/umich-brcf-bioinf/Katana --user``

Following the pip install, you may need to adjust your path settings to include home/.local/bin. 
Then you can execute so:

``$ katana primers.txt input.bam output.bam``

If you already have prerequisite modules installed, you can also clone from github and run directly from the source like so:

``$ git clone https://github.com/umich-brcf-bioinf/Katana``

``$ katana/katana-runner.py primers.txt input.bam output.bam``

