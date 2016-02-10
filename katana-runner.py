#!/usr/bin/env python

"""Convenience wrapper for running amplicon soft clipper directly from source tree."""


from katana.clipper import main
import sys

if __name__ == '__main__':
    main(sys.argv)
