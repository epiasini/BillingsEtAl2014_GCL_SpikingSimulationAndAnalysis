#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Simple script to update analysis result data format to include
explicit 'tau' coordinate. Of course, it depends on specific methods
being available in archive.ResultsArchive, that are going to be
removed before merging with the main branch of the project."""
from parameters import PSlice, ParameterSpace

space = ParameterSpace(PSlice(300), PSlice(6), PSlice(2.), PSlice(4), PSlice(5.), PSlice(.1, 1.0, .1), PSlice(-30., 5., 5.), PSlice(120), PSlice(30), PSlice(10, 80, 10), PSlice(10), PSlice(20), PSlice(200), PSlice(40), PSlice(0.), PSlice(0), PSlice(5.), PSlice(2.))

for k,point in enumerate(space.flat):
    result = point.results_arch._old_load_from_disk()
    if result:
        point.results_arch._write_from_memory()
        point.results_arch._delete_old_datasets()
    print "{0}: {1}/{2} {3}".format(result, k+1, len(space.flat), point.results_arch.path)

