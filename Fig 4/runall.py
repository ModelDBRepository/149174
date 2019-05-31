#!/usr/bin/python

import os, sys

if os.name == 'nt':
    # Whadya know, windows still can't multitask properly
    try:
        from win32process import SetPriorityClass, GetCurrentProcess, \
             HIGH_PRIORITY_CLASS, NORMAL_PRIORITY_CLASS, \
             REALTIME_PRIORITY_CLASS, IDLE_PRIORITY_CLASS
        # Drop priority so that we can still use the computer
        SetPriorityClass(GetCurrentProcess(), IDLE_PRIORITY_CLASS)
    except:
        pass


def submit(c):
    print c
    os.system(c)

def IsFileDone(freq, drug, concentration, Eleak):
    fn = 'Data/M2I10sp_%d_%s_%g_%g.txt' % (freq, drug, concentration, Eleak)
    if not os.path.exists(fn):
        return False
    fh = open(fn)
    lines = fh.readlines()
    line = lines[-1]
    words = line.split()
    try:
        tm = float(words[0])
    except ValueError:
        return False
    if tm > 190:
        return True
    else:
        return False

drugnames = ['none', 'phenytoin', 'carbamazepine']
nrniv = '/Volumes/NEURON-7.2/NEURON-7.2/nrn/i386/bin/nrniv -dll /Users/dugsy/Dropbox/Projects/NaDrugs/Async/i386/.libs/libnrnmech.so'
for freq in [10, 20, 50]:
    for Eleak in [-70, -60, -50]:
        for drugcombo in [[0, 0], [1, 100], [2, 200]]:
            drugopt = drugcombo[0]
            concentration = drugcombo[1]
            drug = drugnames[drugopt]
            if not IsFileDone(freq, drug, concentration, Eleak):
                c = nrniv + ' -nobanner -c freq=%d -c drugopt=%d -c concentration=%g -c Eleak=%g async.hoc' % \
                (freq, drugopt, concentration, Eleak)
                submit(c)

