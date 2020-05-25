'''
Created on 31/gen/2011

@author: Castrense Savojardo
'''

import numpy

class InvalidCheckpointFileError(Exception):
    def __init__(self):
        pass

def BlastCheckPointProfile(checkpointFile, newFormat = True):
    try:
        checkpointFile = open(checkpointFile).readlines()
    except IOError:
        raise
    profile = None
    if newFormat:
        try:
            profile = _profileParseNew(checkpointFile)
        except:
            raise
    return profile

def _profileParseNew(checkpoint, inAAOrder="ARNDCQEGHILKMFPSTWYV",
                     outAAOrder="VLIMFWYGAPSTCHRKQEND"):
    headerSize = 3
    footerSize = -6
    shift = 22
    aaOrder = checkpoint[2].split()[:20]
    try:
        _check(checkpoint[1])
    except:
        raise

    profile = []

    for line in checkpoint[headerSize:footerSize]:
        line = line.split()
        pos = numpy.zeros(20)
        for j in range(20):
            val = float(line[j + shift]) / 100.0
            if not val == 0.0:
                index_j = outAAOrder.index(inAAOrder[j])
                pos[index_j] = val
        if numpy.sum(pos) == 0.0:
            aa = line[1]
            try:
                # clipping the primary sequence
                pos[aaOrder.index(aa)] = 0.001
            except ValueError:
                pass
        profile.append(pos)
    return numpy.array(profile)

def _check(line):
    import re
    if not re.search('Last position-specific scoring matrix computed', line):
        raise InvalidCheckpointFileError
