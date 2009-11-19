"""
:NAME:
motifparser

:DESC:
parser for MEME info-gibbs Gibbs MotifSampler ...


:HELP:
Format  [seqid, strand, position, word] where position starts form 0
Example ['1', '+', 12, 'ATTGT']

:LIMITATION:
seqid must be an integer

"""

import sys
import os
import re

import dna


TYPE = {
'info-gibbs'    : '.info',
'GAME'          : '.game',
'MEME'          : '.meme',
'AlignAce'      : '.alignace',
'BioProspector' : '.bioprospector',
'Gibbs'         : '.gibbs',
'MotifSampler'  : '.motifsampler',
'words'         : '.words',
}

EXT = dict([(TYPE[name], name) for name in TYPE])

########################################
#                                      #
# PARSERS
#                                      #
########################################
def info_gibbs(data):
    predicted_motif = re.findall(r'; (\d+)\t([\+\-])\t(\d+)\t([ACGT]+)', data)
    return [(p[0], p[1], int(p[2]), p[3]) for p in predicted_motif]


def words(data):
    predicted_motif = re.findall(r'([ACGT]+)', data)
    return [(0, '+', 0, p) for p in predicted_motif]


def GAME(data):
    predicted_motif = re.findall(r'>([\w\d]+)\s+([+-])(\d+)\n([ACGTNacgtn]+)', data)
    return [(p[0], p[1], int(p[2])-1, p[3]) for p in predicted_motif]


def MEME(data):
    BLOKS = re.findall(r'^BL.+?$(.*?)^//$', data,  re.MULTILINE | re.DOTALL)
    if len(BLOKS) == 0:
        return []

    predicted_motif = re.findall(r'([\d\w]+)\s+\(\s+(\d+)\)\s([ACGT]+)\s+(\d)', BLOKS[0])

    motif = []
    for m in predicted_motif:
        if m[3] == '1':
            strand = '+'
        else:
            strand = '-' 
        motif += [ (m[0], strand, int(m[1])-1, m[2]) ]
    return motif


def MotifSampler(data):
    motifs = []
    for motifstr in data.split('#'):
        if not motifstr.startswith('id'):
            continue
        try:
            #ic  = float(re.findall(r'ic: (\d+\.\d+)', motifstr)[0])
            llr = float(re.findall(r'll: (\d+\.\d+)', motifstr)[0])
        except IndexError: #-inf case
            continue
            
        motif = re.findall(r'([\d\w]+)\s+MotifSampler.*?(\d+).*([+-]).*site "([ACGT]+)"', motifstr)

        for i in range(len(motif)):
            site = list(motif[i])
            site[1], site[2] = site[2], site[1]
            site[2] = int(site[2])-1
            motif[i] = tuple(site)
        motifs += [(llr, motif)]

    motifs.sort()
    if len(motifs) == 0:
        return []
    else:
        return motifs[-1][1]


def Gibbs(data):
    MOTIF = re.findall(r'MAP MAXIMIZATION RESULTS.*^Num Motif(.*?)^Column', data,  re.MULTILINE | re.DOTALL)
    if len(MOTIF) == 0:
        return []

    predicted_motif = re.findall(r'(\d+),\s+(\d+)\s+(\d+)\s+[acgt]*\s+([ACGT]+).*?([FR])', MOTIF[0])
    motif = []
    for m in predicted_motif:
        if m[4] == 'F':
            strand = '+'
        else:
            strand = '-'

        motif += [ (m[0], strand, int(m[2])-1, m[3]) ]
    return motif


def AlignAce(data):
    MOTIF = re.findall(r'^Motif 1$(.*?)^MAP.+$', data,  re.MULTILINE | re.DOTALL)
    if len(MOTIF) == 0:
        return []
    predicted_motif = re.findall(r'([ACGT]+)\s+(\d+)\s+(\d+)\s+(\d+)', MOTIF[0])

    motif = []
    for m in predicted_motif:
        if m[3] == '1':
            strand = '+'
            word = m[0]
        else:
            strand = '-'
            word =  dna.reverse_complement(m[0])

        motif += [ (str(int(m[1])+1), '+', int(m[2]), word) ]
    return motif


def BioProspector(data):
    predicted_motif = re.findall(r'>([\d\w]+)\s+len \d+\s+site\s+#\d+\s+([rf])\s+(\d+)\s+([ACGT]+)', data)
    motif = []
    for m in predicted_motif:
        if m[1] == 'f':
            strand = '+'
        else:
            strand = '-'

        motif += [ (m[0], strand, int(m[2])-1, m[3]) ]
    return motif


def known_sites(data):
    known_motif = re.findall(r'\[(\d+)\t([+-])\t(\d+)\t([ACGT]+)\]', data)
    return [(p[0], p[1], int(p[2]), p[3]) for p in known_motif]


def parse_words(filename, format='auto'):
    '''return list of predicted sites (words only)
    optional format -- string BioProspector, GAME, MEME
    '''
    return [site[3] for site in parse(filename, format)]


def parse(filename, format='auto'):
    '''return list of predicted sites
    optional format -- string BioProspector, GAME, MEME
    '''
    if format == 'auto':
        format = EXT.get('.'+filename.split('.')[-1], '')

    if format not in TYPE.keys():
        sys.stderr.write('WARNING: unknonw file extension or unsupported format\n')
        return []

    format = format.replace('-', '_')    
    converter = globals()[format]
    data = open(filename).read()
    predicted_motif = converter(data)
    return predicted_motif


########################################
#                                      #
# STATS
#                                      #
########################################
def true_positive(known_motif, predicted_motif, distance=0):
    #if distance == 0:
    #    return len(predicted_motif.intersection(known_motif))

    k = set([ (int(i[0]), int(i[2])) for i in known_motif ])
    p = set([ (int(i[0]), int(i[2])) for i in predicted_motif ])

    TP = 0
    for s,p in p:
        for d in range(distance+1):
            if (s,p+d) in k or (s,p-d) in k:
                TP += 1
                continue
    return TP


def stats(knownfile, predictedfile, format, distance):

    # read known motifs
    f = open(knownfile)
    known_motif = known_sites(f.read())
    f.close()

    # read predicted motifs
    format = format.replace('-', '_')
    converter = globals()[format]
    predicteddata = open(predictedfile).read()
    predicted_motif = converter(predicteddata)

    TP = true_positive(known_motif, predicted_motif, distance)
    predicted = len(predicted_motif)
    real = len(known_motif)
    
    if  predicted == 0:
        PPV = 0.0
    else:
        PPV = TP / float(predicted)
    Sn = TP / float(real)
    PC = TP / float(predicted + real - TP)

    if PPV == 0 and Sn == 0:
        F = 0.0
    else:
        F = 2 * PPV * Sn / (PPV + Sn)

    R = {'real' : real, 'predicted' : predicted, 'TP' : TP, 'PPV' : PPV, 'Sn' : Sn, 'F': F, 'PC': PC}
    return R




