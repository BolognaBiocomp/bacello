#----------- ANDREA PIERLEONI--------------#
import sys
from . import config as cfg
import numpy

def read1Fasta(fname):
    try:
        lines=open(fname,'r').readlines()
    except:
        sys.stderr.write('cannot open file '+fname+'\n')
        sys.exit()
    seq=''
    seqid=""
    for line in lines:
        if line[0]!='>':
            seq+=''.join(line.split())
        else:
            seqid=line.split()[0][1:]
    return seq, seqid

def cal_freq_seq(sequenza, fin=0):
    Laa = ['K','R','H',   'D','E',   'N','Q','S','T','Y',   'A','V','L','I','P','F','M','W','G','C']
    if fin > 0:
        sequenza=sequenza[:fin]
    elif fin < 0:
        sequenza=sequenza[fin:]
    else:
        pass
    Lrisultato=[]
    tot=float(len(sequenza))
    for aa in Laa:
        Lrisultato.append(sequenza.count(aa)/tot)
    return Lrisultato

def one_hot_encoding(sequence, alph="VLIMFWYGAPSTCHRKQEND"):
    profile = numpy.zeros((len(sequence), 20))
    for (i, aa) in enumerate(sequence):
        try:
            j = alph.index(aa)
            profile[i,j] = 1.0
        except:
            pass
    return profile

def cal_freq_prof(profilo, fin=0,baco=0):
    '''calcola la frequenza (come vettore da 20 aminoacidi) di un profilo
      modificare questa funzione se l'input del profilo e' in un formato diverso
       FORMATO PROFILO ATTUALE:
    array numpy:
    primo livello: posizioni nella sequenza dall Nter al Cter
    secondo livello: frequenza aminocidica nella posizione, secondo l'ordine stabilito da Laa
    fin:    definisce la porzione di sequenza di cui calcolare la frequenza aminoacidica
        0---> sequenza intera
        numero positivo intero---> da N-ternimale a posizione specificata
        numero negativo intero--->da C-ternimale a posizione specificata

    baco= implementa il calcolo della frequenza usato nella prima versione di bacello'''
    tot=float(len(profilo))
    if fin > 0:
        profilo=profilo[:fin]
    elif fin < 0:
        profilo=profilo[fin:]
    else:
        pass
    Lrisultato=[]
    if baco:
        Lrisultato.extend(profilo.sum(0))
        for i in range(len(Lrisultato)):
            Lrisultato[i]=Lrisultato[i]/tot
    else:
        Lrisultato.extend(profilo.mean(0))
    return Lrisultato

def get_data_cache(cache_dir):
    import os
    from . import datacache
    ret = None
    if cache_dir is not None:
        if os.path.isdir(cache_dir):
            ret = datacache.DataCache(cache_dir)
    return ret

def write_gff_output(acc, sequence, output_file, localization, prob):
    l = len(sequence)
    print(acc, "BaCelLo", cfg.locmap[localization][0], 1, l, prob, ".", ".",
    "Ontology_term=%s;evidence=ECO:0000256" % cfg.locmap[localization][1],
    file = output_file, sep = "\t")

def get_json_output(acc, sequence, localization, prob):
    acc_json = {'accession': acc, 'dbReferences': [], 'comments': []}
    acc_json['sequence'] = {
                              "length": len(sequence),
                              "sequence": sequence
                           }
    c = cfg.locmap[localization]
    go_info = cfg.GOINFO[c[1]]
    acc_json['dbReferences'].append({
        "id": c[1],
        "type": "GO",
        "properties": {
          "term": go_info['GO']['properties']['term'],
          "source": "IEA:BaCelLo",
          "score": round(float(prob),2)
        },
        "evidences": [
          {
            "code": "ECO:0000256",
            "source": {
              "name": "SAM",
              "id": "BaCelLo",
              "url": "http://gpcr.biocomp.unibo.it/bacello/",
            }
          }
        ]
    })
    acc_json['comments'].append({
        "type": "SUBCELLULAR_LOCATION",
        "locations": [
          {
            "location": {
              "value": go_info["uniprot"]["location"]["value"],
              "score": round(float(prob),2),
              "evidences": [
                {
                  "code": "ECO:0000256",
                  "source": {
                    "name": "SAM",
                    "id": "BaCelLo",
                    "url": "http://gpcr.biocomp.unibo.it/bacello/",
                  }
                }
              ]
            }
          }
        ]
    })
    return acc_json
