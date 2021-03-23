#!/usr/bin/env python
import os
BACELLO_HOME=os.environ.get("BACELLO_HOME")
import sys
sys.path.append(BACELLO_HOME)
import numpy
import argparse
import json
from bacellolib import utils
from bacellolib import svm as svmlib
from bacellolib import cpparser
from bacellolib import workenv
from bacellolib import blast
from bacellolib import config
from Bio import SeqIO

def bacello(sequenza, profilo, classe='A'):
    ''' versione light di BaCelLo:
        prende in input una sequenza aminoacidica il suo profilo e la classe di appartenza per il regno eucariotico:
        'A'===>Animali(Metazoa)
        'F'===>Funghi(Fungi)
        'P'===>Piante(Viridiplantae)

        da in output la localizzazione predetta'''
    if sequenza[0]=='M':#do not consider initial met
        sequenza=sequenza[1:]
    if classe=='P':
        Lmodelli=[os.path.join(BACELLO_HOME, 'MOD/M1'), os.path.join(BACELLO_HOME, 'MOD/M2'), os.path.join(BACELLO_HOME, 'MOD/M3'), os.path.join(BACELLO_HOME, 'MOD/M4')]
        Lcodifiche=[[0,20,40,60],[0,20,40,60,-20,-50,-100],[0],[0]]
        Lsoglia_med=[0.442, -0.201, 0.375, -0.327]
    elif (classe=='A' or classe=='F'):
        Lmodelli=[os.path.join(BACELLO_HOME, 'MOD/M1'), os.path.join(BACELLO_HOME, 'MOD/M2'), os.path.join(BACELLO_HOME, 'MOD/M3')]
        Lcodifiche=[[0,20,40,60],[0,20,40,60,-20,-50,-100],[0],[0]]
        if classe=='A':
            Lsoglia_med=[0.122, 0.467, -0.568]
        else:
            Lsoglia_med=[0.571, 0.344, -0.644]
    #Preparazione input (calcolo frequenze aminoacidiche di sequenza e profilo)
    input_completo=[]
    for codice in [0,20,40,60,-20,-50,-100]:
        input_completo.extend(utils.cal_freq_prof(profilo,codice,baco=1))
        input_completo.extend(utils.cal_freq_seq(sequenza,codice))
    Llimiti=[160,280,40,40]#lunghezza dei vettori di input per i vari modelli
    Linput=[]#contiene gli input dei vari modelli
    for k in range(len(Lmodelli)):
        Linput.append(numpy.array(input_completo[:Llimiti[k]]))
    #Predizione con Kernel RBF (Classe di piero)
    Lrisultati_svm=[]#contiene il risultato delle predizioni per ogni modello
    for k in range(len(Lmodelli)):
        svm=svmlib.getSVMLight(Lmodelli[k]+'_'+classe)
        svmout=svm.predict(Linput[k])+Lsoglia_med[k]
        Lrisultati_svm.append(svmout)#calcolo con kernel RBF + aggiunta soglia di bilanciamento
    #Scelta localizzazione predetta
    if classe=='P':
        if Lrisultati_svm[0]>=0:
            prediction='Secretory'
            score=Lrisultati_svm[0]
        else:
            if Lrisultati_svm[1]>=0:
                if Lrisultati_svm[3]>=0:
                    prediction='Mitochondrion'
                    score=Lrisultati_svm[3]
                else:
                    prediction='Chloroplast'
                    score=-1*Lrisultati_svm[3]
            else:
                if Lrisultati_svm[2]>=0:
                    prediction='Nucleus'
                    score=Lrisultati_svm[2]
                else:
                    prediction='Cytoplasm'
                    score=-1*Lrisultati_svm[2]
    else:#animali e funghi
        if Lrisultati_svm[0]>=0:
            prediction='Secretory'
            score=Lrisultati_svm[0]
        else:
            if Lrisultati_svm[1]>=0:
                prediction='Mitochondrion'
                score=Lrisultati_svm[1]
            else:
                if Lrisultati_svm[2]>=0:
                    prediction='Nucleus'
                    score=Lrisultati_svm[2]
                else:
                    prediction='Cytoplasm'
                    score=-1*Lrisultati_svm[2]
    return prediction, score

#-----------MAIN--------------#

def run_json(ns):
    workEnv = workenv.TemporaryEnv()
    i = 0
    ifs = open(ns.i_json)
    input_json = json.load(ifs)
    ifs.close()
    protein_jsons = []
    for i_json in input_json:
        prefix = "seq%s" % i
        seq = i_json['sequence']['sequence']
        seqid = i_json['accession']
        fastaSeq  = workEnv.createFile(prefix+".", ".fasta")
        fsofs=open(fastaSeq,'w')
        #SeqIO.write([fasta], fsofs, 'fasta')
        print(">%s" % seqid, file=fsofs)
        print(seq, file=fsofs)
        fsofs.close()

        pssmFile = blast.runPsiBlast(prefix, config.BLASTDB, fastaSeq, workEnv)
        prof = cpparser.BlastCheckPointProfile(pssmFile)
        loc, score = bacello(seq, prof, ns.kingdom)
        acc_json = utils.get_json_output(i_json, loc, round(score,2))
        protein_jsons.append(acc_json)
    ofs = open(ns.outf, 'w')
    json.dump(protein_jsons, ofs, indent=5)
    ofs.close()
    workEnv.destroy()
    sys.exit(0)

def main():
    DESC = "BaCelLo: Prediction of subcellular localization"
    parser = argparse.ArgumentParser(description = DESC, prog = "bacello.py")
    parser.add_argument("-i", "--i-json", help = "The input JSON file name", dest = "i_json", required = True)
    parser.add_argument("-o", "--outf", help = "The output file", dest = "outf", required = True)
    parser.add_argument("-k", "--kingdom", help = "The kingdom the sequence belongs to: P=Plant, A=Animal, F=Fungi", choices=['P', 'A', 'F'], dest = "kingdom", required = True)
    ns = parser.parse_args()
    run_json(ns)



if __name__=='__main__':
    main()
