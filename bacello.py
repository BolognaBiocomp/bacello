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

def run_pssm(ns):
    seq, seqid = utils.read1Fasta(ns.fasta)
    prof = cpparser.BlastCheckPointProfile(ns.pssm)
    loc, score = bacello(seq, prof, ns.kingdom)
    ofs = open(ns.outf, 'w')
    if ns.outfmt == "gff3":
        utils.write_gff_output(seqid, seq, ofs, loc, round(score,2))
    else:
        acc_json = utils.get_json_output(seqid, seq, loc, round(score,2))
        json.dump([acc_json], ofs, indent=5)
    ofs.close()
    sys.exit(0)

def run_multifasta(ns):

    data_cache = utils.get_data_cache(ns.cache_dir)
    i = 0
    ofs = open(ns.outf, 'w')
    protein_jsons = []
    for record in SeqIO.parse(ns.fasta, 'fasta'):
        workEnv = workenv.TemporaryEnv()
        prefix = "seq%s" % i
        fastaSeq  = workEnv.createFile(prefix+".", ".fasta")
        SeqIO.write([record], fastaSeq, 'fasta')
        seq, seqid = utils.read1Fasta(fastaSeq)
        pssmFile = blast.runPsiBlast(prefix, ns.dbfile, fastaSeq, workEnv, data_cache=data_cache,
                                     num_alignments=ns.pbnalign, num_iterations=ns.pbniter, evalue=ns.pbeval,
                                     threads=ns.threads)
        prof = cpparser.BlastCheckPointProfile(pssmFile)
        loc, score = bacello(seq, prof, ns.kingdom)
        if ns.outfmt == "gff3":
            utils.write_gff_output(seqid, seq, ofs, loc, round(score,2))
        else:
            acc_json = utils.get_json_output(seqid, seq, loc, round(score,2))
            protein_jsons.append(acc_json)
        workEnv.destroy()
    if ns.outfmt == "json":
        json.dump(protein_jsons, ofs, indent=5)
    ofs.close()
    sys.exit(0)

def main():
    DESC = "BaCelLo: Prediction of subcellular localization"
    parser = argparse.ArgumentParser(description = DESC, prog = "bacello.py")
    subparsers   = parser.add_subparsers(title = "subcommands",
                                         description = "valid subcommands",
                                         help = "additional help")
    multifasta  = subparsers.add_parser("multi-fasta",
                                          help = "Multi-FASTA input module",
                                          description = "BaCelLo: Multi-FASTA input module.")
    pssm  = subparsers.add_parser("pssm", help = "PSSM input module (one sequence at a time)",
                                    description = "BaCelLo: PSSM input module.")
    multifasta.add_argument("-f", "--fasta", help = "The input sequence in FASTA format", dest = "fasta", required = True)
    multifasta.add_argument("-d", "--dbfile", help = "The PSIBLAST DB file", dest = "dbfile", required= True)
    multifasta.add_argument("-k", "--kingdom", help = "The kingdom the sequence belongs to: P=Plant, A=Animal, F=Fungi", choices=['P', 'A', 'F'], dest = "kingdom", required = True)
    multifasta.add_argument("-o", "--output", help = "The output file name", dest="outf", required = "True")
    multifasta.add_argument("-m", "--outfmt", help = "The output format: json or gff3 (default)", choices=['json', 'gff3'], required = False, default = "gff3")
    multifasta.add_argument("-c", "--cache-dir", help="Cache dir for alignemnts", dest="cache_dir", required=False, default=None)
    multifasta.add_argument("-a", "--threads", help="Number of threads (default 1)", dest="threads", required=False, default=1, type=int)
    multifasta.add_argument("-j", "--psiblast-iter", help="Number of PSIBLAST iterations (default 3)", dest="pbniter", required=False, default=3, type=int)
    multifasta.add_argument("-n", "--psiblast-nalign", help="PSIBLAST num_alignments parameter (default 5000)", dest="pbnalign", required=False, default=5000, type=int)
    multifasta.add_argument("-e", "--psiblast-evalue", help="PSIBLAST evalue parameter (default 0.001)", dest="pbeval", required=False, default=0.001, type=float)
    multifasta.set_defaults(func=run_multifasta)
    pssm.add_argument("-f", "--fasta", help = "The input sequence in FASTA format", dest = "fasta", required = True)
    pssm.add_argument("-p", "--pssm", help = "PSSM file as generated by PSI-BLAST", dest = "pssm", required = True)
    pssm.add_argument("-k", "--kingdom", help = "The kingdom the sequence belongs to: P=Plant, A=Animal, F=Fungi", choices=['P', 'A', 'F'], dest = "kingdom", required = True)
    pssm.add_argument("-o", "--output", help = "The output file name", dest="outf", required = "True")
    pssm.add_argument("-m", "--outfmt", help = "The output format: json or gff3 (default)", choices=['json', 'gff3'], required = False, default = "gff3")
    pssm.set_defaults(func=run_pssm)
    if len(sys.argv) == 1:
      parser.print_help()
    else:
      ns = parser.parse_args()
      ns.func(ns)



if __name__=='__main__':
    main()
