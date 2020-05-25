#----------- ANDREA PIERLEONI--------------#
import sys
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
