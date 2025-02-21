
from itertools import combinations

import os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import random
import numpy as np
import scipy
from scipy import stats
from collections import Counter
import datetime
from sklearn.neighbors import KernelDensity
from sklearn.cluster import HDBSCAN
#from hdbscan import HDBSCAN

from sklearn.preprocessing import normalize
from itertools import product
import umap.umap_ as umap
import subprocess
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
from sklearn.decomposition import PCA
from sklearn import svm
import joblib

#dependencies = > porechop, seqtk, mafft, fastp,


def chkpol(liste,sz,strand):
    ls = []
    for x in liste:
        tag = 0
        id = x[0]
        if strand == 'r':
            read = "".join(base for base in reversed(x[1]))
            qual = "".join(base for base in reversed(x[2]))
        else:
            read = x[1]
            qual = x[2]
        if 'A' * 10 in read[:100]: ind, let, tag = read[:100].index('A' * 10), 'A', 1
        elif 'C' * 10 in read[:100]: ind, let, tag = read[:100].index('C' * 10), 'C', 1
        elif 'G' * 10 in read[:100]: ind, let, tag = read[:100].index('G' * 10), 'G', 1
        elif 'T' * 10 in read[:100]: ind, let, tag = read[:100].index('T' * 10), 'T', 1
        elif 'AC' * 20 in read[:100]: ind, let, tag = read[:100].index('AC' * 20), 'AC', 2
        elif 'AG' * 20 in read[:100]:ind, let, tag = read[:100].index('AG' * 20), 'AG', 2
        elif 'AT' * 20 in read[:100]:ind, let, tag = read[:100].index('AT' * 20), 'AT', 2
        elif 'CG' * 20 in read[:100]:ind, let, tag = read[:100].index('CG' * 20), 'CG', 2
        elif 'CT' * 20 in read[:100]:ind, let, tag = read[:100].index('CT' * 20), 'CT', 2
        elif 'GT' * 20 in read[:100]:ind, let, tag = read[:100].index('GT' * 20), 'GT', 2

        if tag == 1:
            while ind <= len(read):
                try:
                    if read[ind] != let:
                        read = read[ind:]
                        qual = qual[ind:]
                        ind = len(read)
                except:
                    pass
                ind = ind + 1
        elif tag == 2:
            while ind <= len(read):
                try:
                    if read[ind:ind+2] != let:
                        read = read[ind+1:]
                        qual = qual[ind+1:]
                        ind = len(read)
                except:
                    pass
                ind = ind + 2
        if len(read) >= sz:
            if strand == 'r':
                ls.append([id, "".join(base for base in reversed(read)), "".join(base for base in reversed(qual))])
            else:
                ls.append([id,read ,qual])
    return ls

def abnPoly(minsize):
    with open("temp_/trimmed_filtered_fastp.fastq") as f:
        rds = f.readlines()
    i = 0
    ls = []
    while i + 2 <= len(rds):
        sls = []
        sls.append(rds[i].replace('\n',''))
        sls.append(rds[i+1].replace('\n',''))
        sls.append(rds[i+3].replace('\n',''))
        ls.append(sls)
        i = i + 4
    fout = open("temp_/trimmed_filtered.fastq",'w')
    nls = chkpol(ls,minsize,'')
    nls = chkpol(nls, minsize, 'r')
    for x in nls:
        fout.write(x[0] + '\n' + x[1]  + '\n+\n' + x[2] + '\n')
    fout.close()


def cleanT(mins,maxs,loc,fastq,locus):   #min_size, max_size, location folder, fastqfilename,locus name = ITS or EF1
    with open(fastq) as f:
        dicParam['Nb. of ON reads before cleaning :\t\t'] = str(int(len(f.readlines())/4))
    os.system('mkdir temp_')
    os.system('./Porechop/porechop-runner.py -t 12 -i ' + fastq + ' --check_reads 1000 --discard_middle -o temp_/trimmed.fastq -v 1')    # trim barcodes and adapters with porechop
    #os.system('porechop -t 12 -i ' + fastq + ' --check_reads 1000 --discard_middle -o temp_/trimmed.fastq -v 1')  # trim barcodes and adapters with porechop
    with open('temp_/trimmed.fastq') as f:
        dicParam['porechop -t 12 -i ' + fastq + ' --check_reads 1000 --discard_middle -o temp_/trimmed.fastq -v 1' + '\t' + 'Nb. reads after porechop :\t'] = str(int(len(f.readlines())/4))

    try:
        os.system('./fastp -i temp_/trimmed.fastq -q 8 -l ' + str(mins) + ' --length_limit ' + str(maxs) + ' -y -g -x -o temp_/trimmed_filtered_fastp.fastq -j temp_/trimmed_filtered_fastp.out.json -h temp_/trimmed_filtered_fastp.out.html')  #run fastp
    except:
        os.system('cp temp_/trimmed.fastq temp_/trimmed_filtered_fastp.fastq')

    os.system('rm temp_/trimmed.fastq')


    os.system('./seqtk/seqtk seq -a temp_/trimmed_filtered_fastp.fastq > temp_/trimmed_filtered_fastp.fas')
    with open('temp_/trimmed_filtered_fastp.fastq') as f:
        dicParam['./fastp -i temp_/trimmed.fastq -q 8 -l ' + str(mins) + ' --length_limit ' + str(maxs) + ' -y -g -x -o temp_/trimmed_filtered_fastp.fastq -j temp_/trimmed_filtered_fastp.out.json -h temp_/trimmed_filtered_fastp.out.html' + '\t' + 'Nb. reads after fastp :\t'] = str(int(len(f.readlines())/4))
    abnPoly(mins)
    with open('temp_/trimmed_filtered.fastq') as f:
        dicParam['Trim PolyX - abnPoly(' + str(mins) +')\t' + 'Nb. reads after trimming polyX :\t'] = str(int(len(f.readlines())/4))
    os.system('./seqtk/seqtk seq -a temp_/trimmed_filtered.fastq > temp_/trimmed_filtered.fas')

    print(dicParam)
    cleanFastq()


def revcomp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rx = "".join(complement.get(base, base) for base in reversed(seq))
    return rx


def create_kmers():
    # 1 - generate a canonical list of 5-mers
    ls = [''.join(c) for c in product('ACGT', repeat=5)]
    lf = []
    for x in ls:
        rev_x = revcomp(x)
        if x < rev_x:
            if x not in lf:
                lf.append(x)
        else:
            if rev_x not in lf:
                lf.append(rev_x)
    return lf

def count_kmers(lsk,locus):
    # 2 - get the kmer count for each read and put it in a list
    if locus == 'ITS':
        clf = joblib.load("modelreads")  # loading model that discreminate between fungal ITS and contaminent reads

    else:
        clf = joblib.load("modelreads_EF")  #loading model that discreminate between fungal ITS and contaminent reads
    frej = open('temp_/Conta_reads.fas','w')
    listeOfk = [] #final list that contains the kmer freq for each read in the library
    with open('temp_/trimmed_filtered.fas') as f:
        ls_rds = (f.read().split('>'))
    ls_name = []
    dicseq = {} #dictionnary with the original sequences from the fasta
    for read in ls_rds[1:]:
        _kmers_ = {}
        for x in lsk:
            _kmers_[x] = 0
        kmer = []
        i = 0
        read_name = read.split(' ',1)[0]
        #ls_name.append(read_name)
        seq = read.split('\n',1)[1].replace('\n','')
        dicseq[read_name]  = seq
        while i < len(seq) - 5 :
            km = seq[i:i+5]
            try:
                _kmers_[km] = _kmers_[km] + 1
            except:
                km = revcomp(km)
                _kmers_[km] = _kmers_[km] + 1
            i = i + 1
        for x in _kmers_:
            kmer.append(_kmers_[x])
        y_pred = clf.predict([kmer])
        if y_pred[0] == 'F':
                listeOfk.append(kmer)
                ls_name.append(read_name)
        else:
                frej.write('>' + read_name + '\n' + dicseq[read_name] + '\n')

        frej.close()
        dicParam['No. of fungal reads :\t'] = str(len(ls_name))
        return ls_name, listeOfk, dicseq


def runHDBSCAN(kmers,nn,md,mcs):
    print("******", nn,md,mcs)
    l2D = np.array(kmers)
    l2DN = normalize(l2D, axis=1, norm='l1')   #produces normalized frequencies for kmers
    print(sum(l2DN[2]))
    X_embedded = umap.UMAP(n_neighbors=int(nn), min_dist=float(md), verbose=2).fit_transform(l2DN)     #default parameters 15, 0.1
    #l2DN = PCA(n_components=2).fit_transform(l2DN)
    clusterer = HDBSCAN(min_cluster_size=int(mcs),allow_single_cluster=False)  #initial cluster size of 5 worked
    cluster_labels = clusterer.fit_predict(X_embedded)
    cluster_labels = cluster_labels.tolist()
    cluster_labels = map(str, cluster_labels)
    cluster_labels = list(map(lambda x: x.replace('-1', 'out'), cluster_labels))
    return cluster_labels

def build_clusters(nm,cl,dseq,_inpt_):
    dcl = {}  #dseq is a dic that contains the original sequences
    dicPam, dicSeq = {}, {}
    for x in cl:
        if x not in dcl:
            dcl[x] = []
    print(dcl)
    print(len(dcl) -1, 'clusters')
    final = open('temp_/polished_seq.fas','w')
    for x,y in zip(nm,cl):
        dcl[y].append(x)
    for x in dcl:
        print(x)
        flistreads = open("reads.lst","w")
        fout = open('temp_/cluster_' + str(x) + '.fas', 'w')
        for id in dcl[x]:
            fout.write('>' + id + '\n' + dseq[id] + '\n')
            flistreads.write(id + '\n')
        fout.close()
        flistreads.close()
        dicPam[str(x)] = [str(len(dcl[x]))]   #dicPam contains all the parameters of the clusters starting with the number of  reads
        dicSeq[str(x)] = []
        #if str(x) != "-1":
        os.system("./seqtk/seqtk subseq temp_/trimmed_filtered_c.fastq reads.lst > " + 'temp_/cluster_' + str(x) + ".fq")  #retrieve fastq reads for cluster x
        os.system("mafft --adjustdirection temp_/cluster_" + str(x) + ".fas > temp_/cluster_" + str(x) + "_aln.fas")  #align fasta reads for cluster x
        #reads correction
        #os.system("./canu2.2/bin/canu -correct -p cluster_" + str(x) + " -d Clusters -nanopore cluster_" + str(x) + ".fq")
        _cons_, longeur = consensus("temp_/cluster_" + str(x) + "_aln.fas", str(x))
        print(_cons_)
        polseq = polish(_cons_,'temp_/cluster_' + str(x) + ".fq",str(x))
        final.write(polseq)
        dicPam[str(x)].append(str(longeur))
        try:
            dicSeq[str(x)].append(polseq.split('\n',1)[1].replace('\n',''))
        except:
            dicSeq[str(x)].append('')
    final.close()
    return dicPam, dicSeq


def polish(cs,fq,nb):
    fx = open('inp.fasta','w')
    fx.write(cs.replace('-',''))
    fx.close()
    n = 0
    while n < 5:
        if n == 0:
            feed = 'inp.fasta'
        print(n, feed)
        os.system('minimap2/minimap2 -ax map-ont ' + feed + ' ' + fq + ' -o mapping.sam --secondary=no')
        os.system('./racon/build/bin/racon -m 8 -x -6 -g -8 -w 200 -t 14 ' + fq + ' mapping.sam ' + feed +  ' > racon' + str(n) + '.fasta')
        #os.system('rm mapping.sam')
        feed = 'racon' + str(n) + '.fasta'
        n = n + 1
    with open('racon4.fasta') as f:
        ps = f.read()
    os.system("mv racon4.fasta temp_/cluster_" + nb + "_polished.fasta")
    os.system('rm racon*.fasta')
    return ps



def cleanFastq():
    out = open("temp_/trimmed_filtered_c.fastq","w")
    with open('temp_/trimmed_filtered.fastq') as f:
        ls = f.readlines()
    for x in ls:
        if x[0] == '@':
            out.write(x.split(' ')[0] + '\n')
        else:
            out.write(x)
    out.close()

def Most_Common(lst):
    data = Counter(lst)
    return data.most_common(1)[0][0]

def consensus(alig,nb):
    fx= open(alig)
    ct = fx.read().split('>')
    fx.close()
    ls = []
    for x in ct[1:]:
        ls.append(x.split('\n',1)[1].replace('\n',''))
    lss = [list(i) for i in zip(*ls)]
    cons = ''
    for ssls in lss:
        cons = cons + Most_Common(ssls)
    fx = open(alig,'a')
    fx.write('>consensus\n' + cons + '\n')
    fx.close()
    lg= len(cons.replace('-',''))
    cons = ">cluster_" + nb + '\n' + cons + '\n'
    return cons, lg


def annotate(dico,locus):
    if locus == 'ITS':
        #os.system('touch temp_/bnres.txt')
        #os.system('touch temp_/bnres.txt')
        os.system('ncbi-blast-2.9.0+/bin/blastn -db SeqDB/AllFungi_ITS.fasta -query temp_/polished_seq.fas -outfmt 6 -num_threads 24 -max_target_seqs 3 > temp_/bnres.txt')
        os.system('ncbi-blast-2.9.0+/bin/blastn -db SeqDB/ITS_RefSeq.fasta -query temp_/polished_seq.fas -outfmt 6 -num_threads 24 -max_target_seqs 3 > temp_/refSres.txt')
        #os.system('ncbi-blast-2.9.0+/bin/blastn -db SeqDB/AllFungi_ITS_S.fasta -query temp_/polished_seq.fas -outfmt 6 -num_threads 24 -max_target_seqs 3 > temp_/bnres.txt')
        #os.system('ncbi-blast-2.9.0+/bin/blastn -db SeqDB/ITS_RefSeq_S.fasta -query temp_/polished_seq.fas -outfmt 6 -num_threads 24 -max_target_seqs 3 > temp_/refSres.txt')


        dico = parseB(dico,'temp_/bnres.txt')
        dico = parseB(dico, 'temp_/refSres.txt')
    else:
        os.system('ncbi-blast-2.9.0+/bin/blastn -db SeqDB/EF1_ALLfungi.fas -query temp_/polished_seq.fas -outfmt 6 -num_threads 24 -max_target_seqs 3 > temp_/bnres.txt')
        dico = parseB(dico, 'temp_/bnres.txt')
    return dico

def parseB(dico, infile):
    with open(infile) as f:
        lsb = f.readlines()
    dicdone = {}
    for line in lsb:
        if line.split('\t')[1] not in dicdone:
            sp = subprocess.check_output("efetch -db nucleotide -id " + line.split('\t')[1] + " -format gbc | xtract -pattern INSDSeq -element INSDSeq_source -element INSDSeq_taxonomy",shell = True)
            sp = sp.decode("utf-8")
            #print(sp)
            dicdone[line.split('\t')[1]] = (sp.replace('\n',''))
        else:
            sp = dicdone[line.split('\t')[1]]
        dico[line.split('\t')[0].replace('cluster_','')].append((line.split('\t')[1] + '\t'+ sp.replace('\n','') + '\t'+ line.split('\t')[2] +'\t'+ line.split('\t')[3] +'\t'+ line.split('\t')[10]))
    return dico


def writout(dico,dicSEQ,_loc_,_out_,cmdl):
    print("***",dicSEQ,"******")
    final = open(_out_,'w')
    final.write(cmdl + '\n\n')
    final.write(out + '\t\tBlastn on nt fungi\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tBlastn on Refseq fungi\n')
    tg = '\ttaxon\ttaxonomy\tsimilarity\toverlap\te-value'
    final.write('Cluster #\tConsensus Seq.\tNb. of reads\tCluster size'+ tg * 6 + '\n')
    for x in dico:
        #if x != '-1':
            print(dico[x])
            final.write('cluster_' + x + '\t' + dicSEQ[str(x)][0] + '\t' + '\t'.join(dico[x]) + '\n') #clusterID, nb of reads, cluster size and blast annot
#        except:
#            out = 'cluster_' + x + '\t' + '\t'.join(dico[x][0]) + '\n'  # clusterID, nb of reads, cluster size and blast annot
    final.write("\n")
    for x in dicParam:
        final.write(x + dicParam[x] + '\n')
    final.close()
    os.system('rm -r ' + _loc_ + '/temp_')
    os.system('mv temp_ ' + _loc_ + '/.')

################################################## MAIN #######################################################
#usage ==> ~/miniconda3/bin/python3  409_40790_1/barcode01 ITS #ITS or EF1 will trigger the checklocus function


try:
    cmdline = sys.argv[0] + ' ' + sys.argv[1] + ' ' + sys.argv[2] + ' ' + sys.argv[3] #+ ' ' + sys.argv[4]
except:
    cmdline = sys.argv[0] + ' ' + sys.argv[1] + ' ' + sys.argv[2]
location = sys.argv[1]
os.system('rm ' + sys.argv[1] + '/' + sys.argv[1].split('/')[-1] + '.fastq')
os.system('cat ' + sys.argv[1] + '/*.fastq > ' + sys.argv[1] + '/' +  sys.argv[1].split('/')[-1] + '.fastq')
inpt = './' + sys.argv[1] + '/' + sys.argv[1].split('/')[-1] + '.fastq'

try:
    min_size = int(sys.argv[3])
except:
    min_size = 250
    pass
max_size = 1000 #int(sys.argv[4])

try:
    min_size = int(sys.argv[3])
except:
    min_size = 250
    pass
max_size = 1000 #int(sys.argv[4])

try:
    _nn_,_md_,_mcs_ = sys.argv[4], sys.argv[5], sys.argv[6]
except:

    _nn_, _md_, _mcs_ = 15, 0.1, 20     #try 15 0.01 5 for detection

print(str(_nn_), str(_md_), str(_mcs_))
dicParam = {'n_neighbors':str(_nn_), 'min_dist':str(_md_), 'min_cluster_size':str(_mcs_)}


barcode = location.split('/')[-1]
out = location + '/' + barcode + '_' + dicParam['n_neighbors'] + '_' + dicParam['min_dist'] + '_' + dicParam['min_cluster_size'] + '_report.xls'

cleanT(min_size, max_size, location, inpt, sys.argv[2])
tk = create_kmers()
rNames, kmerFreq, lstseq = count_kmers(tk,sys.argv[2])
clusterLs = runHDBSCAN(kmerFreq,_nn_,_md_,_mcs_)
dicP,dicS = build_clusters(rNames,clusterLs,lstseq,inpt)  #dicP contains the parameters of the cluster i.e. here after build_clusters its {cluster_number : [nb of reads, lenght of the cluster]}
print(dicP)
dicP = annotate(dicP,sys.argv[2])
writout(dicP,dicS,location,out,cmdline)
