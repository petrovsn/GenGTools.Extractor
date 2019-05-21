from sklearn import model_selection as cross_validation
from sklearn import cross_validation, datasets, metrics, tree 
from sklearn import ensemble,  learning_curve
import sklearn

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt 
import copy
import re
import math
from itertools import product

def joinStrings(stringList):
    return ''.join(string for string in stringList)
def get_all_substrings(input_string, l):
    length = len(input_string)-l
    return [input_string[i:i+l] for i in xrange(length+1)]
def getfeautures(read):
    #print read
    rep = 2
    fs = {joinStrings(s):0 for s in product('ATCG', repeat=rep)}
    
    for sr in get_all_substrings(read,rep):
        if ('N' not in sr and 
        'W' not in sr and
        'Y' not in sr and
        'S' not in sr and
        'M' not in sr and
        'R' not in sr):
            fs[sr] = fs[sr] + 1
    for k in fs:
        fs[k] = fs[k]*1.0/len(read)
    return fs.values()


import re
def getcoords(header):
    T = re.findall(r'\d{1,}\.\.\d{1,}',  header)
    return [(int(t[0]),int(t[1])) for t in [q.split('..') for q in T]]


def getCDS(file):
    return loadfasta(file)


def getClobCoords(CDSdict):
    Res = []
    for k in CDSdict.keys():
        Res = Res+getcoords(k)
    return Res

def getNonCDS(CDSdict, ref):
    NonCDS = {}
    L = list(set(getClobCoords(CDSdict)))
    L.sort()
    
    start = 0
    NonCDScords = []
    for e in L:
        NonCDScords.append((start,e[0]-2))
        start = e[1]
    print len(NonCDScords), len(L)
    for cds in NonCDScords:
        r = ref[cds[0]:cds[1]+1]
        if r!="":
            NonCDS['r:'+str(cds[0])+'-'+str(cds[1]+1)] = r
    return NonCDS

def clean_read(s):
    res = ""
    R = s.split('\n')
    for r in R:
        res = res+r
    return res;

def getrate(clf, read, l= 150):
    l = min(len(read),l)
    frames = [read[k:k+l] for k in range(0,len(read),l)]
    #print frames
    #print frames
    objs = [getfeautures(k) for k in frames]
    #feas = getfeautures(D[R7[4]])
    res = clf.predict(objs)
    #print res
    return sum(res)*1.0/len(res)

class Graph:
    D_str = {}
    D_edge = {}
    D_rate = {}
    def __init__(self, path):
        f = open(path,'r')
        L = []
        D1 = {}
        D_str = {}

        last_node = ""
        for s in f.readlines():
            if s[0]=='>':
                last_node = s[:-2].split(':')[0][1:]
                L.append(s)
                D1[last_node] = 0
                D_str[last_node] = ""
            else:
                D_str[last_node] = D_str[last_node] + s;
                D1[last_node] = D1[last_node]+len(s)-1

        for a in D_str:
            D_str[a] = clean_read(D_str[a])
        D = {}
        vertices = []
        edges = []

        D_edges = {}

        for i,s in enumerate(L):
            s_node = s[:-2].split(':')[0][1:]
            D[s_node] = i
            vertices.append(s_node)
            K = s[:-2].split(':')
            D_edges[s_node] = []
            if (len(K)>1):
                E = s[:-2].split(':')[1].split(',')
                for e in E:
                    D_edges[s_node].append(e)
        
        self.D_str = D_str
        self.D_edges = D_edges
        
    def getmaxedge(self, E,  D, p = True,):
        mx = -1
        es = ""
        if p:
            for e in E:
                if D[e]>mx:
                    mx = D[e]
                    es = e
        else:
            es = self.getmaxedge(E,D)
            mx = D[es]
            for e in E:
                if D[e]<mx:
                    mx = D[e]
                    es = e
        return es

    def extracttracks(self, name, length, D, p = True, S_t = [],   c_len_t = 0):
        S = copy.deepcopy(S_t)
        #print S
        D_str = self.D_str
        c_len = c_len_t
        if S == []:
            S = [[name,D_str[name][:length].upper()]]
        else:
            for i,k in enumerate(S):
                S[i][0] = S[i][0]+'&'+name
                S[i][1] = S[i][1]+D_str[name][:(length-c_len)].upper();

        c_len = len(S[0][1])
        if (c_len == length):
            return S

        Rs = []
        if len(self.D_edges[name])>0:
            ans = self.getmaxedge(self.D_edges[name], D, p)
            Rs2 = self.extracttracks(ans, length, D, p, S,  c_len)
            Rs = Rs+Rs2
        if (Rs == []):
            return [['-','-']]
        return Rs
    
    def countrate(self, clf):
        for d in self.D_str.keys():
            self.D_rate[d]=getrate(clf, self.D_str[d])
    
    def extractreads(self,p, l):
        res = []
        for t in self.D_str.keys():
                res = res+self.extracttracks(t,l,self.D_rate,p)
        res2 = {}
        print len(res)
        count = 0
        for t in res:
            if t[1]!='-':
                #res2['extracted'+str(count)] = t[1]
                res2[t[0]] = t[1]
                count = count+1
        return res2
    def extractnames(self,p, l):
        res = []
        for t in self.D_str.keys():
                res = res+self.extracttracks(t,l,self.D_rate,p)
        res2 = {}
        print len(res)
        count = 0
        for t in res:
            if t[1]!='-':
                res2['extracted'+str(count)] = t[0]
                count = count+1
        return res2



def reverse(r):
    s = ""
    Damins = {'A':'T',
              'T':'A',
              'C':'G',
              'G':'C'}
    for i in range(len(r)-1,-1,-1):
        print i
        s = s+Damins[r[i]]
    return s

def extrcounter(gtf_df, ref):
    D_anno2 = {}
    for g in gtf_df.genes.values:
        D_anno2[g]=gtf_df[gtf_df.genes == g][['name','start','end', 'strand']].values
    D_anno3 = {}
    for k in D_anno2:
        s = ""
        for c in D_anno2[k]:
            s = s+ref[c[1]:c[2]]
        if D_anno2[k][0][3]==-1:
            s = reverse(s)
        D_anno3[k] = s
    return D_anno3

def extrcounter_dataset(gtf_df, fasta_dict):
    D_anno2 = {}
    for g in gtf_df.genes.values:
        D_anno2[g]=gtf_df[gtf_df.genes == g][['name','start','end', 'strand']].values
    D_anno3 = {}
    for k in D_anno2:
        s = ""
        name = D_anno2[k][0][0]
        ref = fasta_dict[name]
        for c in D_anno2[k]:
            s = s+ref[c[1]:c[2]]
        if D_anno2[k][0][3]==-1:
            s = reverse(s)
        D_anno3[k] = s
    return D_anno3

def savefasta(Dfa, name, pref = ""):
    f = open(name,'w')
    fnames = 0
    if pref != "":
        f_names = open(name+'.names','w')
    counter = 0
    for k in Dfa:
        if pref == "":
            f.write('>'+k+'\n')
        else:
            f_names.write(pref+str(counter)+'\t'+k+'\n')
            f.write('>'+pref+str(counter)+'\n')
        f.write(Dfa[k]+'\n')
        counter = counter+1
    f.close()
    if pref != "":
        f_names.close()

def loadgtf(s):
    df = pd.read_csv(s, sep = '\t',header=None)
    genes = [df[8][i].split('"')[1] for i,k in enumerate(df[8])]
    df.columns = ['name','source','type','start','end','p1','strand','p2','genes']
    df = df.drop(df[df.type == 'exon'].index)
    return df

def loadfasta(file):
    fin = open(file,'r')
    currkey = ""
    Dres = {}
    L = fin.readlines()
    print len(L)
    for l in L:
        if l[0]=='>':
            currkey = l[1:-1]
            Dres[currkey] = ""
        else:
            Dres[currkey]=Dres[currkey]+l[:-1]
    return Dres

def loadreffasta(file):
    fin = open(file,'r')
    Ref = ""
    L = fin.readlines()
    print len(L)
    wasenterpoint = False
    for l in L:
        if l[0]=='>':
            if wasenterpoint:
                print "No single sequence"
                return
            else:
                wasenterpoint = True
        else:
            Ref=Ref+l[:-1]
    return Ref

def savereads(resf,out):
        f = open(out,'w')
        D_corr = {}
        D_outf = {}
        counter = 0
        for t in resf[:]:
            f.write('>r'+str(counter)+'\n')
            D_corr['r'+str(counter)] = t[0][0]
            D_outf['r'+str(counter)] = t[0][1]
            counter = counter+1
            f.write(t[0][1]+'\n')
        f.close()
        return D_outf

def loadnames(filenames):
    f = open(filenames,'r')
    Dres = {}
    for k in f.readlines():
        tmp = k.split('\t')
        Dres[tmp[0]]=tmp[1][:-1]
    return Dres
    
def savefirstnames(Dfa, name, pref = ""):
    dR = {}
    counter = 0
    for k in Dfa:
        dR[pref+str(counter)] = k
        counter = counter+1
    return dR

import sys

GRAPH_IN = sys.argv[1]
REF_IN = sys.argv[2]
CDS_IN = sys.argv[3]
RATION_IN = float(sys.argv[4])
TRACKS_OUT = sys.argv[5]
TRACK_LEN = int(sys.argv[6])


DM2L_graph = Graph(GRAPH_IN)
DM2L_chrom = loadreffasta(REF_IN)
DM2L_refCDS = getCDS(CDS_IN)
DM2L_refNonCDS = getNonCDS(DM2L_refCDS, DM2L_chrom)
X = [getfeautures(r) for r in DM2L_refCDS.values()]
X_neg = [getfeautures(r) for r in DM2L_refNonCDS.values()]
k = int(len(X)*RATION_IN)
X_train = X[:k]+X_neg[:k]
X_test = X[k:2*k]+X_neg[k:2*k]
y_train = [1]*len(X[:k])+[0]*len(X_neg[:k])
y_test = [1]*len(X[k:2*k])+[0]*len(X_neg[k:2*k])
from sklearn.model_selection import cross_val_score
from sklearn.metrics import precision_recall_fscore_support as ttest
clf = ensemble.RandomForestClassifier(n_estimators = 25, max_depth = 4, random_state = 1)
print "Cross-validation:", cross_val_score(clf, X_train,y_train, cv = 20).mean()
clf.fit(X_train,y_train)


DM2L_graph.countrate(clf)
DM2L_graph_tracks = DM2L_graph.extractreads(True, TRACK_LEN)
savefasta(DM2L_graph_tracks, TRACKS_OUT, "graph_tracks")

