from Bio import SeqIO
import numpy as np
'''
k_mer(read, k) returns a dictionary of k_mers of the 
input read with key being the index of k_mer and value being the k_mer
'''
def k_mer(read, k):
    size = len(read)-k+1
    output = {}
    for i in range(size-1):
        k_mer = read[i:k+i]
        output[i] = k_mer
    return output
        
'''
parser function to extract all reads from fasta file
and store into a list
'''
def parser(path):
    sequences = [] 
    for record in SeqIO.parse(path, "fasta"):
        read=str(record.seq)
        sequences.append(read)
    return sequences

'''
Generate kmers of size k for all reads, in order
function returns a dicti where key is readId and value is
another dict with key being index, value being kmers
'''
def kmers_for_all_reads(reads, k):
    kmers = {}
    for index, read in enumerate(reads):
        # readID starts at 1
        kmers[index+1]=k_mer(read,k)
    return kmers

'''
Generate kmers dictionary where key is kmer, value is
list of readIDs
D1 = {kmer:readID}
D2 = {kmer:{readID:[kmer_position]}}
'''
def kmerDict(D):
    D1 = {}
    D2 = {}
    for read in D.keys():
        for index, kmer in enumerate(D[read].values()):
            # if kmer in Dict, append its corresponding readID and kmer index
            if kmer in D1:
                # append readID if not been added
                D1[kmer].append(read)
                if read not in D2[kmer]:
                    D2[kmer] = {read: [index]}
                # if readID added, but not kmer index, add kmer index
                else:
                    D2[kmer][read].append(index)
            else:
                D1[kmer] = [read]
                D2[kmer] = {read:[index]}
    return D1, D2
                
# Generate a common kmers frequency dictionary where key is read pairs and value is number of common kmers
def kmerFreqPerPair(D,s):
    FrequencyTable = pair(1, 1056)
    index = 0
    for kmer in D.keys():
        L = D[kmer]
        pairs = pair(D[kmer][0], D[kmer][-1]+1)
        # remove (a,b) where a,b not in pairs since pair generate all possible pairs
        for key in pairs.keys():
            if key not in L:
                del pairs[key]
        for p in pairs.keys():
            print(index)
            FrequencyTable[p]+=1
            index+=1
    # remove readIDs whose common kmers is < some threshold
    return {key : val for key, val in FrequencyTable.items() if val > s}
                          
# function to generate all possible pairs between two numbers n1, n2
def pair(n1, n2):
    # O(n)
    pairDict = {}
    p = np.mgrid[n1:n2,n1:n2]
    p = np.rollaxis(p,0,3)        
    p = p.reshape(((n2-1)*(n2-1), 2))
    p = p[p[:,0] != p[:,1]]
    for index in range(p.shape[0]):
        forward_key = str(p[index][0])+" "+str(p[index][1])
        reversed_key = str(p[index][1])+" "+str(p[index][0])
        if forward_key not in pairDict and reversed_key not in pairDict:
            pairDict[forward_key] = 0  
    return pairDict


path = "C:/Users/Connie/Desktop/COMP561/comp561-project/readsMappingToChr1.fa.txt"
R = parser(path)
D = kmers_for_all_reads(R, 10)
D1,D2 = kmerDict(D)
F = kmerFreqPerPair(D1, 10)

