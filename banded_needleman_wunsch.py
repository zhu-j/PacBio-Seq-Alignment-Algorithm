# -*- coding: utf-8 -*-
"""Banded Needleman Wunsch.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1tPM4j6Hp0ApWZtA_Wr8prwq0M3Qa9E4i
"""

import numpy as np

#this function does the banded Needleman-Wunsch alignment
#X
def banded_needleman_wunsch(read1,read2,X,gap_penalty,match_score,mismatch_score):

  length_read1=len(read1)
  length_read2=len(read2)

  #initialization
  matrix=np.zeros((length_read1+1,length_read2+1))
  #Fill the score for the first slot
  #matrix[1,0]=gap_penalty
  matrix[0,0]=0
  #initilaize the first row and first column within X
  print("alignment start")
  for i in range(1,length_read1+1):
    matrix[i,0]=gap_penalty*(i)
    #Allow gaps at the beginning of read1
  for j in range(1,length_read2+1):
    matrix[0,j]=0
  
  score=np.zeros(3)
  for i in range(1,length_read1+1):
    for d in range(-X,X):
       j=i+d     
       if 1 <= j and j<= length_read2:
         if read1[i-1]==read2[j-1]:
           matrix[i,j]=matrix[i-1,j-1]+match_score
         else:
           matrix[i,j]=matrix[i-1,j-1]+mismatch_score
         if abs((i-1)-j) <= X:
           matrix[i,j]=max(matrix[i,j],matrix[i-1,j]+gap_penalty)
         if abs(i-(j-1)) <= X:
           matrix[i,j]=max(matrix[i,j],matrix[i,j-1]+gap_penalty)
  
  last_columns=list(matrix[:,length_read2])
  optimal_score=max(last_columns)
  starting=last_columns.index(optimal_score)
  match=0
  a=starting
  b=len(read2)
  print(matrix)
  score=matrix[starting][b]
  while(a>0 or b > 0):
    if a>0 and b>0 and ((matrix[a][b] == matrix[a-1][b-1]+match_score) or (matrix[a][b]==matrix[a-1][b-1]+mismatch_score)):
      if matrix[a][b]==matrix[a-1][b-1]+match_score and read1[a-1]==read2[b-1]:
        match +=1
      a-=1
      b-=1
    
    elif b>0 and (matrix[a][b]==matrix[a][b-1]+gap_penalty):

      b-=1
    elif a>0 and (matrix[a][b]==matrix[a-1][b]+gap_penalty):
      a-=1
  print("traceback ended")
  return match

"""read1=['A','C','C']
read2=['A','A']
a=banded_needleman_wunsch(read1,read2,2,-2,1,-1)
b=banded_needleman_wunsch(read2,read1,2,-2,1,-1)"""
