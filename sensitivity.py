from banded_needleman_wunsch import *
import pickle
import ast
import time

file = open("read_pairs.txt", 'rb')
readPairs = pickle.load(file)

with open('true_pairs.txt', 'rb') as handle:
    data = handle.read()
truePairs = pickle.loads(data)
a = [key[0] for key in truePairs.keys()]
b = [key[1] for key in truePairs.keys()]
truePairs_list = np.unique(a+b)
read0 = [key[0] for key in readPairs.keys()]
read1 = [key[1] for key in readPairs.keys()]
readPairs_list = np.unique(read0+read1)
total_truePairs = len(truePairs.keys())
wrong_list = []
wrong_list1 = []
incorrect_list = []
correct_list = []
miss = 0
correct = 0
incorrect = 0
total_readPairs = len(readPairs.keys())
total_truePairs = len(truePairs.keys())
cor_list = []
correct = 0
wrong_match = 0
incorrect = 0
miss = 0

for pair in readPairs.keys():
    reverse = (pair[1], pair[0])
    if pair in truePairs.keys() or reverse in truePairs.keys():
        correct += 1
        if pair in truePairs.keys():
            correct_list.append(pair)
        elif reverse in truePairs.keys():
            correct_list.append(reverse)
    else:
        if pair[0] in truePairs_list or pair[1] in truePairs_list:
            wrong_match +=1
            wrong_list.append(pair[0])
            wrong_list1.append(pair[1])
        else:
            incorrect += 1
            incorrect_list.append(pair)
for pair in truePairs.keys():
    if pair not in correct_list:
            miss +=1
print(correct)
print(miss)
print(correct_list)

with open('sensitivity_full.txt', 'w') as sens:
    sens.write("total truePairs: {}".format(total_truePairs) + '\n')
    sens.write("total readPairs: {}".format(total_readPairs) + '\n')
    sens.write('correct: {}'.format(correct) + '\n')
    sens.write('wrong match: {}'.format(wrong_match) + '\n')
    sens.write('incorrect: {}'.format(incorrect) + '\n')
    sens.write('miss: {}'.format(miss) + '\n')

