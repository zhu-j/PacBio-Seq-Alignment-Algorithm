with open('true_pairs.txt', 'rb') as handle:
    data = handle.read()
truePairs = pickle.loads(data)
a = [key[0] for key in truePairs.keys()]
b = [key[1] for key in truePairs.keys()]
truePairs_list = np.unique(a+b)
total_truePairs = len(truePairs.keys())

miss = 0
correct = 0
incorrect = 0
total_readPairs = len(readPairs.keys())
total_truePairs = len(truePairs.keys()) # truePairs = {(R1,R2): (R1 seq, R2 seq)}

correct = 0
incorrect = 0
miss = 0
for pair in readPairs.keys():
    if pair in truePairs.keys() or (pair[1], pair[0]) in truePairs.keys():
        correct += 1
    else:
        if pair[0] in truePairs_list or pair[1] in truePairs_list:
            incorrect +=1
        else:
            miss +=1
    correct_msg = 'correct: {}'.format(correct)
    incorrect_msg = 'incorrect: {}'.format(incorrect)
    miss_msg = 'miss: {}'.format(miss)
    print("total readPairs: {}, total truePairs: {}".format(total_readPairs, total_truePairs))
    print('correct: {}'.format(correct))
    print('incorrect: {}'.format(incorrect))
    print('miss: {}'.format(miss))
    print('\n')
