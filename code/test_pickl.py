import _pickle as pickle
#file = open("true_pairs.txt", 'rb')
file2 = open("read_pairs.txt", 'rb')
output2 = pickle.load(file2)
#output = pickle.load(file)
#print(output.keys())
print(output2.keys())
