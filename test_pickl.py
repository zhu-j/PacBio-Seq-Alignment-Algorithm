import _pickle as pickle
file = open("read_pairs.txt", 'rb')
print(pickle.load(file))