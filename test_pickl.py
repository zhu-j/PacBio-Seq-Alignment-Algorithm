import _pickle as pickle
file = open("readPairs_900.txt", 'rb')
print(pickle.load(file))