import csv
from pprint import pprint

def import_tagged(filename):
	sentences = []
	with open(filename) as f:
		reader = csv.reader(f, delimiter=' ')
		for row in reader:
			current_row = row[:-1]
			tags = []
			for pair in current_row:
				tags.append(pair.split("_",1)[-1])
			sentences.append(tags)

	return sentences

def accuracy(pred, real):
	pred = [item for sublist in pred for item in sublist]
	real = [item for sublist in real for item in sublist]

	correct = 0
	for i in range(len(pred)):
		if pred[i] == real[i]:
			correct = correct + 1
	print "accuracy = " + str(100*float(correct)/len(real))

def normalize_1d(matrix):
	tot = sum(matrix.values())
	#maps a function that divides everything by the sum
	normalized = matrix
	if tot != 0:
		normalized = dict(zip(matrix, map(lambda val: val/float(tot), matrix.values())))
	return normalized

def normalize_matrix(matrix):
    #normalizes the probabilities
    #first check for 1d
    if isinstance(matrix.itervalues().next(), dict):
    #map the normalize_1d function on the each value of matrix
        return dict(zip(matrix, map(normalize_1d, matrix.values())))	
    else:
        return normalize_1d(matrix)

def confusion_matrix(pred, real):
	pred = [item for sublist in pred for item in sublist]
	real = [item for sublist in real for item in sublist]
	alltags = pred + real
	alltags = list(set(alltags))
	matrix = dict(zip(alltags, [dict(zip(alltags, [0 for end_state in range(len(alltags))])) for start_state in range(len(alltags))]))
	pairs = zip(pred, real)
	for pair in pairs:
		(first, second) = pair
		matrix[first][second] = matrix[first][second] + 1
	matrix = normalize_matrix(matrix)
	print matrix.keys()
	for key in matrix.keys():
		s = str(key) + " "
		for k2 in matrix[key].keys():
			num = '%.2f' % (matrix[key][k2]*100)
			s = s + num + " "
		print s


if __name__ == "__main__":
	pred = import_tagged("tagged.txt")
	real = import_tagged("files/gold.txt")
	accuracy(pred, real)
	"""
	alltags = [item for sublist in real for item in sublist]
	print "confusion matrix"
	confusion_matrix(pred, real)
	"""

