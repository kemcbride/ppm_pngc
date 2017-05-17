file = open("scaffold_distances", "r")

file.readline()
line = file.readline().strip()
results = []
while line != "":
	if line.startswith("#"):
		name = line
		line = file.readline().strip()
	else:
		current = line.split()
		id = current[0]
		result = current
		while line.startswith(id):
			if float(current[2]) < float(result[2]):
				result = current
			line = file.readline().strip()
			current = line.split()
		results.append([name, result])
for x in results:
	print (x[0])
	print ("{0}\t{1}\t{2}".format(x[1][0], x[1][1], x[1][2]))
file.close()

