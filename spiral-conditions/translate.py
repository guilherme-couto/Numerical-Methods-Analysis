new_file = open('spiralV.txt', 'w')
with open('./lastSpiralV.txt', 'r') as f:
	for line in f:
		line = line.split()
		for value in line:
			new_file.write(value)
			new_file.write('\n')

new_file = open('spiralW.txt', 'w')
with open('./lastSpiralW.txt', 'r') as f:
	for line in f:
		line = line.split()
		for value in line:
			new_file.write(value)
			new_file.write('\n')
