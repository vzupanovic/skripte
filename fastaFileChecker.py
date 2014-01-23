def checkFile(fastaFile):
	validAlpha = ['A', 'T', 'G', 'C']
	stream = open(fastaFile, 'r')
	data = stream.readlines()
	headers = []
	for line in data:
		line = line.strip();
		if line[0] == ">":
			headers.append(line)
		else:
			for letter in line:
				if letter not in validAlpha:
					return -1
			
	if len(headers) < 1:
		return -1
		
	return 1

	
