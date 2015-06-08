import sys
import re
import math

def usage():
	print ""
	print "Usage:"
	print "    python change_light_pos <scene> | ./bin/pbrt"
	print ""

# if len(sys.argv) < 2:
# 	usage()
# 	sys.exit()
#filename = sys.argv[1]

filename = 'final_png.pbrt'

dim = 0 #in [x, y, z]
start = 0
end = 3
steps = 10
step = float(end-start)/steps
attr = 'point from'

filedata = None
with open (filename, "r") as myfile:
	filedata = myfile.read()


def setLightPos(pbrtData, dim, value):
	
	coordinatesPattern = '\nLightSource.+"' + attr + '"\s+\[([0-9\-]+)\s+([0-9\-]+)\s+([0-9\-]+)\s?]'
	res = re.search(coordinatesPattern, pbrtData)
	pointAttr = res.group(0)
	newPointAttr = pointAttr[:pointAttr.rfind('[')+1]
	for i in xrange(0,3):
		if i == dim:
			newPointAttr += str(value) + " "
		else:
			newPointAttr += str(res.group(i+1)) + " "
	newPointAttr+= ']'
	
	return re.sub('\nLightSource.+"' + attr + '"\s+\[[0-9\-\s]+\]', newPointAttr, pbrtData)

def appendToOutputFileName(pbrtData, extension):
	fileNameLine = re.search('Film.+"string filename".+', pbrtData).group()
	newFileNameLine = re.sub('.png', extension+".png", fileNameLine)
	return re.sub('Film.+"string filename".+', newFileNameLine, pbrtData)

for i in xrange(steps):
	value = start + i * step
	pbrtInput = filedata

	extension = '_light_' + 'xyz'[dim] + '=%.2f' % value

	pbrtInput = appendToOutputFileName(pbrtInput, extension)
	pbrtInput = setLightPos(pbrtInput, dim, value)
	
	print pbrtInput
