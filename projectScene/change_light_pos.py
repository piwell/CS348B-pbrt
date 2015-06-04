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

filename = 'volumescene_png.pbrt'

dim = 1 #in [x, y, z]
start = -5
end = 5
steps = 10
step = float(end-start)/steps

filedata = None
with open (filename, "r") as myfile:
	filedata = myfile.read()


# def setNumSamplesPerPixel(pbrtData, samplesPerPixel):
# 	samplerLine = re.search('Sampler.+', pbrtData).group()
# 	newSamplerLine = re.sub('[\d+]', str(samplesPerPixel), samplerLine)
# 	return re.sub('Sampler.+', newSamplerLine, pbrtData)

# def setNumSamplesPerLight(pbrtData, samplesPerLight):
# 	samplerLine = re.search('AreaLightSource.+', pbrtData).group()
# 	newSamplerLine = re.sub('\[(\d+)]', '['+str(samplesPerLight)+']', samplerLine)
# 	return re.sub('AreaLightSource.+', newSamplerLine, pbrtData)

def setLightPos(pbrtData, dim, value):
	coordinatesPattern = '\nLightSource.+"point from"\s+\[([0-9\-]+)\s+([0-9\-]+)\s+([0-9\-]+)\s?]'
	res = re.search(coordinatesPattern, pbrtData)
	pointFrom = res.group(0)
	newPointFrom = pointFrom[:pointFrom.rfind('[')+1]
	for i in xrange(0,3):
		if i == dim:
			newPointFrom += str(value) + " "
		else:
			newPointFrom += str(res.group(i+1)) + " "
	newPointFrom+= ']'
	
	#print re.search(pointFrom, pbrtData)
	return re.sub('\nLightSource.+"point from"\s+\[[0-9\-\s]+\]', newPointFrom, pbrtData)

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
