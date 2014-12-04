import argparse
import subprocess
import time
import json
import matplotlib.pyplot as plt
from itertools import groupby
filenames = {
			 'insert' : './insert_throughput.out',
			 'read'   : './read_throughput.out',
			 'resizing': './insert_throughput.out',
			}
tempdir = 'tmp/'
plotdir = 'plot/'
power = 22
numtrials = 3
def generateResizingData():
	test_type = 'resizing'
	startpower = power
	endpower = 10
	filename = tempdir + test_type + str(int(time.time()))
	load = 80
	fp = open(filename, 'w')
	for n in xrange(numtrials):
		for i in xrange(startpower, endpower, -1):
			subprocess.call(
				[
					filenames[test_type], '--end-load', str(load),
					'--power', str(i) 
				],
				stdout = fp
				)
			load *= 2
	return filename

def generateInsertData():
	test_type = 'insert'
	startload = 0
	endload = 90 
	increment = 5
	filename = tempdir + test_type + str(int(time.time()))
	fp = open(filename, 'w')
	for n in xrange(numtrials):
		for i in xrange(startload,endload,increment):
			subprocess.call(
				[
					filenames[test_type], '--begin-load', str(i), '--end-load', str(i+increment),
					'--power', str(power) 
				],
				stdout = fp
				)
	return filename

def generateReadData():
	test_type = 'read'
	startload = 5
	endload = 90
	increment = 5
	filename = tempdir + test_type + str(int(time.time()))
	fp = open(filename, 'w')
	for n in xrange(numtrials):
		for i in xrange(startload,endload,increment):
			subprocess.call(
				[
					filenames[test_type], '--load', str(i),
					'--power', str(power) 
				],
				stdout = fp
				)
	return filename

def generateMixedWorkloadData():
	return

def parseJson( filename ):
	content = []
	with open(filename) as f:
	    for line in f:
	        while True:
	            try:
	                jfile = json.loads(line)
	                break
	            except ValueError:
	            	try:
	            		line += next(f)
	            	except StopIteration:
	            		return content
	                # Not yet a complete JSON value
	        content.append( jfile )
	return content

def parseLatency( filename ):
	info = []
	with open(filename) as f:
		for line in f:
			if line[0] == '{':
				break
			s = line.split(",")
			x = {}
			x["ts"] = int(s[0])
			x["latency"] = int(s[1])
			info.append(x)
	return info

def grouping(data, groupfunc, aggregatefunc):
	groups = []
	uniquekeys = []
	sorteddata = sorted( data, key=groupfunc )
	for k, g in groupby(sorteddata, groupfunc):
   		groups.append(list(g))    # Store group iterator as a list
   		uniquekeys.append(k)
   	ans = [ aggregatefunc(x) for x in groups]
   	return uniquekeys, ans

def generateGraph( info, xdata, ydata, labels):
	fig = plt.figure()
	for i in zip(xdata,ydata,labels):
		plt.plot(i[0],i[1], label=i[2])
	plt.xlabel(info['xlabel'])
	plt.ylabel(info['ylabel'])
	plt.ylim(ymin = 0, ymax=1.4*max(max( ydata, key=lambda r: max(r))))
	plt.legend(loc='upper left')
	fig.suptitle( info['title'] )
	fig.savefig( plotdir + info['filename'].replace(tempdir,"") + '.png')

def generateInsertGraph( filename ):
	content = parseJson(filename)
	xvals, yvals = grouping( content, lambda x: (int) (100*x['end_load']), lambda x: sum(r['throughput'] for r in x)/len(x))
	info = { 'xlabel' : 'load(%)', 'ylabel' : 'throughput(ops/s)', 'title' : "Insert throughput vs load", 'filename' : filename}
	generateGraph( info, [xvals], [yvals], ["Cuckoo table"])

def generateReadGraph( filename ):
	content = parseJson(filename)
	xvals, yvals = grouping( content, lambda x: (int) (100*x['load']), lambda x: sum(r['throughput'] for r in x)/len(x))
	info = { 'xlabel' : 'load(%)', 'ylabel' : 'throughput(ops/s)', 'title' : "Read throughput vs load", 'filename' : filename}
	generateGraph( info, [xvals], [yvals], ["Cuckoo table"])

def generateResizingGraph( filename ):
	content = parseJson(filename)
	xvals, yvals = grouping( content, lambda x: x['initial_table_size'], lambda x: sum(r['throughput'] for r in x)/len(x))
	info = { 'xlabel' : 'Initial table size', 'ylabel' : 'throughput(ops/s)', 
			 'title' : "Resizing throughput to insert {} keys".format(content[0]['num_inserts']), 'filename' : filename}
	generateGraph( info, [xvals], [yvals], ["Cuckoo table"])

def generateLatencyGraph( filename ):
	content = parseLatency( filename )
	xvals, yvals = grouping( content, lambda x: x['ts'], lambda x: sum(r['latency'] for r in x)/len(x))
	info = { 'xlabel' : 'Time', 'ylabel' : 'Latency', 
			 'title' : "Latency while inserting with two resizes",
			 'filename' : filename }
	generateGraph( info, [xvals], [yvals], ["Cuckoo table"])

def main():
	parser = argparse.ArgumentParser(description='Run benchmarks and generate graphs')
	parser.add_argument('-g', '--graphs', action='store_true', help='Whether to generate graphs')
	parser.add_argument('-i', '--inserts', action='store_true')
	parser.add_argument('-z', '--resizing', action='store_true')
	parser.add_argument('-r', '--reads', action='store_true')
	args = vars(parser.parse_args())
	
	#generateLatencyGraph('tmp/constructor_latencyv2.txt')
	#generateLatencyGraph('tmp/mmap_latencyv2.txt')
	#generateLatencyGraph('tmp/calloc_latencyv2.txt')
	#generateInsertGraph('tmp/insert1409090923')
	#generateReadGraph('tmp/read1409091097')
	#generateResizingGraph('tmp/resizing1409091759')
	
	if( args['inserts']):
		insertData = generateInsertData()
		if( args['graphs']):
			generateInsertGraph(insertData)

	if( args['reads']):
		readData = generateReadData()
		if( args['graphs']):
			generateReadGraph(readData)	

	if( args['resizing'] ):
		resizingData = generateResizingData()
		if( args['graphs']):
			generateResizingGraph(resizingData)
main()
