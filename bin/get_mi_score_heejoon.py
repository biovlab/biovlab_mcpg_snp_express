#!/usr/bin/env python
import math
import sys
#from my_argparse import *
#from make_csv2 import *
import random
import numpy as np
#from random_permutation_mi import *


pseudo_cnt = 0.000000001
def mi(vec1, vec2):
	len1 = len(vec1)
	len2 = len(vec2)
	#p1 = {0: 0.5, 1: 0.5}
	p1 = {}
	p2 = {}
	for elem in vec1:
		if not elem in p1:
			p1[elem] = 0.0
		p1[elem] += 1.0
	for key in p1.keys():
		p1[key] = (p1[key] / len1)
	for elem in vec2:
		if not elem in p2:
			p2[elem] = 0.0
		p2[elem] += 1.0
	for key in p2.keys():
		p2[key] = (p2[key] / len2)
	
	p12 = {}
	for k1 in p1.keys():
		for k2 in p2.keys():
			k12 = str(k1)+','+str(k2)
			if not k12 in p12:
				p12[k12] = pseudo_cnt
	for i in range(len1):
		try:
			key = str(vec1[i])+','+str(vec2[i])
			p12[key] += 1.0
		except IndexError:
			print "i : " + str(i)
			print "vec1 : "
			print vec1
			print "vec2 : "
			print vec2
			print "p12 : "
			print p12
	#print len1
	#print p12
	for key in p12.keys():
		p12[key] = p12[key] / len1
	
	#print p12
	mi = 0.0
	for k1 in p1.keys():
		for k2 in p2.keys():
			k12 = str(k1)+','+str(k2)
			mi += p12[k12]*math.log((p12[k12]/(p1[k1]*p2[k2])),2)
	#print mi
	return mi

def rand_permut(clsList, bin_num):
	numSample = len(clsList)
	cls2idx = {}
	idx=0
	for cls in clsList:
		try:
			temp=cls2idx[cls]
		except KeyError:
			cls2idx[cls]=idx
			idx+=1
	
	ranks = []
	#for cls in clsList:
	#	ranks.append(cls2idx[cls])
	for i in range(len(clsList)):
		r = (i / bin_num) + 1
		ranks.append(r)

	res = []
	for i in range(100000):
		res.append(mi(clsList, ranks))
		random.shuffle(ranks)
		#print ranks
		#print res[-1]
	hist, bin_edges = np.histogram(res, bins=50)
	return hist, bin_edges

def get_pvalue(hist, bin_edges, score):
	sum_hist = sum(hist)
	pre_pval = 0
	for i in range(0, len(hist)):
		if bin_edges[i+1] >= score:
			pre_pval += hist[i]
	pvalue = float(pre_pval) / sum_hist
	return pvalue



def mi2(vec1, vec2):
	len1 = len(vec1)
	len2 = len(vec2)
	#p1 = {0: 0.5, 1: 0.5}
	p1 = {}
	p2 = {}
	for elem in vec1:
		try:
			p1[elem] += 1.0
		except KeyError:
			p1[elem] = 0.0
	for key in p1.keys():
		p1[key] = (p1[key] / len1)
	for elem in vec2:
		try:
			p2[elem] += 1.0
		except KeyError:
			p2[elem] = 0.0
	for key in p2.keys():
		p2[key] = (p2[key] / len2)
	
	p12 = {}
	for k1 in p1.keys():
		for k2 in p2.keys():
			k12 = str(k1)+','+str(k2)
			p12[k12] = pseudo_cnt

	for i in range(len1):
		key = str(vec1[i])+','+str(vec2[i])
		p12[key] += 1.0
	#print len1
	#print p12
	for key in p12.keys():
		p12[key] = p12[key] / len1
	
	#print p12
	mi = 0.0
	for k1 in p1.keys():
		for k2 in p2.keys():
			k12 = str(k1)+','+str(k2)
			mi += p12[k12]*math.log((p12[k12]/(p1[k1]*p2[k2])),2)
	#print mi
	return mi

def discretize_values(values, bin_num):
	sorted_vals = sorted(values)
	length = len(sorted_vals)
	bin_edges = []
	for i in range(bin_num):
	#	try:
		bin_edges.append(sorted_vals[length * i / bin_num])
	#	except IndexError:
	#		print "length : " + str(length)
	#		print "i : " + str(i)
	#		print "bin_num : " + str(bin_num)
	#		print "values : " 
	#		print values

	bin_edges.append(float('inf'))
	#print bin_edges
	ranks = []
	for val in values:
		for i in range(bin_num):
			if bin_edges[i] <= val and val < bin_edges[i+1]:
				ranks.append(i)
				break
	#print ranks
	return ranks


if __name__=='__main__':

	# replace with first line of stdin, that is put these value to stdin from previos file
	#clsList = ["a","a","a","b","b","b","c","c"]
	#clsList = sys.stdin.readline()[:-1].split(',')
	#clsNum = int(sys.stdin.readline()[:-1])

	clsList_fd  = open("temp1.txt", "r")
	clsList = clsList_fd.readline()[:-1].split(',')
	clsList_fd.close()
	clsNum_fd  = open("temp2.txt", "r")
	clsNum = int(clsNum_fd.readline()[:-1])
	clsNum_fd.close()
	#print clsNum

	hist, bin_edges = rand_permut(clsList, clsNum)

	for line in sys.stdin:
		tokens = line[:-1].split('\t')
		id = tokens[0]
		values=[]
		for value in tokens[1:]:
			values.append(float(value))
		#print id
		#print values
		try:
			disc_values = discretize_values(values, clsNum)
		except IndexError:
			print "line : "
			print line 
			print "clsList : "
			print clsList
			print "clsNum : "
			print clsNum
			sys.exit(0)
		#print disc_values
		#print str(mi(clsList, disc_values))
		score=mi(clsList, disc_values)

		try:
			pvalue = get_pvalue(hist, bin_edges, float(score))
		except ValueError:
			print "score : " + str(score)
		print str(id) + "\t" + str(score) + "\t" + str(pvalue)

	#output_tuples.sort(key=lambda tuple: tuple[2], reverse=True)


	#print 'Order\tName\t\tScore'
	#i=1
	#for tup in output_tuples[0:info.ntop]:
	#for tup in output_tuples:
	#	print str(i)+'\t'+tup[0]+'\t\t'+str(tup[1])+'\t'+str(tup[2])
		#print tup[0]
	#	i+=1
	
	#make_csv_in_outdir(info.outdir, output_tuples, info.fin_name)
	#if(info.outdir != None):
	#	make_csv(info.outdir, output_tuples, info.fin_name, fmat.clsList, 50)
