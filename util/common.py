#!/usr/bin/env python

import numpy as np

def mkraisincuts(fr):

	iCut = np.where(fr.SNRMAX1 > 10)[0]
	for k in fr.__dict__.keys():
		fr.__dict__[k] = fr.__dict__[k][iCut]
	return fr

def concat_simple(fr1,fr2):
	for k in fr1.__dict__.keys():
		fr1.__dict__[k] = \
			np.append(fr1.__dict__[k],
		  			  fr2.__dict__[k])
	return fr1

def sort_fitres(fr,cols):
	for k in fr.__dict__.keys():
		fr.__dict__[k] = fr.__dict__[k][cols]
	return fr
