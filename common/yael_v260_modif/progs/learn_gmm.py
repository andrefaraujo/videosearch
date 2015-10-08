#! /usr/bin/python

import sys, array, random, os

import numpy
from yael import ynumpy



def desc_size(f, fmt):
  """ return matrix size of an open file """
  
  fsize = os.fstat(f.fileno()).st_size
  if fsize == 0: return -1, 0  
  
  if fmt == 'fvecs':      pass # readable at beginnning of file
  elif fmt == 'siftgeo':  f.seek(4 * 9) # size of header
  else:                   assert False, "unknown format %s" % fmt
  
  a = array.array('i')
  a.fromfile(f, 1)
  d = a[0]

  if fmt == 'fvecs':     unitsize = 4 * (d + 1)
  elif fmt == 'siftgeo': unitsize = 4 * 9 + 4 + d
  
  assert fsize % unitsize == 0, "weird size (size %d, dim %d)" % (fsize, d)

  return d, fsize / unitsize


def descs_size(infiles, fmt):
  """ total size of matrices in a set of files """
  n, d = 0, -1  
  for fname in infiles:
    ni, di = desc_size(open(fname, 'r'), fmt)
    if d == -1:    d = di
    else:          assert d == di, "file %s has dim %d, expected %d" % (d2, d)
    n += n2    
  return d, n

def read_descs(infiles, fmt):
  """ read and concatenate matrices from a list of files"""
  vl = []

  for fname in infiles:
    print "reading", fname, "\r",
    sys.stdout.flush()
    
    if fmt == 'fvecs':
      v = ynumpy.fvecs_read(fname)
    elif fmt == 'siftgeo':
      v, meta = ynumpy.siftgeo_read(fname)
      v = v.astype(numpy.float32)
    else: assert False, "unknown format %s" % informat

    if v.shape[1] != 0: vl.append(v.T)

  return numpy.vstack(vl).T
    
def read_descs_subset(infiles, fmt, d, subset):
  """ read and concatenate a subset of columns from a set of files. Subset shoud be sorted."""
  n = len(subset)
  v = numpy.zeros((d, n), dtype = numpy.float32, order = 'FORTRAN')

  # n_begin:n_end: indices of the elements in current file
  n_end, fno = 0, 0 
  for i, j in enumerate(subset): 

    while j >= n_end:
      # need to load new file
      print "reading", infiles[fno], "\r",
      sys.stdout.flush()
      f = open(infiles[fno], 'r')
      d2, n2 = desc_size(f, fmt)
      if n2 >= 0: assert d == d2
      fno += 1
      n_begin, n_end = n_end, n_end + n2

    # ok, we are in the correct file -> seek & read
    if fmt == 'siftgeo':
      f.seek((4 * d + 4) * (j - n_begin) + 4)
      col = numpy.frombuffer(f.read(4 * d), dtype = numpy.float32)      
    elif fmt == 'siftgeo':
      f.seek((4 * 9 + 4 + d) * (j - n_begin) + 4 * 9 + 4)
      col = numpy.frombuffer(f.read(4 * d), dtype = numpy.uint16).astype(numpy.float32)
    v[:, i] = col
  
  return v

def usage():
  print >>sys.stderr, """
learn_gmm.py [options] infile1.fvecs .... infileN.fvecs

Learns a GMM (and optionaly a PCA matrix) from a set of files.

options:

-h --help        this help
-siftgeo         input files are in siftgeo instead of fvecs
-nlearn n        sample n pts from infiles (default 100 * k)
-k k             nb of Gaussians (default 16)
-o out.gmm       output GMM
""" 
  sys.exit(1)

if __name__ == '__main__': 

  infiles = []
  informat = 'fvecs'
  nlearn = -1
  k = 16
  outfile = None
  outfile_pca = None
  pca = -1

  args = sys.argv[:1]
  
  while args:
    a = args.pop(0)

    if a in ('-h', '--help'):   usage()
    elif not a.startswith('-'): infiles.append(a)
    elif a == '-siftgeo':       informat = 'siftgeo'
    elif a == '-k':             k = int(args.pop(0))
    elif a == '-o':             outfile = args.pop(0)    
    else:
      print >> sys.stderr, "unknown arg",a
      sys.exit(1)

  if nlearn < 0: nlearn = 100 * k

  print "scanning %d files (format %s)" % (len(infiles), informat)

  d, n = descs_size(infiles, informat) 
  
  print "found %d pts in %dD" % (n, d)

  if n <= nlearn:
    print "using all points as learning data"

    v = read_descs(infiles, informat)

  else:
    print "sampling %d points"
    subset = random.sample(xrange(n), nlearn)
    subset.sort()
    print "reading files"

    if nlearn * 5 < n:
      print "load all files..."
      vi = read_descs(infiles, informat)
      print "subsample"
      v = numpy.zeros((d, nlearn), dtype = numpy.float32, order = 'FORTRAN')
      for i, j in enumerate(subsample):
        v[:,i] = vi[:,j]
    else:
      print "load subset..."
      vi = read_descs_subset(infiles, informat, subset)

  print "Learning from %d points in %d D" % (n, d)
  
    
    
    
