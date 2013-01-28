import numpy as np

def index_1d(n, row, col):
    '''Return linear index for an upper triangular matrix element.'''
    if col <= row:
	raise ValueError, '{0}, {1}'.format(row ,col)
    return n*row + col - (row+1)*(row+2)/2

def index_2d(n,t):
    '''Return upper triangular matrix index pair representation.'''
    # fill matrix from the bottom
    tot = n*(n-1)/2
    entries = 0
    row_bottom = 0
    while entries <= (tot - (t+1)):
	row_bottom += 1
	entries += row_bottom
	#print entries, row_bottom
    row = (n-1)  - row_bottom
    col = (row + 1) + t - (tot - entries)
    return row,col

def test(n=5):
    # matrix -> vector -> matrix
    for r in range(n):
	for c in range(r+1, n):
	    t = index_1d(n,r,c)
	    print '{}, -> {} -> {}'.format((r, c), t, index_2d(n, t))
	    assert (r,c) == index_2d(n, t)
    # vector -> matrix -> vector
    for t in range(n*(n-1)/2):
	r, c = index_2d(n, t)
	print '{} -> {} -> {}'.format(t, (r, c), index_1d(n, r, c))
	assert t == index_1d(n, r, c)
