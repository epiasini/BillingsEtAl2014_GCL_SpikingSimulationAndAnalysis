import numpy as np

def index_1d(n, row, col):
    '''
    Return linear index for an upper triangular matrix element. If a
    lower triangular element (i,j: j<i) is given, the linear index for
    (j,i) will be returned instead.
    '''
    if col < row:
	col_old = col
	col = row
	row = col_old
    elif col == row:
	raise ValueError, "Entries on the diagonal are not represented."
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
	for c in range(r) + range(r+1, n):
	    t = index_1d(n,r,c)
	    print '{}, -> {} -> {}'.format((r, c), t, index_2d(n, t))
	    if c < r:
		# if on the lower triangle, transforming to vector
		# notation and back should bring us on the upper
		# triangle
		assert (c, r) == index_2d(n, t)
	    else:
		assert (r, c) == index_2d(n, t)
    # vector -> matrix -> vector
    for t in range(n*(n-1)/2):
	r, c = index_2d(n, t)
	print '{} -> {} -> {}'.format(t, (r, c), index_1d(n, r, c))
	assert t == index_1d(n, r, c)
