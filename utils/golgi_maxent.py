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

def matrix_form(n, v):
    a = np.zeros((n,n), dtype=np.bool)
    for (r,c) in [(r,c) for r in range(n) for c in range(n) if r!=c]:
	a[r,c] = v[index_1d(n,r,c)]
    return a

def vector_form(a):
    n = a.shape[0]
    v = np.zeros(n*(n-1)/2, dtype=np.bool)
    for t in range(v.shape[0]):
	v[t] = a[index_2d(n, t)]
    return v

def vector_from_string(s):
    return np.array([int(x) for x in s], dtype=np.bool)

def vector_from_int(n, k):
    return vector_from_string(np.binary_repr(k, width=n*(n-1)/2))

def gamma(n, v):
    '''Variance of the degree distribution of the network.'''
    return matrix_form(n, v).sum(axis=0).var()

def pdf(n, v, lamb, mu, nu):
    # mu*v + nu*(1-v) = f_t(v_t) (remember that v is binary)
    # mu and nu are n(n-1)/2 dimensional parameter vectors. lambda is scalar.
    return (mu*v + nu*(1-v)).prod() * np.exp(lamb*gamma(n, v))

def tth_marginal(n, t, lamb, mu, nu):
    full_space = np.arange(2**(n*(n-1)/2))
    print [np.binary_repr(v, width=n*(n-1)/2) for v in full_space]
    probabilities = np.array([pdf(n,vector_from_int(n, k),lamb,mu,nu) for k in full_space])
    print probabilities
    print probabilities.sum()
    #value_in_one = (np.bitwise_and(space, 2**t)/2**t).dot(probabilities)
    #return value_in_one
    def marginal(n, v_t):
	pass

def marginal_contraint_exact(required_marginal_v):
    pass

def constraint_expression(x, n, required_marginal_v, required_gamma):
    lamb, mu, nu = x
    pass

if __name__ == "__main__":
    n = 5
    marginal_v = np.array([0.59 , 0.38, 0.33, 0.39, 0.11, 0.70, 0.78, 0.36, 0.56, 0.74])
    # our initial guess for the solution is just a product of marginals
    initial_lamb = 0
    initial_mu = marginal_v
    initial_nu = 1 - initial_mu
    # v = np.random.randint(low=0, high=2, size=n*(n-1)/2)
    v1 = np.array([1, 1, 0, 0, 0, 1, 1, 0, 0, 0])
    v2 = np.array([1, 0, 0, 0, 0, 1, 1, 0, 1, 1])
    print matrix_form(n, marginal_v)
    for v in [v1, v2]:
	print v
	print gamma(n, v)
	print pdf(n, v, initial_lamb, initial_mu, initial_nu)


def test_index(n=5):
    # matrix -> vector -> matrix
    for r in range(n):
	for c in range(n):
	    if c == r:
		pass
	    else:
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

def test_transformation(n=5):
    # matrix -> vector -> matrix
    a = np.random.random(size=(n,n), dtype=np.bool)
    a = np.abs((a - np.transpose(a))/2) # symmetrise
    print a
    print vector_form(a)
    print matrix_form(n, vector_form(a))
    assert (a == matrix_form(n, vector_form(a))).all()

    # vector -> matrix -> vector
    v = np.random.random(size=n*(n-1)/2, dtype=np.bool)
    print v
    print matrix_form(n, v)
    print vector_form(matrix_form(n, v))
    assert (v == vector_form(matrix_form(n, v))).all()
