import numpy as np
import scipy.stats as st

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

def matrix_from_vector(n, v):
    a = np.zeros((n,n), dtype=np.bool)
    for (r,c) in [(r,c) for r in range(n) for c in range(n) if r!=c]:
	a[r,c] = v[index_1d(n,r,c)]
    return a

def vector_from_matrix(a):
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
    return matrix_from_vector(n, v).sum(axis=0).var()

def pdf(n, v, lamb, mu, nu):
    # mu[t]*(1-v[t]) + nu[t]*v[t] = f_t(v_t) (remember that v is binary)
    # mu and nu are n(n-1)/2 dimensional parameter vectors. lambda is scalar.
    return (mu*(1-v) + nu*v).prod() * np.exp(lamb*gamma(n, v))

def constraint_expression(n, lamb, mu, nu, required_marginals, required_gamma):
    assert mu.shape[0] == nu.shape[0]
    assert mu.shape[0] == n*(n-1)/2
    # each 'adjacency vector' (vector form of an adjacency matrix)
    # corresponds to the binary representation of a number between 0
    # and 2**(n*(n-1)/2)
    full_space = np.arange(2**(n*(n-1)/2))
    probabilities = np.array([pdf(n,vector_from_int(n, m),lamb,mu,nu) for m in full_space]).reshape([2 for each in range(n*(n-1)/2)])
    # for c in zip([np.binary_repr(m) for m in full_space], probabilities.flat):
    # 	print c

    # evaluate marginals
    marginals = np.array([np.squeeze(p) for p in st.contingency.margins(probabilities)])

    # evaluate expectation value of gamma
    e_gamma = (probabilities * np.array([gamma(n, vector_from_int(n, m)) for m in full_space]).reshape([2 for each in range(n*(n-1)/2)])).sum().reshape(1,)

    # return flattened array of functions that we want to find the root of
    return np.concatenate((e_gamma-required_gamma, np.ravel(marginals-required_marginals)))

if __name__ == "__main__":
    n = 5
    marginal_v = np.array([0.59 , 0.38, 0.33, 0.39, 0.11, 0.70, 0.78, 0.36, 0.56, 0.74])
    # our initial guess for the solution is just a product of marginals
    initial_lamb = 0
    initial_nu = marginal_v
    initial_mu = 1 - initial_nu
    # v = np.random.randint(low=0, high=2, size=n*(n-1)/2)
    v1 = np.array([1, 1, 0, 0, 0, 1, 1, 0, 0, 0])
    v2 = np.array([1, 0, 0, 0, 0, 1, 1, 0, 1, 1])
    print matrix_from_vector(n, marginal_v)
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
    a = np.asarray(np.random.randint(low=0, high=2, size=(n,n)), dtype=np.bool)
    for r in range(n):
	for c in range(0,r):
	    a[r,c] = a[c,r]
	a[r,r] = False
    print a
    print vector_from_matrix(a)
    print matrix_from_vector(n, vector_from_matrix(a))
    assert (a == matrix_from_vector(n, vector_from_matrix(a))).all()

    # vector -> matrix -> vector
    v = np.asarray(np.random.randint(low=0, high=2, size=n*(n-1)/2), dtype=np.bool)
    print v
    print matrix_from_vector(n, v)
    print vector_from_matrix(matrix_from_vector(n, v))
    assert (v == vector_from_matrix(matrix_from_vector(n, v))).all()

def test_constraint_expression():
    goodness = constraint_expression(n=3,
				     lamb=0,
				     mu=np.array([0.4, 0.3, 0.8]),
				     nu=np.array([0.6, 0.7, 0.2]),
				     required_marginals=np.array([[0.4,0.6], [0.3,0.7], [0.8,0.2]]),
				     required_gamma=0.18222222222222223)
    print goodness
    assert np.square(goodness).sum() < 1.e-30
