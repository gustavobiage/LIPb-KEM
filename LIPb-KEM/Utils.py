

# This file was *autogenerated* from the file Utils.sage
from sage.all_cmdline import *   # import sage library

_sage_const_3 = Integer(3); _sage_const_4 = Integer(4); _sage_const_0 = Integer(0); _sage_const_1 = Integer(1); _sage_const_2 = Integer(2); _sage_const_128 = Integer(128)
def LLL(B, delta=_sage_const_3 /_sage_const_4 ):
	return B.transpose().LLL().transpose()

def GramSchmidt(B):
	return B.transpose().gram_schmidt(orthonormal=False)[_sage_const_0 ].transpose()

def GaussianHeuristic(B):
	N, M = B.dimensions()
	return abs(B.det())**(_sage_const_1 /N) * sqrt(N / (_sage_const_2  * pi * e))

def IntegralLatticeGap(B, shortest_vector_length):
	return GaussianHeuristic(B) / shortest_vector_length

def InnerProduct(v1, v2):
	return (v1.transpose() * v2)[_sage_const_0 ]

def Norm2(v):
	return v.dot_product(v)

def BasisNorm(B):
	return max(v.norm() for v in B.columns())

def Coeff(B, v):
	N, N = B.dimensions()
	return vector(ZZ, N, [ZZ(x) for x in B * BackslashOperator() * v])

def RoundZZ(r):
	return ZZ(round(r))

def SwapColumns(B, i, j):
	aux = B.column(i)
	InsertColumn(B, i, B.column(j))
	InsertColumn(B, j, aux)

def InsertColumn(A, c, v):
	N, M = A.dimensions()
	assert _sage_const_0  <= c and c < M
	for i in range(N):
		A[i,c] = v[i]

def SortByColumnNorm(B):
	N, N = B.dimensions()
	columns = B.columns()
	columns.sort(key=lambda v : Norm2(v))
	return block_matrix([[matrix(column).transpose() for column in columns]])

def Rho(x, s, c):
	return exp(-pi * Norm2(x - c)/ s**_sage_const_2 )
N = _sage_const_128 
def T(n):
	return log(n, _sage_const_2 )

def SampleZ(s, c):
	lower_bound = ceil(c - s * T(N))
	upper_bound = floor(c + s * T(N))
	x = ZZ.random_element(lower_bound, upper_bound)
	x_ = vector(ZZ, _sage_const_1 , [x])
	c_ = vector(QQ, _sage_const_1 , [c])
	return x if random() <= Rho(x_, s, c_) else SampleZ(s, c)

def SampleD(B, s, c):
	N, N  = B.dimensions()
	b = B.columns()

	assert c in VectorSpace(RR, N)
	assert all(Norm2(b[i]) <= Norm2(b[i+_sage_const_1 ]) for i in range(N-_sage_const_1 ))
	
	vi = zero_vector(ZZ, N)
	ci = c
	b_ = GramSchmidt(B).columns()

	for i in range(N-_sage_const_1 , -_sage_const_1 , -_sage_const_1 ):
		bi_ = b_[i]
		c_ = ci.dot_product(bi_) / Norm2(bi_)
		s_ = s / bi_.norm()

		zi = SampleZ(s_, c_)
		bi = b[i]
		ci = ci - zi * bi
		vi = vi + zi * bi
	return vi

def isLatticeVector(B, v):
	N, N = B.dimensions()
	sol = B * BackslashOperator() * v # Solve system
	# A vector belongs to the lattice if it is a linear combination of the basis
	# WITH INTEGER COEFFFICIENTS.
	return all(x in ZZ for x in sol)

