import Utils as Utils

def Norm(Q, u):
	return sqrt(Norm2(Q, u))

def Norm2(Q, u):
	return InnerProduct(Q, u, u)

def InnerProduct(Q, u, v):
	return u * Q * v

def Proj(Q, u, v):
	return (InnerProduct(Q, u, v) / InnerProduct(Q, u, u)) * u

def GramSchmidt(Q, B):
	N, N = Q.dimensions()
	B_ = [zero_vector(N) for _ in range(N)]

	for i in range(N):
		B_[i] = B[i]
		for j in range(i):
			B_[i] = B_[i] - Proj(Q, B_[j], B_[i])

	return B_

def SampleD(Q, s):
	N, N = Q.dimensions()
	B = identity_matrix(N)
	B = SortByColumnNorm(Q, B).columns()

	vi = zero_vector(ZZ, N)
	ci = zero_vector(ZZ, N)

	b_ = GramSchmidt(Q, B)

	for i in range(N-1, -1, -1):
		bi_ = b_[i]
		c_ = InnerProduct(Q, ci, bi_) / Norm2(Q, bi_)
		s_ = s / Norm(Q, bi_)

		zi = Utils.SampleZ(s_, c_)
		bi = B[i]
		ci = ci - zi * bi
		vi = vi + zi * bi

	return vi

def SortByColumnNorm(Q, B):
	N, N = B.dimensions()
	columns = B.columns()
	columns.sort(key=lambda v : Norm2(Q, v))
	return block_matrix([[matrix(column).transpose() for column in columns]])
