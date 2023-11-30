os.system('sage --preparse Utils.sage')
os.system('mv Utils.sage.py Utils.py')
import Utils as Utils

# Gaussian Integers
GG = ZZ[I]
GQ = QQ[I]
PHI = GG(1 + I)

def KroneckerProduct(A, B):
	m, n = A.dimensions()
	p, q = B.dimensions()
	return block_matrix([[A[i, j] * B for j in range(n)] for i in range(m)], subdivide=False)

def BarnesWall(n):
	BW = matrix(GG, 2, 2, [[1,	 1],
					   	   [0, 1 + i]]);
	BWN = BW
	for _ in range(n-1):
		BWN = KroneckerProduct(BWN, BW)

	N, N = BWN.dimensions()
	return (BWN, sqrt(N), sqrt(N)/2)

def BWvector(B, LOWER_BOUND=-10, UPPER_BOUND=10):
	N, N = B.dimensions()
	coeff = vector(GG, N, [(ZZ(randint(LOWER_BOUND, UPPER_BOUND)) + I*ZZ(randint(LOWER_BOUND, UPPER_BOUND))) for _ in range(N)])
	return B*coeff

def BasisOverZZ(B):
	N, N = B.dimensions()
	return block_matrix([[ matrix(ZZ, 2, 2, [[ZZ(   R(B[i, j])), ZZ(Img(B[i, j]))],
											 [ZZ(-Img(B[i, j])), ZZ(  R(B[i, j]))]]) for j in range(N)] for i in range(N)], subdivide=False);

def VectorOverCC(v):
	N = len(v)
	return vector(CC, [CC(v[i+1] + v[i] * I) for i in range(0, N-1, 2)])

def VectorOverGG(v):
	N = len(v)
	return vector(GG, [GG(v[i+1] + v[i] * I) for i in range(0, N-1, 2)])

def VectorOverZZ(v):
	N = len(v)
	return vector(ZZ, block_matrix([ [ matrix(ZZ, 1, 2, [[ZZ(Img(x)), ZZ(  R(x))]]) for x in v]])[0])

def isBWvector(B, v):
	N, N = B.dimensions()
	BZ = BasisOverZZ(B)
	vz = VectorOverZZ(v)
	return Utils.isLatticeVector(BZ, vz)

def Error(B):
	N, N = B.dimensions()
	
	dmin2 = N/4 # Minimum correction distance squared
	errorbound = RR(sqrt(dmin2)/3)

	# First we sample a simple vector in the form (error, 0, ..., 0).
	simple_error = vector(QQ, 2*N, [0 if i != 0 else random()*errorbound for i in range(2*N)])
	return vector(CC, N, [(simple_error[i+1] + I * simple_error[i]) for i in range(0, 2*N, 2)])

	# Then, we sample a random rotation using the QR decomposition of a random matrix.
	G = MatrixSpace(RDF, 2*N)
	M = G.random_element()
	Q, R = M.QR()
	# We compute the real erro by "rotating" the simples error.
	real_vector = Q*simple_error

	# Finally, we convert the error on RR^N to CC^N.
	return vector(CC, N, [(real_vector[i+1] + I * real_vector[i]) for i in range(0, 2*N, 2)])

def ParBW(p, s):
	N = len(s)
	if p < 4 or N == 1:
		return SeqBW(0, s)

	terms = list(s)
	s0 = vector(GQ, terms[0:N/2])
	s1 = vector(GQ, terms[N/2:N])
	PHI_2 = GQ(PHI/2)
	sm = PHI_2 * (s0 - s1)
	sp = PHI_2 * (s0 + s1)

	z0 = ParBW(p/4, s0)
	z1 = ParBW(p/4, s1)
	zm = ParBW(p/4, sm)
	zp = ParBW(p/4, sp)

	INVERSE_PHI = PHI.inverse()
	z0m = vector(GG, N, list(z0) + list(z0 - GG(2*INVERSE_PHI) * zm));
	z0p = vector(GG, N, list(z0) + list(GG(2*INVERSE_PHI)*zp - z0))
	z1m = vector(GG, N, list(GG(2*INVERSE_PHI)*zm + z1) + list(z1))
	z1p = vector(GG, N, list(GG(2*INVERSE_PHI)*zp - z1) + list(z1))

	arr = [z0m, z0p, z1m, z1p]
	return min(((s - v).norm(2), v) for v in arr)[1]

def RoundZZ(r):
	return ZZ(round(r))

def RoundGG(c):
	return GG(RoundZZ(c[0]) + I * RoundZZ(c[1]))

def Frac(r):
	return QQ(abs(r - RoundZZ(r)))

def R(c):
	return c[0]

def Img(c):
	return c[1]

def SeqBW(r, s):
	N = len(s)

	if N <= 2^r:
		return vector(GG, [RoundGG(t) for t in s])

	b = vector(GF(2), [GF(2)(RoundZZ(R(t)) + RoundZZ(Img(t))) for t in s])
	rho = vector(QQ, [QQ(1 - 2*max(Frac(R(t)), Frac(Img(t)))) for t in s])

	t = list([(b[i], rho[i]) for i in range(N)])
	PSIc = RMDec(r, t)

	INVERSE_PHI = GQ(PHI.inverse())
	v = SeqBW(r + 1, (s - PSIc) * INVERSE_PHI)
	return PSIc + (PHI * v)

def B(t):
	return t[0]

def Rho(t):
	return t[1]

def Eval(exp):
	return GF(2)(1) if exp else GF(2)(0)

def RMDec(r, t):
	# print("RMDec", t, r)
	N = len(t)
	if r == 0 or N == 2^r:
		return vector(ZZ, N, [ZZ(B(x)) for x in t])
	else:
		t0 = t[0:N/2]
		t1 = t[N/2:N]
		tp = [0]*int(N/2)
		for j in range(N/2):
			tp[j] = (B(t0[j]) + B(t1[j]), min(Rho(t0[j]), Rho(t1[j])))
		v = RMDec(r-1, tp)

		tm = [0]*int(N/2)
		for j in range(N/2):
			if (B(t0[j]) + B(t1[j])) == (v[j] % 2):
				tm[j] = (B(t0[j]), QQ((Rho(t0[j]) + Rho(t1[j])) / 2))
			else:
				tm[j] = (B(t0[j]) + Eval(Rho(t0[j]) < Rho(t1[j])), QQ(abs(Rho(t0[j]) - Rho(t1[j]))/2))

		u = RMDec(r, tm)
		return vector(ZZ, list(u) + list(u + v))
