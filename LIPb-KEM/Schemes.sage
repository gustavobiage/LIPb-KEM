import os
from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler
from sage.crypto.mq.rijndael_gf import RijndaelGF
import sage.matrix.matrix_integer_dense_hnf as HNF

from sage.misc.trace import trace

os.system('sage --preparse Utils.sage')
os.system('mv Utils.sage.py Utils.py')
import Utils as Utils

os.system('sage --preparse QFUtils.sage')
os.system('mv QFUtils.sage.py QFUtils.py')
import QFUtils as QFUtils

os.system('sage --preparse BW.sage')
os.system('mv BW.sage.py BW.py')
import BW as BW

class LatticePair:

	def __init__(self, B, shortest_vector_length, decoding_distance):
		self.shortest_vector_length = shortest_vector_length
		self.decoding_distance = decoding_distance
		self.BQ, self.BS, self.g = self.__build_pair(B, shortest_vector_length)
		self.s, self.q = self.__compute_parameters(self.BQ, shortest_vector_length, decoding_distance)
		self.SIS = self.__compute_SIS_parameters(self.BQ, self.q, decoding_distance)

	def __build_pair(self, B, shortest_vector_length):
		gap = Utils.IntegralLatticeGap(B, shortest_vector_length)
		g = ZZ(ceil(min(gap^2, 16 * sqrt(2))))

		BQ = block_matrix([[B,				 0],
						   [0,   g*(g + 1) * B]])
		BS = block_matrix([[g * B,			 0],
						   [	0, (g + 1) * B]])

		return (BQ, BS, g)

	def __compute_parameters(self, B, shortest_vector_length, decoding_distance):
		N, N = B.dimensions()
		S = B.transpose() * B
		BSstar = Utils.GramSchmidt(B)

		s_ = max(v.norm() for v in BSstar) * sqrt(ln(2 * N + 4)/pi)
		s = ZZ(ceil(max(shortest_vector_length, s_)))

		q_ = (s * N) / decoding_distance * sqrt(ln(2 * N + 4) / pi)
		q = ZZ(ceil(q_))

		return (s, q)

	def __compute_SIS_parameters(self, B, q, decoding_distance):
		N, N = B.dimensions()

		# extractor function dimension 'm'. Values were chose so that 'l' is 'theta(n)'' bouded.
		M = ZZ(ceil(sqrt(N)))

		# extractor function moduli (TODO)
		# q_ = next_prime(ceil(2*q*rho * sqrt(n * log(n))))
		moduli = next_prime(floor(2^(N/M + M)))

		beta = 2 * q * decoding_distance
		return (M, N, moduli, beta)

	def PublicKeySize(self):
		N, N = self.BQ.dimensions()
		# This is the maximum length of the vectors in the generating set S
		lenS = self.s * sqrt(N)

		inBytes = 0
		for i in range(N):
			for j in range(i+1):
				# These are the length of the i-th vector and j-th vector in the basis R that generates the Gram Matrix of the public key
				lenA = max(sqrt(i)/2 * lenS, lenS)
				lenB = max(sqrt(j)/2 * lenS, lenS)
				# By definition of Gram Matrix, G(i, j) = < R(i), R(j) > <= |R(i)| * |R(j)|
				inBytes += ZZ(ceil(log(lenA * lenB, 2)))

		print("PK", inBytes, "b")
		inBytes = inBytes / 8
		print("PK", inBytes, "B")
		inBytes = inBytes / 1024
		print("PK", inBytes, "KB")
		inBytes = inBytes / 1024
		print("PK", inBytes, "MB")
		inBytes = inBytes / 1024
		print("PK", inBytes, "GB")
		return inBytes

	def SecretKeySize(self):
		N, N = self.BQ.dimensions()
		# This is the maximum length of the vectors in the generating set S
		lenS = self.s * sqrt(N)

		inBytes = N * ZZ(ceil(log(lenS, 2)))
		print("SK", inBytes, "b")
		inBytes = inBytes / 8
		print("SK", inBytes, "B")
		inBytes = inBytes / 1024
		print("SK", inBytes, "KB")
		inBytes = inBytes / 1024
		print("SK", inBytes, "MB")
		inBytes = inBytes / 1024
		print("SK", inBytes, "GB")

	def EncapsulatedKeySize(self):
		N, N = self.BQ.dimensions()
		# This is the maximum length of the vectors in the generating set S

		inBytes = N * ZZ(ceil(log(self.q, 2))) + 256
		print("CH", inBytes, "b")
		inBytes = inBytes / 8
		print("CH", inBytes, "B")
		inBytes = inBytes / 1024
		print("CH", inBytes, "KB")
		inBytes = inBytes / 1024
		print("CH", inBytes, "MB")
		inBytes = inBytes / 1024
		print("CH", inBytes, "GB")

	def EstimateBDDHardness(self):
		N, N = self.BQ.dimensions()
		det = self.LargestDeterminantInSUVPReduction(self.BQ, self.decoding_distance)
		beta = self.EstimateBetaForSUVP(det, N + 1)
		print("BETA", beta)

	def LargestDeterminantInSUVPReduction(self, B, rho):
		return B.det() * rho

	def EstimateBetaForSUVP(self, det, dim):
		for b in range(50, dim):
			delta = ((pi * b)^(1/b) *(b/(2 * pi * e)))^(1/(2*b - 2))
			eq1 = sqrt(b)
			eq2 = delta^(2*b - dim - 1) * det^(1 / dim)
			if self.Compare(eq1, eq2, 10e-7) < 1:
				return b
		return dim

	def Compare(self, a, b, epsilon):
		if abs(a - b) <= epsilon:
			return 0

		if a < b:
			return -1

		return 1

class IND_CPA_KEM:

	def __init__(self, latticePair):
		self.latticePair = latticePair

	def NearestPlane(B, t):
		B = Utils.LLL(B, 3/4) # Delta = 3/4
		N, M = B.dimensions()
		b = t
		for j in range(M-1, -1, -1):
			c = round(b.dot_product(B.column(j)) / Utils.Norm2(B.column(j)))
			b = b - c * B.column(j)
		return t - b

	def SampleQuadraticForm(B, s):
		print("SampleQuadraticForm")
		N, N  = B.dimensions()
		C = 1 - (1 + exp(-pi))^(-1)
		m = ceil(2*N / C)
		c = zero_vector(ZZ, N)

		#time S = matrix(ZZ, m, N, [KEM.SampleD(B, s, c) for _ in range(m)]).transpose()
		D = DiscreteGaussianDistributionLatticeSampler(B.transpose(), s, zero_vector(N))
		S = matrix(ZZ, m, N, [D() for _ in range(m)]).transpose()

		if S.rank() < N:
			return KEM.SampleQuadraticForm(B, s)

		S = KEM.Simplify(S)
		S = Utils.SortByColumnNorm(S)
		R, U = KEM.Extract(B, S)
		return (R, U, S)

	def Simplify(S):
		"""
			Let S be an set of N-dimensional vectors in L(B) ordered by norm
			and with rank N.

			We select the first N linearly independent of vectors from S such
			that |s1| <= |s2| <= ... <= |sN|.
		"""
		N, M = S.dimensions()
		assert N > 0
		S_ = matrix(ZZ, N, 0)
		
		for column in S.columns():
			aux = block_matrix([[S_, matrix(column).transpose()]])

			if aux.rank() > S_.rank():
				S_ = aux

			if S_.rank() == N:
				return S_

		assert S.rank() == N

	def Extract(B, S):
		print("Extract")
		N, M = S.dimensions()

		Q = matrix(ZZ, B.inverse() * S)
		T, U_ = HNF.hnf_with_transformation(Q)

		U = matrix(ZZ, U_.inverse())
		R = B * U

		assert S == R * T
		S_ = Utils.GramSchmidt(S)

		for k in range(N):
			s = S.column(k)
			r = R.column(k)

			if Utils.Norm2(r) <= max((k+1)/4, 1) * Utils.Norm2(s):
				continue

			R_ = Utils.GramSchmidt(R)
			assert all(R_.column(i).norm() <= S_.column(i).norm() for i in range(N))

			r_ = R_.column(k)
			s_ = S_.column(k)

			if r_ == s_ or r_ == -s_:
				Utils.InsertColumn(R, k, s)
			else:
				x = r - r_ # Projection
				BR = block_matrix([[matrix(R.column(i)).transpose() for i in range(k)]])

				v = KEM.NearestPlane(BR, x)
				assert Utils.Norm2(x - v) <= (N-1) * max(Utils.Norm2(s) for s in S.columns())/4
				assert Utils.Norm2(r - v) <= max(((k+1)/4), 1) * Utils.Norm2(s)

				Utils.InsertColumn(R, k, r - v)

		# Since we modified vectors in R, we must recompute U such that R = B * U
		U = matrix(ZZ, N, N, B.inverse() * R)
		return (R.transpose() * R, U)

	def GenerateKeyPair(self):
		print("GENERATE")
		B = self.latticePair.BQ
		s = self.latticePair.s
		R, U, S = KEM.SampleQuadraticForm(B, s)
		assert U.det() == 1 or U.det() == -1
		assert R.rank() == B.rank()
		return (R, S)

	def SampleA(seed, N, M, Ring):
		set_random_seed(seed)
		M = random_matrix(Ring, N, M)
		set_random_seed()
		return M

	def EncapsulateKey(self, public_key):
		print("ENCAPSULATE")
		P = public_key
		q = self.latticePair.q
		rho = self.latticePair.decoding_distance
		SISm, _, SISq, _ = self.latticePair.SIS

		N, N = P.dimensions()
		s_ = (q * rho) / sqrt(N)
		e_ = QFUtils.SampleD(P, s_)
		e = e_ / q
		assert QFUtils.Norm2(P, e) <= rho*rho

		c = vector(QQ, N, [QQ((x.numerator() % x.denominator()) / x.denominator()) for x in e])

		seed = randint(0, 256)

		Ring = IntegerModRing(SISq)
		A = KEM.SampleA(seed, SISm, N, Ring)

		# We must multiply the result by the inverse modular of
		# the determinant of BP^T (mod SISq).
		determinant = sqrt(P.det())
		determinant = Ring(determinant)
		k = ((A * P) * e_) / determinant

		return ((c, seed), k)

	def DecapsulateKey(self, secret_key, encapsulatedKey):
		print("DECAPSULATE")
		c, seed = encapsulatedKey
		S = secret_key

		# Utils information
		BQ = self.latticePair.BQ
		q = self.latticePair.q
		g = self.latticePair.g
		rho = self.latticePair.decoding_distance
		SISm, _, SISq, _ = self.latticePair.SIS

		# Extract secret uniform transformation
		P, U = KEM.Extract(BQ, S)

		# U*c returns the coefficients of the vector we must decode.
		# By multiplying it by BQ, we obtain the actual vector.
		v = (BQ * U) * c

		# Since BW is implemented using Gaussian Integers, we must convert
		# the vectors back and forward between ZZ and GG
		# (Analogue, error vectors contain rationals and must be converted to CC)
		v1 = BW.VectorOverCC(v[0:N/2])
		v2 = BW.VectorOverCC(vector(v[N/2:N]) / (g * (g + 1)))


		y1 = BW.VectorOverZZ(BW.SeqBW(0, v1))
		y2 = BW.VectorOverZZ(BW.SeqBW(0, v2))

		# Compute decoded vector, as defined by the lattice pair.
		y = vector(ZZ, N, list(y1) + list(y2 * (g * (g + 1))))
		
		assert Utils.isLatticeVector(BQ, y)

		# Reconstruct encapsulated key from extractor
		x = c - (U.inverse() * Utils.Coeff(BQ, y))
		x_ = vector(ZZ, q * x)
		assert Utils.Norm2(y - v) <= rho*rho

		Ring = IntegerModRing(SISq)
		A = KEM.SampleA(seed, SISm, N, Ring)

		Q = BQ.transpose() * BQ
		P = U.transpose() * Q * U
		
		determinant = sqrt(P.det())
		determinant = Ring(determinant)
		k = ((A * P) * x_) / determinant

		# Return key
		return k


class IND_CPA_PKE:

	def __init__(self, IND_CPA_KEM):
		self.IND_CPA_KEM = IND_CPA_KEM

	def GenerateKeyPair(self):
		return self.IND_CPA_KEM.GenerateKeyPair()

	def Encrypt(self, pk, message):
		encap, k = self.IND_CPA_KEM.EncapsulateKey(pk)


B, shortest_vector_length, decoding_distance = BW.BarnesWall(4)
B = BW.BasisOverZZ(B)
latticePair = LatticePair(B, shortest_vector_length, decoding_distance)
N, N = latticePair.BQ.dimensions()
kem = KEM(latticePair)
while True:
	pk, sk = kem.GenerateKeyPair()
	ch, k = kem.EncapsulateKey(pk)
	k2 = kem.DecapsulateKey(sk, ch)
	assert k == k2
