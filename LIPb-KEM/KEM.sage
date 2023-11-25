import os
from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler
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

os.system('sage --preparse LaTeX.sage')
os.system('mv LaTeX.sage.py LaTeX.py')
import LaTeX as LaTeX


class LatticePair:

	def __init__(self, B, shortest_vector_length, decoding_distance, g):
		self.shortest_vector_length = g * shortest_vector_length
		self.decoding_distance = g * decoding_distance
		self.BS, self.BQ, self.g = self.__build_pair(B, g)
		self.s, self.q = self.__compute_parameters(self.BS, self.shortest_vector_length, self.decoding_distance)
		self.SIS = self.__compute_SIS_parameters(self.BS, self.q, self.decoding_distance)

	def __build_pair(self, B, g):
		BS = block_matrix([[g * B,			 0],
						   [	0, (g + 1) * B]])
		BQ = block_matrix([[B,				 0],
						   [0,   g*(g + 1) * B]])

		return (BS, BQ, g)

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
		M = ZZ(ceil(sqrt(N)/2))
		#M = ZZ(ceil(sqrt(N)/1.5))

		# extractor function moduli (TODO)
		# q_ = next_prime(ceil(2*q*rho * sqrt(n * log(n))))

		# THIS IS THE RIGHT ONE
		moduli = previous_prime(ceil(2^(N/(2*M)) / 3))
		#moduli = next_prime(floor(2^((N + log(M, 2))/(2*M))))

		beta = 2 * q * decoding_distance
		return (M, N, moduli, beta)

	def PublicKeySize(self):
		N, N = self.BS.dimensions()
		# This is the maximum length of the vectors in the generating set S
		lenS = self.s * sqrt(N)

		inBits = 0
		for i in range(N):
			for j in range(i+1):
				# These are the length of the i-th vector and j-th vector in the basis R that generates the Gram Matrix of the public key
				lenA = max(sqrt(i)/2 * lenS, lenS)
				lenB = max(sqrt(j)/2 * lenS, lenS)
				# By definition of Gram Matrix, G(i, j) = < R(i), R(j) > <= |R(i)| * |R(j)|
				inBits += ZZ(ceil(log(lenA * lenB, 2)))
		return inBits

	def SecretKeySize(self):
		N, N = self.BS.dimensions()
		# This is the maximum length of the vectors in the generating set S
		lenS = self.s * sqrt(N)

		inBits = N * ZZ(ceil(log(lenS, 2)))
		return inBits

	def EncapsulatedKeySize(self, seed_length=256):
		N, N = self.BS.dimensions()
		inBits = N * ZZ(ceil(log(self.q, 2))) + seed_length
		return inBits

	def EstimateBDDHardness(self):
		N, N = self.BS.dimensions()
		det = self.LargestDeterminantInSUVPReduction(self.BS, self.decoding_distance)
		beta = self.EstimateBetaForSUVP(det, N + 1)
		coreSVPHardness = log(beta, 2) + 0.292 * beta
		return beta, coreSVPHardness
		#print("BETA", beta)

	def KeySize(self):
		SISm, _, SISq, _ = self.SIS
		print(self.SIS)
		print("k =", SISm * ceil(log(SISq, 2)))
		return SISm * ceil(log(SISq, 2)) 

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

	def SmoothingParameterAssumption(self):
		N = self.BS.dimensions()[0]
		return self.decoding_distance / (2 * sqrt(N))

	def Compare(self, a, b, epsilon):
		if abs(a - b) <= epsilon:
			return 0

		if a < b:
			return -1

		return 1

	def SanityCheck(self, fail_on_assert = True):
		SISm, SISn, SISq, SISbeta = self.SIS

		if fail_on_assert:
			"""
				There are two main restrictions for the randomness extractor.
				(i) SISq must be greater than 2 * SISbeta, so that the proof of universial hash function holds.
				(ii) SISq must be a prime.

				More presisely, every integer in {-2 * beta, ..., -1, 1, ..., 2 * beta} must have an inverse (mod SISq)
			"""
			assert SISq > SISbeta
			assert SISq in Primes()

			n, n = self.BS.dimensions()
			"""
			For the statistical distance from the random bits to the uniform distribution to be negligible,
			we require SISn >= 2*SISM * log(SISq, 2)
			"""
			assert SISn >= 2 * SISm * log(SISq, 2)

			"""
			The key secretly shared through the encapsulation mechanism must be smaller than r - log2(3), where r is the rank
			of the underlying Barnes-Wall latttice
			"""
			r = n/2
			k = self.KeySize()
			assert k <= r - log(3, 2)
			return True

		try:
			return self.SanityCheck(fail_on_assert=True)
		except AssertionError:
			return False

class KEM:

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
		B = self.latticePair.BS
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
		BS = self.latticePair.BS
		q = self.latticePair.q
		g = self.latticePair.g
		rho = self.latticePair.decoding_distance
		SISm, _, SISq, _ = self.latticePair.SIS

		# Extract secret uniform transformation
		P, U = KEM.Extract(BS, S)

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

		Q = BS.transpose() * BS
		P = U.transpose() * Q * U
		
		determinant = sqrt(P.det())
		determinant = Ring(determinant)
		k = ((A * P) * x_) / determinant

		# Return key
		return k

SanityCheck = True

latex_security = LaTeX.SecurityOutput()
latex_parameters = LaTeX.ParameterOutput()

securityCaption = "Comparison of security levels and key sizes between our KEM and NewHope, where $| \\pk |$ is the public key size, $| \\sk |$ is the secret key size, $| \\textit{ch} |$ is the encapsulated key size, $\\beta$ is the blocksize, and $\\lambda$ is the number of operations needed to break the scheme using classical algorithms."
parameterCaption = "Suggested parameter sets for our concrete KEM based on LIP."

latex_security.BeginTable(securityCaption)
latex_security.MakeHeader()

latex_parameters.BeginTable(parameterCaption)
latex_parameters.MakeHeader()


dimensions = [512, 1024]
gs = [2, 4, 8, 12]

import logging

for g in gs:

	for dim in dimensions:
		logging.warning("Begin computing parameter set (" + str(dim) + ", " + str(g) + ")")
		i = int(log(dim/4, 2))

		B, shortest_vector_length, decoding_distance = BW.BarnesWall(i)
		B = BW.BasisOverZZ(B)
		latticePair = LatticePair(B, shortest_vector_length, decoding_distance, g)

		#dim, dim = latticePair.BS.dimensions()
		parameter_set = "OurBW" + str(int(dim/2)) + "g" + str(g)
		
		if SanityCheck:
			if not latticePair.SanityCheck(fail_on_assert=False):
				print(parameter_set + " failed SanityCheck")

		pk_in_bits = latticePair.PublicKeySize()
		sk_in_bits = latticePair.SecretKeySize()
		ch_in_bits = latticePair.EncapsulatedKeySize()
		k_in_bits = latticePair.KeySize()

		beta, security = latticePair.EstimateBDDHardness()
		if g >= 12:
			assumption = None
		else:
			assumption = latticePair.SmoothingParameterAssumption()
		
		BWdim, BWdim = B.dimensions()

		lattice = "\\BW_{" + str(BWdim) + "}"

		SISm, SISn, SISq, SISbeta = latticePair.SIS

		latex_security.MakeRow(parameter_set, dim, pk_in_bits, sk_in_bits, ch_in_bits, beta, security, assumption)
		latex_parameters.MakeRow(parameter_set, lattice, dim, g, latticePair.s, latticePair.q, SISq, SISm, SISbeta, k_in_bits)

		logging.warning("End computing parameter set (" + str(dim) + ", " + str(g) + ")")

latex_security.EndTable()
latex_security.Print()

latex_parameters.EndTable()
latex_parameters.Print()

# B, shortest_vector_length, decoding_distance = BW.BarnesWall(4)
# B = BW.BasisOverZZ(B)
# latticePair = LatticePair(B, shortest_vector_length, decoding_distance)
# N, N = latticePair.BQ.dimensions()
# kem = KEM(latticePair)
# while True:
# 	pk, sk = kem.GenerateKeyPair()
# 	ch, k = kem.EncapsulateKey(pk)
# 	k2 = kem.DecapsulateKey(sk, ch)
# 	assert k == k2

#trace("LatticePair(B, shortest_vector_length, decoding_distance)")			 # not tested
#time LatticePair(B, shortest_vector_length, decoding_distance)

# N, N = latticePair.BQ.dimensions()
# print("N", N)
# latticePair.PublicKeySize()
# latticePair.SecretKeySize()
# latticePair.EncapsulatedKeySize()
# latticePair.EstimateBDDHardness()

# kem = KEM(latticePair)
# pk, sk = kem.GenerateKeyPair()
# print("pk")
# print(pk)
# print("sk")
# print(sk)
# ch, k = kem.EncapsulateKey(pk)
# k = kem.DecapsulateKey(sk, ch)
# print("+------------------------------------_+")



# s = latticePair.s
# N, N = B.dimensions()
# S = matrix(ZZ, N, N, [[ 1,  1,  0, -1],
# 					  [ 1,  0, -1,  1],
# 					  [ 0, -1,  1,  1],
# 					  [ 0, -1,  1, -1]])
# KEM.Extract(B, S)

# S = matrix(ZZ, N, N, [[-2,-3,-2, 3],
# 					  [ 0, 1, 4,-1],
# 					  [-2, 2, 0,-2],
# 					  [ 0, 0, 0,-4]])
# KEM.Extract(B, S)

# S = matrix(ZZ, N, N, [[ -9, 11, -5, -9],
# 					  [ -6,  8, 18,-13],
# 					  [  3,  6, -7, 17],
# 					  [ -3,  2,  7, 11]])

# KEM.Extract(B, S)
# for _ in range(20000):
#  	R = KEM.SampleQuadraticForm(latticePair.BQ, s)

# kem = KEM(latticePair)

# print(latticePair.BQ)

# pk, sk = kem.GenerateKeyPair()
# encap, k = kem.EncapsulateKey(pk)
# k = kem.DecapsulateKey(sk, encap)
