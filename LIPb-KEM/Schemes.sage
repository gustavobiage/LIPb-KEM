import os
from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler
import sage.matrix.matrix_integer_dense_hnf as HNF
from Crypto.Cipher import AES
from Crypto.Hash import SHAKE256
from sage.misc.trace import trace

os.system('sage --preparse Utils.sage')
os.system('mv Utils.sage.py Utils.py')
import Utils as Utils

os.system('sage --preparse AES.sage')
os.system('mv AES.sage.py AES.py')
from AES import AESEncrypt, AESDecrypt

os.system('sage --preparse QFUtils.sage')
os.system('mv QFUtils.sage.py QFUtils.py')
import QFUtils as QFUtils

os.system('sage --preparse BW.sage')
os.system('mv BW.sage.py BW.py')
import BW as BW

os.system('sage --preparse LaTeX.sage')
os.system('mv LaTeX.sage.py LaTeX.py')
import LaTeX as LaTeX

os.system('sage --preparse KEM.sage')
os.system('mv KEM.sage.py KEM.py')
from KEM import AbstractKEM, AbstractPKE

class ParameterSet:

	def __init__(self, B, shortest_vector_length, decoding_distance, g, security_estimation = None):
		self.shortest_vector_length = g * shortest_vector_length
		self.decoding_distance = g * decoding_distance
		self.BS, self.BQ, self.g = ParameterSet.__BuildLatticePair(B, g)

		self.s, self.q = ParameterSet.__ComputeFrameworkParameters(self.BS, self.shortest_vector_length, self.decoding_distance)
		self.SIS = ParameterSet.__ComputeSISParameters(self.BS, self.q, self.decoding_distance)
		self.bkz_block, self.lambda_ = security_estimation if security_estimation else self.EstimateBDDHardness()

	def __BuildLatticePair(B, g):
		BS = block_matrix([[g * B,			 0],
						   [	0, (g + 1) * B]], subdivide=False)
		BQ = block_matrix([[B,				 0],
						   [0,   g*(g + 1) * B]], subdivide=False)

		return (BS, BQ, g)

	def __ComputeFrameworkParameters(B, shortest_vector_length, decoding_distance):
		N, N = B.dimensions()
		S = B.transpose() * B
		BSstar = Utils.GramSchmidt(B)

		s_ = max(v.norm() for v in BSstar) * sqrt(ln(2 * N + 4)/pi)
		s = ZZ(ceil(max(shortest_vector_length, s_)))

		q_ = (s * N) / decoding_distance * sqrt(ln(2 * N + 4) / pi)
		q = ZZ(ceil(q_))

		return (s, q)

	def __ComputeSISParameters(B, q, decoding_distance):
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

	def EstimateBDDHardness(self, C=0.292):
		N, N = self.BS.dimensions()
		det = self.LargestDeterminantInSUVPReduction(self.BS, self.decoding_distance)
		beta = self.EstimateBetaForSUVP(det, N + 1)
		coreSVPHardness = log(beta, 2) + C * beta
		return beta, floor(coreSVPHardness)
		#print("BETA", beta)

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
			k = SISm * ceil(log(SISq, 2))
			assert k <= r - log(3, 2)
			return True

		try:
			return self.SanityCheck(fail_on_assert=True)
		except AssertionError:
			return False

class IND_CPA_KEM(AbstractKEM):

	def __init__(self, parameter_set):
		super().__init__()
		self.parameter_set = parameter_set

	def GetParameterSet(self):
		return self.parameter_set

	#@property
	def forcePublicKeySize(self):
		parameter_set = self.GetParameterSet()

		N, N = parameter_set.BS.dimensions()
		# This is the maximum length of the vectors in the generating set S
		lenS = parameter_set.s * sqrt(N)

		inBits = 0
		for i in range(N):
			for j in range(i+1):
				# These are the length of the i-th vector and j-th vector in the basis R that generates the Gram Matrix of the public key
				lenA = max(sqrt(i)/2 * lenS, lenS)
				lenB = max(sqrt(j)/2 * lenS, lenS)
				# By definition of Gram Matrix, G(i, j) = < R(i), R(j) > <= |R(i)| * |R(j)|
				inBits += ZZ(ceil(log(lenA * lenB, 2)) + 1)
		return inBits

	#@property
	def forceSecretKeySize(self):
		parameter_set = self.GetParameterSet()

		N, N = parameter_set.BS.dimensions()
		s = parameter_set.s

		# This is the maximum length of the vectors in the generating set S
		lenS = s * sqrt(N)

		inBits = N * N * ZZ(ceil(log(lenS, 2)) + 1)
		return inBits

	#@property
	def forceEncapsulatedKeySize(self):
		parameter_set = self.GetParameterSet()
		
		q = parameter_set.q
		N, N = parameter_set.BS.dimensions()
		seed_length = parameter_set.lambda_

		inBits = N * ZZ(ceil(log(q, 2))) + seed_length
		return inBits
	
	#@property
	def forceSharedSecretSize(self):
		parameter_set = self.GetParameterSet()
		m, _, q, _ = parameter_set.SIS
		return m * ceil(log(q, 2)) 

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

		#time S = matrix(ZZ, m, N, [IND_CPA_KEM.SampleD(B, s, c) for _ in range(m)]).transpose()
		D = DiscreteGaussianDistributionLatticeSampler(B.transpose(), s, zero_vector(N))
		S = matrix(ZZ, m, N, [D() for _ in range(m)]).transpose()

		if S.rank() < N:
			return IND_CPA_KEM.SampleQuadraticForm(B, s)

		S = IND_CPA_KEM.Simplify(S)
		S = Utils.SortByColumnNorm(S)
		R, U = IND_CPA_KEM.Extract(B, S)
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

				v = IND_CPA_KEM.NearestPlane(BR, x)
				assert Utils.Norm2(x - v) <= (N-1) * max(Utils.Norm2(s) for s in S.columns())/4
				assert Utils.Norm2(r - v) <= max(((k+1)/4), 1) * Utils.Norm2(s)

				Utils.InsertColumn(R, k, r - v)

		# Since we modified vectors in R, we must recompute U such that R = B * U
		U = matrix(ZZ, N, N, B.inverse() * R)
		return (R.transpose() * R, U)

	#@property
	def GenerateKeyPair(self):
		print("GENERATE")
		parameter_set = self.GetParameterSet()
		B = parameter_set.BS
		s = parameter_set.s
		R, U, S = IND_CPA_KEM.SampleQuadraticForm(B, s)
		assert U.det() == 1 or U.det() == -1
		assert R.rank() == B.rank()
		return (R, S)

	def SampleA(seed_, N, M, Ring):
		with seed(seed_):
			M = random_matrix(Ring, N, M)
		return M

	#@property
	def EncapsulateKey(self, public_key):
		print("ENCAPSULATE")
		P = public_key
		parameter_set = self.GetParameterSet()
		q = parameter_set.q
		rho = parameter_set.decoding_distance
		SISm, _, SISq, _ = parameter_set.SIS

		N, N = P.dimensions()
		s_ = (q * rho) / sqrt(N)
		e_ = QFUtils.SampleD(P, s_)
		e = e_ / q
		assert QFUtils.Norm2(P, e) <= rho*rho

		c = vector(QQ, N, [QQ((x.numerator() % x.denominator()) / x.denominator()) for x in e])

		seed_length = parameter_set.lambda_
		seed_ = randint(0, 2^seed_length)

		Ring = IntegerModRing(SISq)
		A = IND_CPA_KEM.SampleA(seed_, SISm, N, Ring)

		# We must multiply the result by the inverse modular of
		# the determinant of BP^T (mod SISq).
		determinant = sqrt(P.det())
		determinant = Ring(determinant)
		k = ((A * P) * e_) / determinant

		return ((c, seed_), k)

	#@property
	def DecapsulateKey(self, secret_key, encap):
		print("DECAPSULATE")
		c, seed_ = encap
		S = secret_key

		parameter_set = self.GetParameterSet()
		# Utils information
		N, N = parameter_set.BS.dimensions()
		BS = parameter_set.BS
		q = parameter_set.q
		g = parameter_set.g
		rho = parameter_set.decoding_distance
		SISm, _, SISq, _ = parameter_set.SIS

		# Extract secret uniform transformation
		P, U = IND_CPA_KEM.Extract(BS, S)

		# U*c returns the coefficients of the vector we must decode.
		# By multiplying it by BS, we obtain the actual vector.
		v = (BS * U) * c

		# Since BW is implemented using Gaussian Integers, we must convert
		# the vectors back and forward between ZZ and GG
		# (Analogue, error vectors contain rationals and must be converted to CC)
		v1 = BW.VectorOverCC(vector(v[0:N/2]) /  g     )
		v2 = BW.VectorOverCC(vector(v[N/2:N]) / (g + 1))


		y1 = BW.VectorOverZZ(BW.SeqBW(0, v1))
		y2 = BW.VectorOverZZ(BW.SeqBW(0, v2))

		# Compute decoded vector, as defined by the lattice pair.
		y = vector(ZZ, N, list(y1 * g) + list(y2 * (g + 1)))
		
		assert Utils.isLatticeVector(BS, y)

		# Reconstruct encapsulated key from extractor
		x = c - (U.inverse() * Utils.Coeff(BS, y))
		x_ = vector(ZZ, q * x)
		assert Utils.Norm2(y - v) <= rho*rho

		Ring = IntegerModRing(SISq)
		A = IND_CPA_KEM.SampleA(seed_, SISm, N, Ring)

		Q = BS.transpose() * BS
		P = U.transpose() * Q * U
		
		determinant = sqrt(P.det())
		determinant = Ring(determinant)
		k = ((A * P) * x_) / determinant

		# Return key
		return k

	#@property
	def EncodePublicKey(self, Q2):
		parameter_set = self.GetParameterSet()
		N, N = parameter_set.BS.dimensions()
		# This is the maximum length of the vectors in the generating set S
		lenS = parameter_set.s * sqrt(N)

		encoding = 0
		allBits = 0
		for i in range(N):
			for j in range(i+1):
				# These are the length of the i-th vector and j-th vector in the basis R that generates the Gram Matrix of the public key
				lenA = max(sqrt(i)/2 * lenS, lenS)
				lenB = max(sqrt(j)/2 * lenS, lenS)
				# By definition of Gram Matrix, G(i, j) = < R(i), R(j) > <= |R(i)| * |R(j)|
				bits = int(ceil(log(lenA * lenB, 2)) + 1)

				encoding = (int(encoding) << bits) | int(Q2[i,j] + lenA * lenB)
				allBits += bits

		return int(encoding).to_bytes(ceil(allBits/8), byteorder='big')

	#@property
	def EncodeSecretKey(self, sk):
		parameter_set = self.GetParameterSet()

		N, N = parameter_set.BS.dimensions()
		s = parameter_set.s

		# This is the maximum length of the vectors in the generating set S
		bits = log(s * sqrt(N), 2) + 1
		length = N * N * bits

		S = sk

		encoding = 0

		for vec in S.columns():
			for term in vec:
				encoding = (int(encoding) << bits) | int(term)

		return int(encoding).to_bytes(ceil(length/8), byteorder='big')

	#@property
	def EncodeEncapsulatedKey(self, cc):
		parameter_set = self.GetParameterSet()
		N, N, = parameter_set.BS.dimensions()
		q = parameter_set.q
		seed_length = parameter_set.lambda_

		bits = ceil(log(q, 2))
		length = N * bits + seed_length

		c, seed_ = cc

		encoded = 0

		for term in c:
			encoded = (int(encoded) << bits) | int(term)

		encoded = int(encoded).to_bytes(ceil(length/8), byteorder='big')
		encoded_seed = int(seed_).to_bytes(ceil(seed_length/8), byteorder='big')
		return encoded + encoded_seed

	#@property
	def EncodeSharedSecret(self, ss):
		parameter_set = self.GetParameterSet()
		_, _, q, _ = parameter_set.SIS

		bits = int(ceil(log(q, 2)))

		encoding = 0
		for term in ss:
			encoding = (int(encoding) << bits) | int(term)

		return encoding

class IND_CPA_PKE(AbstractPKE):

	def __init__(self, IND_CPA_KEM):
		super().__init__()
		self.IND_CPA_KEM = IND_CPA_KEM

	def GetParameterSet(self):
		return self.IND_CPA_KEM.GetParameterSet()

	def EncodeKey(self, k, left=None):
		keyLength = self.IND_CPA_KEM.SharedSecretSize()
		encoding = self.IND_CPA_KEM.EncodeSharedSecret(k)

		aes = [128, 192, 256, keyLength - (keyLength % 8)]
		aes = max(filter(lambda x : x <= keyLength, aes))

		if left:
			aes = min(aes, left)

		encoding = encoding % aes
		return int(encoding).to_bytes(aes//8, byteorder='big')

	#@property
	def forcePublicKeySize(self):
		return self.IND_CPA_KEM.PublicKeySize()

	#@property
	def forceSecretKeySize(self):
		return self.IND_CPA_KEM.SecretKeySize()

	def CiphertextSize(self, message_length, SYMMETRIC_KEY_LENGTH=128):
		return message_length + self.IND_CPA_KEM.EncapsulatedKeySize()

	#@property
	def GenerateKeyPair(self):
		return self.IND_CPA_KEM.GenerateKeyPair()

	#@property
	def Encrypt(self, pk, message, SYMMETRIC_KEY_LENGTH=128, seed_=None):
		# set random seed, if requested
		if seed_:
			set_random_seed(seed_)

		key = bytearray()
		encap = []

		while len(key)*8 < SYMMETRIC_KEY_LENGTH: 
			cc, k = self.IND_CPA_KEM.EncapsulateKey(pk)
			local_key = self.EncodeKey(k, left = SYMMETRIC_KEY_LENGTH - len(key)*8)
			key = key + local_key
			encap.append(cc)

		cipher = AES.new(key, AES.MODE_CBC, iv=b"0000000000000000")
		ciphertext = cipher.encrypt(message)

		return (ciphertext, encap)

	#@property
	def Decrypt(self, sk, ch):
		ciphertext, cc = ch

		key = bytearray()
		
		for c in cc:
			k = self.IND_CPA_KEM.DecapsulateKey(sk, c)
			local_key = self.EncodeKey(k, left = 128 - len(key)*8)
			key += local_key

		cipher = AES.new(key, AES.MODE_CBC, iv=b"0000000000000000")
		plain = cipher.decrypt(ciphertext)

		return plain;

	#@property
	def EncodePublicKey(self, pk):
		return self.IND_CPA_KEM.EncodePublicKey(pk)

	#@property
	def EncodeSecretKey(self, sk):
		return self.IND_CPA_KEM.EncodeSecretKey(sk)

	def EncodeCiphertext(self, ch):
		ciphertext, cc = ch
		encoded = ciphertext

		for c in cc:
			encoded = encoded + self.IND_CPA_KEM.EncodeEncapsulatedKey(c)

		return encoded

class IND_CCA_KEM(AbstractKEM):

	def __init__(self, IND_CPA_PKE):
		super().__init__()
		self.IND_CPA_PKE = IND_CPA_PKE

	def GetParameterSet(self):
		return self.IND_CPA_PKE.GetParameterSet()

	#@property
	def forcePublicKeySize(self):
		return self.IND_CPA_PKE.PublicKeySize()

	#@property
	def forceSecretKeySize(self):
		parameter_set = self.GetParameterSet()
		seed_length = parameter_set.lambda_
		return self.IND_CPA_PKE.SecretKeySize() + seed_length

	#@property
	def forceEncapsulatedKeySize(self, AES_BLOCK_SIZE=128):
		parameter_set = self.GetParameterSet()
		message_length = ceil(parameter_set.lambda_ / AES_BLOCK_SIZE) * AES_BLOCK_SIZE
		return self.IND_CPA_PKE.CiphertextSize(message_length)

	#@property
	def forceSharedSecretSize(self):
		parameter_set = self.GetParameterSet()
		return parameter_set.lambda_

	#@property
	def GenerateKeyPair(self):
		parameter_set = self.GetParameterSet()
		pk, sk = self.IND_CPA_PKE.GenerateKeyPair()
		seed_length = parameter_set.lambda_
		skr = randint(0, 2^seed_length)
		return (pk, (sk, skr))

	#@property
	def EncapsulateKey(self, pk, AES_BLOCK_SIZE=128):
		parameter_set = self.GetParameterSet()
		seed_length = parameter_set.lambda_

		# Find the next multiple of the aes blocksize
		message_length = ceil(parameter_set.lambda_ / AES_BLOCK_SIZE) * AES_BLOCK_SIZE
		m = randint(0, 2^message_length).to_bytes(ceil(message_length/8), byteorder='big')

		encoded_public_key = self.EncodePublicKey(pk)

		shake = SHAKE256.new()
		shake.update(encoded_public_key)
		shake.update(m)

		extracted_bytes = shake.read(2*seed_length)

		r = int.from_bytes(extracted_bytes[0:seed_length], byteorder='big')
		k = extracted_bytes[seed_length:]

		ch = self.IND_CPA_PKE.Encrypt(pk, m, seed_=r)
		encoded_ch = self.IND_CPA_PKE.EncodeCiphertext(ch)

		shake = SHAKE256.new()
		shake.update(encoded_ch)
		shake.update(k)

		s = shake.read(seed_length)

		return (ch, s)

	#@property
	def DecapsulateKey(self, secret_key, c):
		S, skr = secret_key
		m_ = self.IND_CPA_PKE.Decrypt(S, c)

		"""
			There is no need to pass the public key to decapsulate since the public key can be derived
			from the secret key.

			The secret key is 'S', which allows computing the transformation 'U' and the basis 'R'.
			The public key is the quadratic form U^T Q U.
		"""
		parameter_set = self.GetParameterSet()
		# Utils information
		BS = parameter_set.BS
		q = parameter_set.q
		g = parameter_set.g
		rho = parameter_set.decoding_distance
		SISm, _, SISq, _ = parameter_set.SIS
		seed_length = parameter_set.lambda_

		# Extract secret uniform transformation and the public key
		P, U = IND_CPA_KEM.Extract(BS, S)

		encoded_public_key = self.EncodePublicKey(P)
		shake = SHAKE256.new()
		shake.update(encoded_public_key)
		shake.update(m_)

		extracted_bytes = shake.read(2*seed_length)

		r_ = int.from_bytes(extracted_bytes[0:seed_length], byteorder='big')
		k_ = extracted_bytes[seed_length:]

		c_ = self.IND_CPA_PKE.Encrypt(P, m_, seed_=r_)
		encoded_c = self.IND_CPA_PKE.EncodeCiphertext(c)

		shake = SHAKE256.new()
		if c_ == c:
			shake.update(encoded_c)
			shake.update(k_)
		else:
			shale.update(encoded_c)
			shake.update(skr)

		return shake.read(seed_length)

	#@property
	def EncodePublicKey(self, pk):
		return self.IND_CPA_PKE.EncodePublicKey(pk)

	#@property
	def EncodeSecretKey(self, sk):
		S, skr = sk
		return self.IND_CPA_PKE.EncodeSecretKey(S) + skr

	#@property
	def EncodeEncapsulatedKey(self, cc):
		return self.IND_CPA_PKE.EncodeCiphertext(cc)

	#@property
	def EncodeSharedSecret(self, ss):
		return ss # Already a bytearray
