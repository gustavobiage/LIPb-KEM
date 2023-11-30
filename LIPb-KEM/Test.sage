os.system('sage --preparse Schemes.sage')
os.system('mv Schemes.sage.py Schemes.py')
from Schemes import ParameterSet, IND_CPA_KEM, IND_CPA_PKE, IND_CCA_KEM

os.system('sage --preparse BW.sage')
os.system('mv BW.sage.py BW.py')
import BW as BarnesWall

os.system('sage --preparse Utils.sage')
os.system('mv Utils.sage.py Utils.py')
import Utils as Utils

os.system('sage --preparse LaTeX.sage')
os.system('mv LaTeX.sage.py LaTeX.py')
import LaTeX as LaTeX

import sys
import pickle
import logging

def GetFilename(dim, g):
	BWdim = dim // 2
	parameter = "OurBW" + str(BWdim) + "g" + str(g)
	return parameter + ".par"

def GetParameterSet(dim, g):
	filename = GetFilename(dim, g)

	if os.path.exists(filename):
		logging.warning("Obtaining local parameter set")
		with open(filename, "rb") as file:
			parameter_set = pickle.load(file)
	else:
		logging.warning("Computing parameter set from scratch")
		i = int(log(dim/4, 2))
		B, shortest_vector_length, decoding_distance = BarnesWall.BarnesWall(i)
		B = BarnesWall.BasisOverZZ(B)
		parameter_set = ParameterSet(B, shortest_vector_length, decoding_distance, g)
		SaveParameterSet(parameter_set)

	return parameter_set

def SaveParameterSet(parameter_set):
	dim, dim = parameter_set.BS.dimensions()
	g = parameter_set.g

	filename = GetFilename(dim, g)

	if not os.path.exists(filename):
		with open(filename, "wb") as file:
			pickle.dump(parameter_set, file)

def GetIND_CPA_KEM(dim, g):
	parameter_set = GetParameterSet(dim, g)
	return IND_CPA_KEM(parameter_set)

def GetIND_CPA_PKE(dim, g):
	ind_cpa_kem = GetIND_CPA_KEM(dim, g)
	return IND_CPA_PKE(ind_cpa_kem)

def GetIND_CCA2_KEM(dim, g):
	ind_cpa_pke = GetIND_CPA_PKE(dim, g)
	return IND_CCA_KEM(ind_cpa_pke)

def TestKEM(kem, TRIES=100):
	for _ in range(TRIES):
		pk, sk = kem.GenerateKeyPair()
		cc, k = kem.EncapsulateKey(pk)
		k_ = kem.DecapsulateKey(sk, cc)
		assert k == k_
		logging.warning("Success on decapsulation!")

def TestPKE(pke, TRIES=100, MESSAGE_LENGTH=3*128 + 64, BLOCK_SIZE=128):
	assert BLOCK_SIZE % 8 == 0
	
	for _ in range(TRIES):
		message = int(randint(0, 2^MESSAGE_LENGTH)).to_bytes(MESSAGE_LENGTH/8, byteorder='big')
		
		bytes_ = BLOCK_SIZE / 8
		# Add padding, if necessary
		m = message.ljust(ceil(len(message)/bytes_) * bytes_, b'\0')

		pk, sk = pke.GenerateKeyPair()
		ch = pke.Encrypt(pk, m)
		m_ = pke.Decrypt(sk, ch)
		assert m == m_
		logging.warning("Success on decryption!")

def TestIND_CPA_KEM(dim=2^6, g=2):
	ind_cpa_kem = GetIND_CPA_KEM(dim, g)
	TestKEM(ind_cpa_kem)

def TestIND_CPA_PKE(dim=2^5, g=2):
	ind_cpa_pke = GetIND_CPA_PKE(dim, g)
	TestPKE(ind_cpa_pke)

def TestIND_CCA2_KEM(dim=2^5, g=2):
	ind_cca2_kem = GetIND_CCA2_KEM(dim, g)
	TestKEM(ind_cca2_kem)

def TestBarnesWall(dim=2^8, g=2, TRIES=100):
	parameter_set = GetParameterSet(2*dim, g)
	logging.warning("Testing BW decoder of Micciancio and Nicolosi")

	for _ in range(TRIES):
		i = int(log(dim/2, 2))
		BW, shortest_vector_length, rho = BarnesWall.BarnesWall(i)
		N, N = BW.dimensions()

		BWZ = Utils.SortByColumnNorm(BarnesWall.BasisOverZZ(BW))

		s = parameter_set.s
		q = parameter_set.q

		standard_deviation = (q * rho) / sqrt(4*N)
		c = zero_vector(ZZ, 2 * N)

		error = BarnesWall.VectorOverCC(Utils.SampleD(BWZ, standard_deviation, c) / q)
		assert error.norm()^2 < rho * rho # Error must be decodable

		bwvector = BarnesWall.BWvector(BW) # Sample random vector in BW
		assert BarnesWall.isBWvector(BW, bwvector)

		decodable_vector = bwvector + error
		bwvector2 = BarnesWall.ParBW(1, decodable_vector)

		assert bwvector == bwvector2
		logging.warning("Success on decoding!")

def LaTeXTable(SANITY_CHECK=True):
	latex_security = LaTeX.SecurityOutput()
	latex_parameters = LaTeX.ParameterOutput()

	securityCaption = "Comparison of security levels and key sizes between our KEM and NewHope, where $| \\pk |$ is the public key size, $| \\sk |$ is the secret key size, $| \\textit{ch} |$ is the encapsulated key size, $\\beta$ is the blocksize, and $\\lambda$ is the number of operations needed to break the scheme using classical algorithms."
	parameterCaption = "Suggested parameter sets for our concrete KEM based on LIP."

	latex_security.BeginTable(securityCaption)
	latex_security.MakeHeader()

	latex_parameters.BeginTable(parameterCaption)
	latex_parameters.MakeHeader()

	dimensions = [32, 64, 128, 256, 512, 1024]
	gs = [2, 4, 8, 12]

	for g in gs:

		for dim in dimensions:
			pair = "(" + str(dim) + ", " + str(g) + ")"
			logging.warning("Begin computation on parameter set " + pair)
			
			BWdim = int(dim/2)

			parameter_set = GetParameterSet(dim, g)

			passes_sanity_check = not SANITY_CHECK or parameter_set.SanityCheck(fail_on_assert=False)

			if passes_sanity_check:
				logging.warning("Building KEM/PKE(s)")
				ind_cpa_kem = IND_CPA_KEM(parameter_set)
				ind_cpa_pke = IND_CPA_PKE(ind_cpa_kem)
				ind_cca2_kem = IND_CCA_KEM(ind_cpa_pke)

				for assymetric, class_ in [(ind_cpa_kem, "IND-CPA KEM"), (ind_cca2_kem, "IND-CCA2 KEM")]:
					parameter_name = "(" + class_ + ") OurBW" + str(BWdim) + "g" + str(g)
					
					logging.warning("Computing public key size")
					pk_in_bits = assymetric.PublicKeySize()
					logging.warning("Computing secret key size")
					sk_in_bits = assymetric.SecretKeySize()
					logging.warning("Computing encapsulated key size")
					ch_in_bits = assymetric.EncapsulatedKeySize()
					logging.warning("Computing shared secret size")
					k_in_bits = assymetric.SharedSecretSize()

					beta, security = (parameter_set.bkz_block, parameter_set.lambda_)
					if g >= 12:
						assumption = None
					else:
						assumption = parameter_set.SmoothingParameterAssumption()

					lattice = "\\BW_{" + str(BWdim) + "}"

					s = parameter_set.s
					q = parameter_set.q
					SISm, SISn, SISq, SISbeta = parameter_set.SIS

					latex_security.MakeRow(parameter_name, dim, pk_in_bits, sk_in_bits, ch_in_bits, beta, security, assumption)
					latex_parameters.MakeRow(parameter_name, lattice, dim, g, s, q, SISq, SISm, SISbeta, k_in_bits)
			else:
				logging.warning(pair + " failed SanityCheck!")

			logging.warning("End computation on parameter set " + pair)

	latex_security.EndTable()
	latex_security.Print()

	latex_parameters.EndTable()
	latex_parameters.Print()

exec_ = sys.argv[1] if len(sys.argv) >= 2 else "ALL"

if exec_ == "CPA_KEM":
	TestIND_CPA_KEM()
elif exec_ == "CPA_PKE":
	TestIND_CPA_PKE()
elif exec_ == "CCA2_KEM":
	TestIND_CCA2_KEM()
elif exec_ == "BW":
	TestBarnesWall()
elif exec_ == "LaTeX":
	LaTeXTable()
elif exec_ == "ALL":
	TestBarnesWall()
	TestIND_CPA_KEM()
	TestIND_CPA_PKE()
	TestIND_CCA2_KEM()