import abc
from abc import ABC, abstractmethod 

class Assymetric(ABC):
	def __init__(self):
		self.public_key_size = None
		self.secret_key_size = None

	def PublicKeySize(self):
		if self.public_key_size is None:
			self.public_key_size = self.forcePublicKeySize()

		return self.public_key_size

	#@abc.abstractproperty
	@abstractmethod
	def forcePublicKeySize(self):
		pass

	def SecretKeySize(self):
		if self.secret_key_size is None:
			self.secret_key_size = self.forceSecretKeySize()

		return self.secret_key_size

	@abstractmethod
	def forceSecretKeySize(self):
		pass

	@abstractmethod
	def GenerateKeyPair(self):
		pass

	@abstractmethod
	def EncodePublicKey(self, pk):
		pass

	@abstractmethod
	def EncodeSecretKey(self, sk):
		pass

class AbstractKEM(Assymetric):

	def __init__(self):
		super().__init__()
		self.encapsulated_key_size = None
		self.shared_secret_size = None

	def EncapsulatedKeySize(self):
		if self.encapsulated_key_size is None:
			self.encapsulated_key_size = self.forceEncapsulatedKeySize()

		return self.encapsulated_key_size

	@abstractmethod
	def forceEncapsulatedKeySize(self):
		pass

	def SharedSecretSize(self):
		if self.shared_secret_size is None:
			self.shared_secret_size = self.forceSharedSecretSize()

		return self.shared_secret_size

	@abstractmethod
	def forceSharedSecretSize(self):
		pass

	@abstractmethod
	def GenerateKeyPair(self):
		pass

	@abstractmethod
	def EncapsulateKey(self, pk):
		pass

	@abstractmethod
	def DecapsulateKey(self, sk, cc):
		pass

	@abstractmethod
	def EncodeEncapsulatedKey(self, cc):
		pass

class AbstractPKE(Assymetric):

	def __init__(self):
		super().__init__()

	@abstractmethod
	def Encrypt(self, pk, message, seed_=None):
		pass

	@abstractmethod
	def Decrypt(self, sk, ch):
		pass

	@abstractmethod
	def EncodeCiphertext(self, cc):
		pass

