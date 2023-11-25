
class Output:
	IDENT = "\t"
	ParameterTabular = ""

	def __init__(self):
		self.output = []

	def BeginTable(self, caption, tabular):
		header = []
		header.append("\\begin{table}[htbp]")
		header.append(Output.IDENT + "\\centering")
		header.append(Output.IDENT + "\\setlength{\\tabcolsep}{3pt}")
		header.append(Output.IDENT + "\\caption{" + caption + "}")
		header.append(Output.IDENT + "\\begin{tabular}{" + tabular + "}")
		header.append(2*Output.IDENT + "\\toprule")
		self.output += header

	def EndTable(self):
		footer = []
		footer.append(2*Output.IDENT + "\\bottomrule")
		footer.append(Output.IDENT + "\\end{tabular}")
		footer.append("\\end{table}")
		self.output += footer

	def Print(self):
		print("\n".join(self.output))

class ParameterOutput(Output):
	tabular = "crrrrrrrrr"

	def __init__(self):
		super().__init__()

	def BeginTable(self, caption):
		super().BeginTable(caption, ParameterOutput.tabular)

	def MakeHeader(self):
		header = 2*Output.IDENT + "Parameter set & $\\LL$ & $\\dim$ & $g$ & $\\sim s$ & $q$ & $\\hat{q}$ & $\\hat{m}$ & $\\beta$ & $|k|$ \\\\"
		self.output.append(header)
		self.output.append(2*Output.IDENT + "\\midrule")

	def MakeRow(self, parameter_set, lattice, dim, g, s, q, q_, m_, beta, k):
		row = 2*Output.IDENT + parameter_set + " & $" + lattice + "$ & $" + str(dim) + "$ & $" + str(g) + "$ & $" + str(s) + "$ & $" + str(q) + "$ & $" + str(q_) + "$ & $" + str(m_) + "$ & $" + str(beta) + "$ & $" + str(k) + "$ b" + " \\\\"
		self.output.append(row)

class SecurityOutput(Output):

	tabular = "crrrrrc"

	def __init__(self):
		super().__init__()

	def BeginTable(self, caption):
		super().BeginTable(caption, SecurityOutput.tabular)

	def MakeHeader(self):
		header = 2*Output.IDENT + "Parameter set & $|\\pk|$ & $|\\sk|$ & $|\\ch|$ & $\\beta$ & $\\lambda$ & Assumption \\\\"
		self.output.append(header)
		self.output.append(2*Output.IDENT + "\\midrule")

	def MakeRow(self, parameter_set, dim, pk, sk, ch, beta, security, assumption):
		if assumption:
			assumption = "$\\eta(\\BW_{" + str(int(dim/2)) + "}) \\leq " + str(assumption) + "$"
		else:
			assumption = "--"

		pk = self.FormatLength(pk)
		sk = self.FormatLength(sk)
		ch = self.FormatLength(ch)
		security = "2^{" + str(int(floor(security))) + "}"
		row = 2*Output.IDENT + parameter_set + " & " + pk + " & " + sk + " & " + ch + " & $" + str(beta) + "$ & $" + security + "$ & " + assumption + " \\\\"
		self.output.append(row)

	def FormatLength(self, length_in_bits, original_as_extension=False):
		metric = ["b", "B", "KB", "MB", "GB", "TB"]
		conversion = [8, 1024, 1024, 1024, 1024]
		idx = 0
		previous = length_in_bits
		next_ = length_in_bits / conversion[idx]
		
		while next_ >= 1:
			idx += 1
			previous = next_
			next_ = next_ / conversion[idx]

		if original_as_extension:
			extension = " ($" + str(length_in_bits) + "$ b)"
		else:
			extension = ""

		if previous not in ZZ:
			return "$\\sim " + "{:.2f}".format(float(previous)) + "$ " + metric[idx] + extension

		return "$" + str(previous) + "$ " + metric[idx] + extension
