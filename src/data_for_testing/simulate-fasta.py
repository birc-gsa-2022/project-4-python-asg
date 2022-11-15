"""
Program for simulating fasta files.
"""

import random
import argparse


DNA = 'ACGT'
LINE_WIDTH = 60


def simulate_string(m):
	"""Simulate a DNA sequence of length m."""
	nucleotides = [random.choice(DNA) for _ in range(m)]
	lines = []
	for i in range(0, m, LINE_WIDTH):
		lines.append(''.join(nucleotides[i:(i+LINE_WIDTH)]))
	return '\n'.join(lines)

def simulate_sequences(n, m):
	"""Simulate n sequences of length m."""
	for i in range(n):
		print("> seq{}".format(i))
		print(simulate_string(m))


if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='simulate-fasta', usage='%(prog)s [options]',
									 description="Simulate random sequences and output them in FASTA format.")
	parser.add_argument('-n', nargs='?', type=int, help='Number of sequences, default 1', default=1)
	parser.add_argument('-m', nargs='?', type=int, help='The length of the sequences, default 500', default=500)

	args = parser.parse_args()

	simulate_sequences(args.n, args.m)
