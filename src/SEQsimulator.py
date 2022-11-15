import random

def simulate_string(m):
    """Simulate a DNA sequence of length m.
    
    Example:
    simulate_string(20)
    
    """
    DNA = 'ACGT'
    nucleotides = [random.choice(DNA) for _ in range(m)]
    return ''.join(nucleotides)


def get_exact_read(DNA_seq:str, n:int):
    """Simulates an exact read of length n from a DNA sequence.
    
    Example:
    get_exact_read('GCTGAATCTATTTATGGTGG',5)

    """
    assert n < len(DNA_seq), 'm must be smaller or equal to length of DNA_seq'
    idx = random.randint(0,len(DNA_seq)-n)
    read = DNA_seq[idx:idx+n]
    return read
