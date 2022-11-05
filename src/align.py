"""A module for translating between alignments and edits sequences."""


def get_edits(p: str, q: str):
    """Extract the edit operations from a pairwise alignment.

    Args:
        p (str): The first row in the pairwise alignment.
        q (str): The second row in the pairwise alignment.

    Returns:
        str: The list of edit operations as a string.

    >>> get_edits('ACCACAGT-CATA', 'A-CAGAGTACAAA')
    ('ACCACAGTCATA', 'ACAGAGTACAAA', 'MDMMMMMMIMMMM')

    """
    assert len(p) == len(q)

    edits = ''  # edits to go from p to q
    for i in range(len(p)):
        if p[i] != '-' and q[i] != '-':
            edits += 'M'
        if p[i] == '-' and q[i] != '-':
            edits += 'I'
        if p[i] != '-' and q[i] == '-':
            edits += 'D'

    if len(p) == 0 and len(q) == 0:
        p_out = ''
        q_out = ''
        edits = ''
    else:
        p_out = p.replace('-','')
        q_out = q.replace('-','')
        edits = edits

    return p_out, q_out , edits



def align(p: str, q: str, edits: str):
    """Align two sequences from a sequence of edits.

    Args:
        p (str): The first sequence to align.
        q (str): The second sequence to align
        edits (str): The list of edits to apply, given as a string

    Returns:
        tuple[str, str]: The two rows in the pairwise alignment

    >>> align("ACCACAGTCATA", "ACAGAGTACAAA", "MDMMMMMMIMMMM")
    ('ACCACAGT-CATA', 'A-CAGAGTACAAA')

    """
    p_align = ''
    q_align = ''
    
    ins=0
    dels=0
    for i in range(len(edits)):
        
        if edits[i:] == 'M'*len(edits[i:]):
            p_align += p[i-ins:]
            q_align += q[i-dels:]
            break
        
        if edits[i] == 'M':
            p_align += p[i-ins]
            q_align += q[i-dels]

        if edits[i] == 'I':
            p_align += '-'
            q_align += q[i-dels]
            ins+=1

        if edits[i] == 'D':
            p_align += p[i-ins]
            q_align += '-'
            dels+=1

    return p_align, q_align


def local_align(p: str, x: str, i: int, edits: str):
    """Align two sequences from a sequence of edits.

    Args:
        p (str): The read string we have mapped against x
        x (str): The longer string we have mapped against
        i (int): The location where we have an approximative match
        edits (str): The list of edits to apply, given as a string

    Returns:
        tuple[str, str]: The two rows in the pairwise alignment

    >>> local_align("ACCACAGTCATA", "GTACAGAGTACAAA", 2, "MDMMMMMMIMMMM")
    ('ACCACAGT-CATA', 'A-CAGAGTACAAA')

    """
    edit_counts = edits.count('I') + edits.count('D')
    ref_subset = x[i:len(p)+edit_counts]

    return align(p,ref_subset,edits)


def edit_dist(p: str, x: str, i: int, edits: str):
    """Get the distance between p and the string that starts at x[i:]
    using the edits.

    Args:
        p (str): The read string we have mapped against x
        x (str): The longer string we have mapped against
        i (int): The location where we have an approximative match
        edits (str): The list of edits to apply, given as a string

    Returns:
        int: The distance from p to x[i:?] described by edits

    >>> edit_dist("accaaagta", "cgacaaatgtcca", 2, "MDMMIMMMMIIM")
    5
    """
    p = local_align(p,x,i,edits)[0]
    string = local_align(p,x,i,edits)[1]
    
    counts=0
    for i in range(len(p)):
        if p[i] != string[i]:
            counts+=1
    
    return counts

