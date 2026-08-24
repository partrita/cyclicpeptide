import re
from .setting import *

# Open the text file and read its contents
file_path = AminoAcids_path  # Replace with your file path
AminoAcids = []
with open(file_path, 'r', encoding='utf-8') as file:
    # Read the file line by line and append each entry to the list
    for line in file:
        # Strip the trailing newline and append the entry to the list
        AminoAcids.append(eval(line.strip()))

################## Sequence Format Conversion ##################
def read_sequence(sequence):

    """
    **read_sequence** reads sequences in a variety of formats, converting the item into nodes(lists of amino acids) and edges(link information).
    you can also use exact function to transform, such as :py:func:`read_iupac_condensed`, :py:func:`read_graph_representation`, and :py:func:`read_one_letter_sequence`

    :param sequence: cyclic peptide sequence.
    :type sequence: *Format* [Graph presentation, IUPAC condensed, Amino acid chain, One letter code]
    :return: seq_format, nodes, edges
    :rtype: string, List[str], List[(int, int)]

    Example::

        seq_format, nodes, edges = SequenceTransformer.read_sequence('aThr,Tyr,dhAbu,bOH-Gln,Gly,Gln,His,Dab,C13:2(t4.t6)-OH(2.3),Lyx,dhAbu @1,5 @6,10 @0,8')
        print(seq_format, nodes, edges)

    Results::
    
        Graph presentation\n
        ['aThr', 'Tyr', 'dhAbu', 'bOH-Gln', 'Gly', 'Gln', 'His', 'Dab', 'C13:2(t4.t6)-OH(2.3)', 'Lyx', 'dhAbu']\n
        [(1, 5), (6, 10), (0, 8)]

    """

    seq_format = ''
    if '@' in sequence:
        seq_format = 'Graph presentation'
        nodes, edges = read_graph_representation(sequence)
    elif ',' in sequence:
        no_parentheses = re.sub(r'\([^)]*\)', '', sequence)  # Strip all parenthesized parts
        if ',' in no_parentheses:
            seq_format = 'Graph presentation'
            nodes, edges = read_graph_representation(sequence)
        else:
            seq_format = 'IUPAC condensed'
            nodes, edges = read_iupac_condensed(sequence, sep='-')
    elif '--' in sequence:
        seq_format = 'Amino acid chain'
        nodes, edges = read_iupac_condensed(sequence, sep='--')
    elif '-' in sequence:
        try:
            seq_format = 'IUPAC condensed'
            nodes, edges = read_iupac_condensed(sequence, sep='-')
        except:
            nodes, edges = '', ''
    else:
        seq_format = 'One letter peptide'
        nodes, edges = read_one_letter_sequence(sequence)
    return seq_format, nodes, edges


def create_sequence(nodes, edges):
    """
    **create_sequence** can generate processed nodes(amino acid lists) and edges(link information) into sequences in many different formats.
    you can also use exact function to create, such as :py:func:`create_iupac_condensed`, :py:func:`create_graph_presentation`, :py:func:`create_amino_acid_chain`, and :py:func:`create_one_letter_peptide`

    :param nodes: lists of amino acids.
    :param edges: link information.
    :type edges: List[(int, int)]
    :return: {'iupac_condensed': iupac_condensed, 'amino_acid_chain': amino_acid_chain,
            'graph_presentation': graph_presentation, 'one_letter_peptide': one_letter_peptide}

    Example::

        seq_list = SequenceTransformer.create_sequence(['aThr', 'Tyr', 'dhAbu', 'bOH-Gln', 'Gly', 'Gln', 'His', 'Dab', 'C13:2(t4.t6)-OH(2.3)', 'Lyx', 'dhAbu'], [(1, 5), (6, 10), (0, 8)])
        print(seq_list)
    
    Results::

        {\n
        'iupac_condensed': 'aThr(3)-Tyr(1)-dhAbu-bOH-Gln-Gly-Gln(1)-His(2)-Dab-C13:2(t4.t6)-OH(2.3)(3)-Lyx-dhAbu(2)',\n
        'amino_acid_chain': 'aThr(3)--Tyr(1)--dhAbu--bOH-Gln--Gly--Gln(1)--His(2)--Dab--C13:2(t4.t6)-OH(2.3)(3)--Lyx--dhAbu(2)',\n
        'graph_presentation': 'aThr,Tyr,dhAbu,bOH-Gln,Gly,Gln,His,Dab,C13:2(t4.t6)-OH(2.3),Lyx,dhAbu @1,5 @6,10 @0,8',\n
        'one_letter_peptide': None\n
        }

    """
    try:
        iupac_condensed = create_iupac_condensed(nodes, edges)
    except:
        iupac_condensed = None
    try:
        amino_acid_chain = create_amino_acid_chain(nodes, edges)
    except:
        amino_acid_chain = None
    try:
        graph_presentation = create_graph_presentation(nodes, edges)
    except:
        graph_presentation = None
    try:
        one_letter_peptide = create_one_letter_peptide(nodes)
    except:
        one_letter_peptide = None
    return {'iupac_condensed': iupac_condensed, 'amino_acid_chain': amino_acid_chain,
            'graph_presentation': graph_presentation, 'one_letter_peptide': one_letter_peptide}


def read_iupac_condensed(sequence, sep='-'):
    """

    :param sequence: cyclic peptide sequence.
    :type sequence: IUPAC condensed
    :param sep: Optional "sep".
    :type sep: Amino acid link in sequences
    :return: amino acids, edges
    :rtype: List[str], List[(int, int)]


    """

    # 0. Check for salts or multiple chains
    if '.' in sequence:
        sequence = sequence.split('.')[0].strip()
    # 1. Check for cyclo[] cyclization info
    iscyclo = True if 'cyclo' in sequence else False
    sequence = sequence.replace('cyclo', '').strip('[]')  # Remove 'cyclo' and the square brackets
    sequence = re.sub(r'\([^)]*\)', replace_hyphen, sequence)  # Replace '-' inside parentheses with '~~' (parentheses hold modifications)
    items = sequence.split(sep)
    # 2. Check terminal (head/tail) modification info
    header = ''
    tail = ''
    if items[0] in ['H', 'NH2', 'Unk']:
        header = items[0]
        items = items[1:]
    if items[-1] in ['H', 'NH2', 'Unk']:
        tail = items[-1]
        items = items[:-1]
    # 3. Extract amino acids; DL-/D-/L- prefixes and terminal modifications are merged into the amino acid names
    amino_acids = []
    edge_marks = []
    pred = ''
    max_edge_mark = 0
    for item in items:
        if item in ['DL', 'D', 'L']:
            pred = item + '-'
            continue
        amino_acid = (pred + item).replace('~~', '-').strip()
        edge_mark = []
        while True:
            edge_mark_match = re.search(r'\(\d+\)$', amino_acid)
            if edge_mark_match:
                amino_acid = amino_acid.replace(edge_mark_match[0], '').strip()
                em = int(edge_mark_match[0][1:-1])
                max_edge_mark = max([max_edge_mark, em])
                edge_mark.append(em)
            else:
                break
        amino_acids.append(amino_acid)
        edge_marks.append(edge_mark)
        pred = ''
    if header:
        amino_acids[0] = header + '-' + amino_acids[0]
    if tail:
        amino_acids[-1] = amino_acids[-1] + '-' + tail
    # 4. Handle an empty trailing item that carries edge labels
    if amino_acids[-1] == '' and edge_marks:
        amino_acids = amino_acids[:-1]
        edge_mark = edge_marks[-1][:]
        edge_marks = edge_marks[:-1]
        edge_marks[-1].extend(edge_mark)
    # 5. Handle amino acid N(1) annotations
    for i, aa in enumerate(amino_acids):
        match = re.search(r'N\((\d)\)', aa)
        if match:
            amino_acids[i] = re.sub(r'N\((\d)\)', 'N-', aa)
            edge_marks[i].append(int(match.group(1)))

    # print(amino_acids)
    # print(edge_marks)
    # 6. Process the edge info
    edges = []
    for i in range(1, max_edge_mark + 1):
        pair = []
        for idx in range(len(edge_marks)):
            if i in edge_marks[idx]:
                pair.append(idx)
        edges.append((pair[0], pair[1]))
    if iscyclo:
        edges.append((0, len(amino_acids) - 1))
    return amino_acids, edges


def read_graph_representation(sequence):
    """

    :param sequence: cyclic peptide sequence.
    :type sequence: graph representation
    :return: amino acids, edges
    :rtype: List[str], List[(int, int)]

    """

    amin_acids = [i.strip() for i in sequence.split(' @')[0].strip().split(',')]
    edges = [tuple([int(j) for j in i.strip().split(',')]) for i in sequence.split('@')[1:]]
    return amin_acids, edges


def read_one_letter_sequence(sequence):
    """

    :param sequence: cyclic peptide sequence.
    :type sequence: one letter code
    :return: amino acids, []
    :rtype: List[str], []

    """
    tail = ''
    if '(NH2)' == sequence.upper()[-5:]:
        sequence = sequence[:-5]
        tail = '(NH2)'
    amino_acid_refs = {i: j for i, j, _, _ in AminoAcids}
    amino_acid_refs['X'] = 'UNK'
    aas = [amino_acid_refs[i] for i in sequence]
    aas[-1] = aas[-1] + tail
    return aas, []


def create_iupac_condensed(nodes, edges):
    """    
    
    :param nodes: lists of amino acids.
    :param edges: link information.
    :type edges: List[(int, int)]
    :return: sequence
    :rtype: IUPAC condensed
    
    """
    edge_marker = 1
    cyclo = False
    amino_acids = nodes[:]
    for i, j in edges:
        if min(i, j) == 0 and max(i, j) == len(amino_acids) - 1:
            cyclo = True
            continue
        amino_acids[i] += f'({edge_marker})'
        amino_acids[j] += f'({edge_marker})'
        edge_marker += 1
    sequence = '-'.join(amino_acids)
    sequence = f'cyclo[{sequence}]' if cyclo else sequence
    return sequence


def create_amino_acid_chain(nodes, edges):
    """    
    
    :param nodes: lists of amino acids.
    :param edges: link information.
    :type edges: List[(int, int)]
    :return: sequence
    :rtype: Amino acid chain
    
    """

    edge_marker = 1
    amino_acids = nodes[:]
    for i, j in edges:
        amino_acids[i] += f'({edge_marker})'
        amino_acids[j] += f'({edge_marker})'
        edge_marker += 1
    sequence = '--'.join(amino_acids)
    return sequence


def create_graph_presentation(nodes, edges):
    """    
    
    :param nodes: lists of amino acids.
    :param edges: link information.
    :type edges: List[(int, int)]
    :return: sequence
    :rtype: Graph presentation
    
    """
    sequence = ','.join(nodes)
    edges_pre = '' if len(edges) == 0 else ' '.join(['@' + ','.join([str(k) for k in pair]) for pair in edges])
    return (sequence + ' ' + edges_pre).strip()


def create_one_letter_peptide(nodes):  # Side-chain interactions are ignored
    """    
    
    :param nodes: lists of amino acids.
    :param edges: link information.
    :type edges: List[(int, int)]
    :return: sequence
    :rtype: one letter peptide
    
    """

    amino_acid_refs = {i: j for j, i, _, _ in AminoAcids}
    amino_acid_refs['UNK'] = 'X'
    tail = ''
    if '(NH2)' == nodes[-1].upper()[-5:]:
        nodes[-1] = nodes[-1][:-5]
        tail = '(NH2)'
    is_essential_aa = [i.upper().strip() in amino_acid_refs.keys() for i in nodes]
    if False in is_essential_aa:
        return None
    else:
        return ''.join([amino_acid_refs[i.upper().strip()] for i in nodes]) + tail


def replace_hyphen(match):
    """
    This function is designed to replace all occurrences of the hyphen character ('-') within a matched string with the tilde character ('~').

    :param match: The match object obtained from a regex matching operation. It represents the portion of the string that matched a particular pattern.
    :type match: re.Match
    :return: The modified version of the matched string where all hyphens have been replaced with tildes.
    :rtype: str

    When called with a valid match object, it accesses the entire matched string using `match.group(0)` and then replaces all hyphens within that string with tildes. This can be useful in scenarios where specific text formatting or substitution within a particular part of a larger text (as identified by the regex match) is required.
    """
    # Replace all hyphens (-) inside parentheses with tildes (~)
    return match.group(0).replace('-', '~~')