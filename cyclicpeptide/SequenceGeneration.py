import csv
import itertools


def generate_similar_sequences(sequence, replacement_rules_file=None, Replacement_ratio=None):
    """
    Generate similar sequences based on a given sequence and optional replacement rules.

    :param sequence: The input sequence for which similar sequences are to be generated.
    :type sequence: str
    :param replacement_rules_file: (Optional) The file path containing the replacement rules. If not provided, default replacement rules will be used.
    :type replacement_rules_file: str or None
    :param Replacement_ratio: (Optional) The ratio of amino acids to be replaced. If not provided, it will be set to one-third of the total amino acid count.
    :type Replacement_ratio: float or None
    :return: A list of possible similar sequences generated based on the input sequence, replacement rules (either default or from file), and the replacement ratio (either default or specified).
    :rtype: list[str]

    This function generates similar sequences to the given input sequence. 
    """

    # Use the default rules if no external file is specified
    if replacement_rules_file is None:
        replacement_rules = {
            'Glu': ['Asp'],  # Acidic amino acids
            'Asp': ['Glu'],
            'Lys': ['Arg', 'His'],  # Basic amino acids
            'Arg': ['Lys', 'His'],
            'His': ['Lys', 'Arg'],
            'Leu': ['Ile', 'Val', 'Ala'],  # Hydrophobic amino acids
            'Ile': ['Leu', 'Val', 'Ala'],
            'Val': ['Leu', 'Ile', 'Ala'],
            'Ala': ['Leu', 'Ile', 'Val'],
            'Phe': ['Tyr', 'Trp'],  # Aromatic amino acids
            'Tyr': ['Phe', 'Trp'],
            'Trp': ['Phe', 'Tyr'],
            'Ser': ['Thr'],  # Polar amino acids
            'Thr': ['Ser'],
            'Gly': ['Ala'],  # Small-volume amino acids
            'Cys': ['Ser'],  # Special functional group (polar replacement)
            'Met': ['Leu', 'Ile'],  # Hydrophobic/non-polar amino acids
            'Asn': ['Gln'],  # Polar amino acids with an amide group
            'Gln': ['Asn'],
        }
    else:
        replacement_rules = read_replacement_rules_from_file(replacement_rules_file)
    total_aa_count = len(sequence.split('--'))
    # If no maximum number of replacements is given, use one-third of the total residue count
    if Replacement_ratio is None:
        max_replacements = total_aa_count // 3
    else:
        max_replacements = total_aa_count * Replacement_ratio

    # Match amino acids and their modifications with regex (modification-specific logic is not handled here yet)
    fragments = sequence.split('--')
    all_combinations = []

    # Iterate over each fragment to build its replacement options
    for fragment in fragments:
        aa = fragment
        options = []
        if aa in replacement_rules and replacement_rules[aa]:
            # Options include keeping the original or replacing it
            options.append(aa)
            for replacement in replacement_rules[aa]:
                options.append(replacement)
        else:
            options.append(aa)

        all_combinations.append(options)

    # Generate all possible replacement combinations
    possible_sequences = []
    for combo in itertools.product(*all_combinations):
        replacements_count = sum(1 for i in range(len(fragments)) if combo[i]!= fragments[i])
        if replacements_count >= 0 and replacements_count <= max_replacements:
            new_sequence = '--'.join(combo)
            possible_sequences.append(new_sequence)

    return possible_sequences


def read_replacement_rules_from_file(file_path):
    """
    Read replacement rules from a specified file.

    :param file_path: The path to the file containing the replacement rules.
    :type file_path: str
    :return: A dictionary where the keys are amino acids and the values are lists of possible replacement amino acids as read from the file.
    :rtype: dict

    This function reads replacement rules from a given file. 
    """

    replacement_rules = {}
    with open(file_path, 'r') as file:
        for line in file.readlines():
            line = line.strip()
            if line:
                key, values_str = line.split(':')
                values = [value.strip() for value in values_str.strip('[]').spanplit(',')]
                replacement_rules[key.strip()] = values
    return replacement_rules


def save_sequences_to_csv(sequences, filename):
    """
    Save a list of sequences to a CSV file.
    
    :param sequences: A list of sequences to be saved to the CSV file.
    :type sequences: list[str]
    :param filename: The name of the CSV file to which the sequences will be saved.
    :type filename: str
    :return: None

    This function takes a list of sequences and saves them to a CSV file with the specified filename. 

    """

    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Sequence'])  # Write the header row
        for sequence in sequences:
            writer.writerow([sequence])