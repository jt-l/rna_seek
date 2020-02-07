import math
import fileinput


UNPAIRED_PENALTY = 0
UNPAIRED_PENALTY_IN_LOOP = 0.1
MULTIPLE_OPTIMAL_STRUCTURES = "NO"


def generate_structure(pairings, sequence):
    "Generate the min free energy structure in dot bracket notation"""

    # initalize all elements to .
    structure = ['.'] * (len(sequence))

    # iterate through the pairings
    # if there is a pairing replace the dots the indicies with brackets
    for pair in pairings:
        (i, j) = pair
        structure[i] = '('
        structure[j] = ')'

    return structure


def score_structure(structure, pairings, sequence):
    """ Calculate the min free energy of a given structure """

    structure_score = 0
    stack = []

    for char in structure:
        if char == "(":
            stack.append(')')

        elif char == ")":
            stack.pop()
            (i, j) = pairings.pop()
            structure_score += score(sequence, i, j)  # loop is finished

        else:
            if stack:
                structure_score += UNPAIRED_PENALTY_IN_LOOP

            structure_score += UNPAIRED_PENALTY

    # if the structure score is > 0, return 0
    return min(structure_score, 0)


def compute_traceback(d_matrix, sequence, pairings):
    """ Compute the traceback to find the structure"""
    global MULTIPLE_OPTIMAL_STRUCTURES

    traceback_stack = []
    traceback_path = []

    # first push (1, L) onto the stack
    traceback_stack.append((0, len(d_matrix)-1))
    traceback_path.append((0, len(d_matrix)-1))

    while traceback_stack:
        (i, j) = traceback_stack.pop()

        if i >= j:
            continue

        elif ((d_matrix[i+1][j] + UNPAIRED_PENALTY) == (d_matrix[i][j]
                                                        + UNPAIRED_PENALTY)):
            traceback_stack.append((i+1, j))
            traceback_path.append((i+1, j))

        elif ((d_matrix[i][j-1] + UNPAIRED_PENALTY) == (d_matrix[i][j]
                                                        + UNPAIRED_PENALTY)):
            traceback_stack.append((i, j-1))
            traceback_path.append((i, j-1))

        elif d_matrix[i+1][j-1] + score(sequence, i, j) == d_matrix[i][j]:
            MULTIPLE_OPTIMAL_STRUCTURES = "YES"
            pairings.append((i, j))

            traceback_stack.append((i+1, j-1))
            traceback_path.append((i+1, j-1))

            lb_l = (len(sequence) - 1) - j  # lower bound for left recurse
            ub_r = len(sequence) - 1  # upper bound for right recurse

            # recurse on left side
            if i - lb_l > 3:
                compute_traceback(
                    d_matrix[lb_l:i][lb_l:i],
                    sequence[lb_l:i], pairings)

            # recurse on right side
            if ub_r - j > 3:
                compute_traceback(
                    d_matrix[j:ub_r][j:ub_r], sequence[j:ub_r], pairings)

        else:
            for k in range(i+1, j-1, 1):
                if d_matrix[i][k] + d_matrix[k+1][j] == d_matrix[i][j]:
                    traceback_stack.append((k+1, j))
                    traceback_stack.append((i, k))
                    traceback_path.append((k+1, j))
                    traceback_path.append((i, k))
                    break


def generate_diagonal_indicies(d_matrix_size):
    """ helper function to generate indicies to be used for traversing
    d_matrix diagonally """

    diagonal_indicies = []

    # first generate the indicies for the diagonal below the half
    # this means i is always greater than 1
    j = 0
    for i in range(1, d_matrix_size):
        diagonal_indicies.append((i, j))
        j += 1

    # now generate the remainding indices
    iterations = 0
    while iterations <= 0:
        for i in range(d_matrix_size):
            diagonal = 0
            for j in range(d_matrix_size):
                if i+diagonal <= d_matrix_size - 1 and j <= d_matrix_size - 1:
                    diagonal_indicies.append((j, i+diagonal))
                    diagonal += 1

        iterations += 1

    return diagonal_indicies


def parse_fasta_file():
    """Parse FASTA formated input sequence the input file may or may not
    contain > SN"""

    lines = fileinput.input()

    # store lines[0] into a variable so that we do not access lines out of order
    line0 = lines[0]

    # The sequence is in fasta format if when we split on '>'
    # there are two items
    if len(line0.split('>')) == 2:
        return lines[1].strip()

    # otherwise we grab the first line
    return line0.strip()


def score(sequence, i, j):
    "Score a base pairing"

    base_0 = sequence[i]
    base_1 = sequence[j]

    # no loops of less than size 3 will be paired
    if j - i > 3:
        if ((base_0 == 'A' and base_1 == 'U') or
                (base_0 == 'U' and base_1 == 'A') or
                (base_0 == 'G' and base_1 == 'U') or
                (base_0 == 'U' and base_1 == 'G')):
            return -2

        elif ((base_0 == 'G' and base_1 == 'C') or
              (base_0 == 'C' and base_1 == 'G')):
            return -3

        else:
            return math.inf

    else:
        return math.inf


def recurrence(sequence, d_matrix, i, j):
    """ Main recurrence for the modified Nussinov Algorithm"""

    min_val = math.inf

    # first find the min for the last value in the recurrence
    for k in range(j):
        if (i < k and k > j):
            candidate_val = d_matrix[i, k] + d_matrix[k+1, j]
            if candidate_val < min_val:
                min_val = candidate_val

    final_min = min((d_matrix[i+1][j] + UNPAIRED_PENALTY),
                    (d_matrix[i][j-1] + UNPAIRED_PENALTY),
                    (d_matrix[i+1][j-1] + score(sequence, i, j)),
                    (min_val))

    return final_min


def main():
    """Main function to handle control and write contents outputs files"""
    sequence = parse_fasta_file()

    sequence_length = len(sequence)

    # initalize dynamic programming matrix
    d_matrix = [[0]*sequence_length for i in range(sequence_length)]

    # initalize the values of the matrix (redundant)
    for i in range(2, sequence_length):
        d_matrix[i][i-1] = 0

    for i in range(1, sequence_length):
        for j in range(1, sequence_length):
            if i == j:
                d_matrix[i][j] = 0

    # we need to traverse d_matrix diagonally, this function returns
    # the indicies we need to visit in order to correctly traverse diagonally
    diagonal_indicies = generate_diagonal_indicies(len(d_matrix))

    for indicies in diagonal_indicies:
        (i, j) = indicies[0], indicies[1]
        if i != j and i != j+1:
            d_matrix[i][j] = recurrence(
                sequence, d_matrix, i, j)

    # traceback
    pairings = []
    compute_traceback(d_matrix, sequence, pairings)

    # get the structure in dot bracket notation
    structure = generate_structure(pairings, sequence)

    # get structure score
    structure_score = score_structure(structure, pairings, sequence)

    # ensure structures with scores of 0 are displayed correctly
    if structure_score == 0:
        structure = ['.' * len(sequence)]

    # write output file
    sequence_file = open("5.o1", "w+")
    sequence_file.write(str(sequence) + "\n")
    sequence_file.write(str(structure_score) + "\n")
    sequence_file.write(str("".join(structure)))

    # write output file for multiple structures
    multiple_optimal_structures_file = open("5.o2", "w+")
    multiple_optimal_structures_file.write(str(MULTIPLE_OPTIMAL_STRUCTURES))


main()
