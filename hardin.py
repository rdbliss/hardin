"""
Compute generating functions and closed forms for H_r(n, k), the number of n x
k Hardinian arrays with parameter r.

The sequence H_r(n, k) is a polynomial in n of degree r for sufficiently large
n > k. (And vice versa, as H_r(n, k) = H_r(k, n).) To compute the polynomial
that H_r(n, k) eventually equals, try this:

    >>> n, r, k = symbol("n"), 2, 3
    >>> rectPoly(n, r, k)

To compute the generating function of H_r(n, k) with respect to n, try this:

    >>> x, r, k = symbol("x"), 2, 3
    >>> rectGF(x, r, k)

To be precise, rectGF computes the generating function

    sum(x^n H_r(n, k), n >= k).
"""
import itertools
from collections import deque
from sympy import SparseMatrix, var, eye
from sympy.series.formal import rational_algorithm
from sympy.polys.matrices import DomainMatrix

def rectGF(x, r, k):
    """
    Return the generating function for H_r(n, k) with respect to n for n >= k.
    (The terms with exponents less than k are set to 0.) The sequence is a
    polynomial of degree r, plus a finite number of correction terms.
    """
    _, invMap, edges, dim = squareTM(r, k)
    m = SparseMatrix(dim, dim, edges)
    m = DomainMatrix.from_Matrix(m)
    m = m**k

    # figure out what the possible final "square part" rows are by looking at
    # the last row of the matrix power
    start_states = []
    for key, count in m[-1, :].to_dok().items():
        index = key[-1]
        state = invMap[index]
        row_index, row = state
        # type conversion because sympy stores entries as mpz items
        start_states.append((row, int(count)))

    rowMap, invMap, edges, dim = rectTM(r, k)
    m = SparseMatrix(dim, dim, edges)

    stopCols = [k for k in range(dim) if is_final_rect(invMap[k])]

    m = SparseMatrix(eye(m.shape[0])) - x * m
    m = DomainMatrix.from_Matrix(m).to_field().inv()

    gf = 0
    for row, count in start_states:
        index = rowMap[row]
        gf += sum(count * m[index, stop].to_sympy() for stop in stopCols)

    return (x**k * gf).factor()

def rectPoly(n, r, k):
    """
    Return the polynomial that H_r(n, k) eventually equals.
    """
    x = var("x", dummy=True)
    f = rectGF(x, r, k)
    formula = rational_algorithm(f, x, n)[0].simplify()

    return formula

def valid_follow(i, r, row, new_row):
    """
    Check if `new_row` could follow `row` as the (i + 1)th row of a Hardinian
    array with parameter r.
    """

    x = new_row[0] + i + 1 - r
    w = row[0] + i - r

    if x - w not in {0, 1}:
        return False

    if new_row[0] < 0:
        return False

    for j in range(len(row) - 1):
        # does new_row satisfy the king-distance rule?

        if i > j:
            maxij = i
            maxij1 = i
            maxi1j1 = i + 1
            maxi1j = i + 1
        elif i < j:
            maxij = j
            maxij1 = j + 1
            maxi1j1 = j + 1
            maxi1j = j
        else:
            maxij = j
            maxij1 = j + 1
            maxi1j1 = j + 1
            maxi1j = i + 1

        x = new_row[j] + maxi1j - r
        y = new_row[j + 1] + maxi1j1 - r
        w = row[j] + maxij - r
        z = row[j + 1] + maxij1 - r

        if new_row[j + 1] < 0:
            return False

        if y - x not in {0, 1}:
            return False
        if y - w not in {0, 1}:
            return False
        if y - z not in {0, 1}:
            return False

    return True

def gen_follows_square(r, point):
    i, row = point
    k = len(row)
    new_row = [0] * k
    row_delta = [row[k + 1] - row[k] for k in range(len(row) - 1)]

    unknown = set(range(k))
    for j, delta in enumerate(row_delta):
        # if j < i, then a step to (i, j + 1) keeps the king-distance the same,
        # so an increase of 1 implies the "hardinian" value increased. this
        # implies that row[j + 1] and new_row[j + 1] have the same hardinian
        # value (but not necessarily the same transformed value). similar
        # argument for j >= i.
        if (j < i and delta == 1) or (j >= i and delta == 0):
            hard_val = row[j + 1] + max(i, j + 1) - r
            new_row[j + 1] = hard_val - max(i + 1, j + 1) + r
            unknown -= {j + 1}

    # we don't know what new_row[j] is for j in unknown.
    # we could probably do better than the below code, but it gets the job
    # done.
    for delta in itertools.product({0, 1}, repeat=len(unknown)):
        for j, delt in zip(unknown, delta):
            if j <= i:
                new_row[j] = row[j] - delt
            else:
                new_row[j] = row[j] + delt

        if valid_follow(i, r, row, new_row):
            yield (i + 1, tuple(new_row))

def gen_inits_square(r, k):
    """
    Return all valid initial rows of a reduced Hardinian array with parameter r
    and k columns.
    """
    # every step to the right increases the king distance by 1. using
    #   T = M - KD + r,
    # we see that T is changing by 0 or -1 each step to the right. it can
    # decrease at most r times since its initial value is r.

    for decreases in range(r + 1):
        choices = itertools.combinations(range(1, k + 1), decreases)
        for choice in choices:
            row = [r - sum(True for j in choice if k >= j) for k in range(k)]
            yield (0, tuple(row))

def squareTM(r, k):
    """
    Return the transition matrix for the "square part" of a reduced Hardinian
    array with parameter r and k rows.
    """
    def is_final(point):
        i, row = point
        return i == k - 1

    def is_initial(point):
        i, _ = point
        return i == 0

    return TM(r, k, gen_inits_square, gen_follows_square, is_final, is_initial=is_initial)

def TM(r, k, gen_inits, gen_follows, is_final, is_initial=None):
    """
    Compute a transition matrix for Hardinian arrays with parameter r and k
    rows, with parameters allowing for square and non-square parts to be
    generated.

    The return is

        (stateToCol, colToState, edges, dimension),

    where stateToCol and colToState are dictionaries mapping states to column
    indicies of the matrix, and vice versa. `edges` is a sparse dictionary
    representation of the edges, and `dimension` is the size of the matrix.

    gen_inits(r, k): function that produces initial states.

    gen_follows(r, row): function that produces valid rows which could follow
                         the current row.

    is_final(row): checks if the given row is final.

    is_initial: optional function that checks whether a given state is initial.
                if this is given, then a special "start" state will be added
                which leads only to initial states.
    """
    to_visit = deque([tuple(row) for row in gen_inits(r, k)])
    rows_seen = set()
    edges = dict()

    while to_visit:
        row = to_visit.popleft()
        if row in rows_seen:
            continue

        rows_seen |= {row}

        if is_initial and is_initial(row):
            edges[(None, row)] = 1

        if is_final(row):
            # mark this as a potential "final row" by giving it a loop.
            edges[(row, row)] = 1
            continue

        follows = gen_follows(r, row)
        for follow in follows:
            if follow not in rows_seen:
                to_visit.append(follow)

            edges[(row, follow)] = 1

    # this is the most efficient map to turn rows to integers.
    intMap = {point: k for k, point in enumerate(rows_seen)}
    invMap = {k: point for k, point in enumerate(rows_seen)}
    newEdges = dict()

    needed = len(rows_seen)

    for edge in edges.keys():
        start, stop = edge

        if start is None:
            # this is a special "start" state. send it beyond the theoretical maximum.
            newEdges[(needed, intMap[stop])] = 1
        else:
            newEdges[(intMap[start], intMap[stop])] = 1

    if is_initial:
        needed += 1

    return intMap, invMap, newEdges, needed

def is_final_rect(row):
    return row[-1] == 0

def rectTM(r, k, *args, **kwargs):
    """
    Return the transition matrix for the "rectangular part" of a reduced
    Hardinian array with parameter r and k rows.
    """
    return TM(r, k, gen_inits_rect, gen_follows_rect, is_final_rect, *args, **kwargs)

def gen_inits_rect(r, k):
    """
    Return all valid kth rows of a reduced Hardinian array with parameter r
    and k columns.
    """
    # the king distance is the same across the entire row. since
    #   T = M - KD + r,
    # T is changing by 0 or 1 each step to the right.

    for init in range(r + 1):
        for increases in range(r + 1 - init):
            choices = itertools.combinations(range(1, k + 1), increases)
            for choice in choices:
                row = [init + sum(True for j in choice if k >= j) for k in range(k)]
                yield row

def gen_follows_rect(r, row):
    """
    Return all rows which could follow `row` in the rectangular part of a
    Hardinian array with parameter r and len(row) columns. "Rectangular part"
    means j >= len(row), so the king distance is constant along every row and
    increases by 1 for every step down.
    """
    k = len(row)
    new_row = [0] * k
    row_delta = [row[k + 1] - row[k] for k in range(len(row) - 1)]

    unknown = set(range(k))
    for j, delta in enumerate(row_delta):
        if delta == 1:
            # T = M - KD + r
            # new_row[j + 1] and row[j] have the same Hardinian value, and king
            # distance went up up by 1 from the step down.
            new_row[j + 1] = row[j + 1] - 1
            unknown -= {j + 1}

    for delta in itertools.product({0, 1}, repeat=len(unknown)):
        for j, delt in zip(unknown, delta):
            new_row[j] = row[j] - delt

        if valid_follow_rect(r, row, new_row):
            yield tuple(new_row)

def valid_follow_rect(r, row, new_row):
    """
    Check if `new_row` could follow `row` in the rectangular part of a
    Hardinian array with parameter r.
    """

    # T = M - max(i, j) + r
    if row[0] - new_row[0] not in {0, 1}:
        return False

    if new_row[0] < 0:
        return False

    for j in range(len(row) - 1):
        if new_row[j + 1] - new_row[j] not in {0, 1}:
            return False

        if row[j + 1] - new_row[j + 1] not in {0, 1}:
            return False

        if row[j] - new_row[j + 1] not in {0, 1}:
            return False

    return True

def printrectable(K):
    for k in range(1, K):
        _, _, edges, dim = rectTM(2, k)
        m = SparseMatrix(dim, dim, edges)
        print(f"{k} & {m.shape[0]} & {m.nnz() / m.shape[0]**2:.3f} \\\\")

def printsquaretable(K):
    for k in range(1, K):
        _, _, edges, dim = squareTM(2, k)
        m = SparseMatrix(dim, dim, edges)
        print(f"{k} & {m.shape[0]} & {m.nnz() / m.shape[0]**2:.3f} \\\\")

def demo():
    n, x = var("n x")
    for r in [2, 3]:
        print(f"r = {r}")
        for k in range(1, 5):
            print(f"k = {k}")
            print("generating function:", rectGF(x, r, k))
            print("eventual polynomial:", rectPoly(n, r, k))

        print()
