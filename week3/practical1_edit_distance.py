import timeit


def edit_distance(x, y):
    # Create distance matrix
    D = []
    for i in range(len(x) + 1):
        D.append([0] * (len(y) + 1))
    # Initialize first row and column of matrix
    for i in range(len(x) + 1):
        D[i][0] = i
    for i in range(len(y) + 1):
        D[0][i] = i
    # Fill in the rest of the matrix
    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            distHor = D[i][j - 1] + 1
            distVer = D[i - 1][j] + 1
            if x[i - 1] == y[j - 1]:
                distDiag = D[i - 1][j - 1]
            else:
                distDiag = D[i - 1][j - 1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    return D[-1][-1]


def edit_dist_recursive(x, y):
    # This implementation is very slow
    if len(x) == 0:
        return len(y)
    elif len(y) == 0:
        return len(x)
    else:
        distHor = edit_dist_recursive(x[:-1], y) + 1
        distVer = edit_dist_recursive(x, y[:-1]) + 1
        if x[-1] == y[-1]:
            distDiag = edit_dist_recursive(x[:-1], y[:-1])
        else:
            distDiag = edit_dist_recursive(x[:-1], y[:-1]) + 1
        return min(distHor, distVer, distDiag)


_cache = {}
def cache(func):
    def wrapper(x, y):
        if (x, y) not in _cache:
            _cache[(x, y)] = func(x, y)
            _cache[(y, x)] = _cache[(x, y)]
        return _cache[(x, y)]

    return wrapper


@cache
def edit_distance_cache(x, y):
    # This implementation is very slow
    if len(x) == 0:
        return len(y)
    elif len(y) == 0:
        return len(x)
    else:
        distHor = edit_dist_recursive(x[:-1], y) + 1
        distVer = edit_dist_recursive(x, y[:-1]) + 1
        if x[-1] == y[-1]:
            distDiag = edit_dist_recursive(x[:-1], y[:-1])
        else:
            distDiag = edit_dist_recursive(x[:-1], y[:-1]) + 1
        return min(distHor, distVer, distDiag)


if __name__ == "__main__":
    setup = '''
    x = 'shake spea'
    y = 'Shakespear'
    '''
    code_to_execute = "edit_dist_recursive(x, y)"

    print("Average time needed to calculate edit distance by recursive function: ",
          timeit.timeit(stmt=code_to_execute, setup=setup.replace('    ', ''), number=10,
                        globals=globals()) / 10)

    setup = '''
    x = 'shake spea'
    y = 'Shakespear'
        '''
    code_to_execute = "edit_distance(x, y)"

    print("Average time needed to calculate edit distance by function without recursion: ",
          timeit.timeit(stmt=code_to_execute, setup=setup.replace('    ', ''), number=1000,
                        globals=globals()) / 1000)

    setup = '''
    x = 'shake spea'
    y = 'Shakespear'
    '''

    code_to_execute = "edit_distance_cache(x, y)"

    print("Average time needed to calculate edit distance by function with cache wrapper: ",
          timeit.timeit(stmt=code_to_execute, setup=setup.replace('    ', ''), number=10,
                        globals=globals()) / 10)
