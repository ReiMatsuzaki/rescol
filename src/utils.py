def flatten(xs):
    return reduce(lambda a, b: a+b, xs)


def uniq(xs):
    return list(set(xs))


def with_index(xs):
    return zip(xs, range(len(xs)))
