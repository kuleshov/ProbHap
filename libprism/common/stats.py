def N50(l):
    l = list(l)
    total_length = sum(l)
    half_length = total_length / 2
    length = 0

    l.sort(reverse=True)

    x=0
    while length < half_length:
        x = l.pop()
        length += x
        # print "%d/%d" % (length, total_length)

    return x

def div(a, b):
    return a/b if b else 0