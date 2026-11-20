"""
qr_engine.py - a from-scratch QR code encoder AND damaged-code decoder.

Runs in the browser via Pyodide (and also under normal CPython for testing).
No libraries: everything here, including the Galois-field arithmetic and the
Reed-Solomon error corrector, is plain Python.

Pipeline (encode):
    text -> bytes -> bit stream -> data codewords
         -> split into blocks -> Reed-Solomon EC codewords (polynomial division in GF(2^8))
         -> interleave -> zig-zag placement -> XOR mask (8 candidates, lowest penalty wins)
         -> format info (BCH-protected) -> grid

Pipeline (decode, what a phone does after it has a grid of black/white squares):
    1. locate the 3 finder patterns
    2. read the 15-bit format info (which EC level + which mask)
    3. un-mask and read codewords in zig-zag order, de-interleave into blocks
    4. Reed-Solomon: syndromes -> Berlekamp-Massey -> Chien search -> Forney
       => find WHERE the errors are and WHAT the right values were
    5. parse the bit stream back into text
"""

import json

# --------------------------------------------------------------------------
# 1. Arithmetic in GF(2^8)  (the "byte field": 256 elements, + is XOR)
# --------------------------------------------------------------------------
PRIM = 0x11D  # x^8 + x^4 + x^3 + x^2 + 1, the polynomial QR codes use
EXP = [0] * 512
LOG = [0] * 256
_x = 1
for _i in range(255):
    EXP[_i] = _x
    LOG[_x] = _i
    _x <<= 1
    if _x & 0x100:
        _x ^= PRIM
for _i in range(255, 512):
    EXP[_i] = EXP[_i - 255]


def gmul(a, b):
    if a == 0 or b == 0:
        return 0
    return EXP[LOG[a] + LOG[b]]


def gdiv(a, b):
    if b == 0:
        raise ZeroDivisionError()
    if a == 0:
        return 0
    return EXP[(LOG[a] - LOG[b]) % 255]


def gpow(a, n):
    return EXP[(LOG[a] * n) % 255]


def ginv(a):
    return EXP[(255 - LOG[a]) % 255]


# Polynomials are lists of coefficients, highest degree first.
def p_scale(p, s):
    return [gmul(c, s) for c in p]


def p_add(p, q):
    r = [0] * max(len(p), len(q))
    for i, c in enumerate(p):
        r[i + len(r) - len(p)] = c
    for i, c in enumerate(q):
        r[i + len(r) - len(q)] ^= c
    return r


def p_mul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for j, b in enumerate(q):
        for i, a in enumerate(p):
            r[i + j] ^= gmul(a, b)
    return r


def p_eval(p, x):
    y = p[0]
    for c in p[1:]:
        y = gmul(y, x) ^ c
    return y


def p_div(dividend, divisor):
    """Synthetic division by a monic divisor. Returns (quotient, remainder)."""
    out = list(dividend)
    for i in range(len(dividend) - (len(divisor) - 1)):
        c = out[i]
        if c:
            for j in range(1, len(divisor)):
                if divisor[j]:
                    out[i + j] ^= gmul(divisor[j], c)
    sep = len(out) - (len(divisor) - 1)
    return out[:sep], out[sep:]


# --------------------------------------------------------------------------
# 2. Reed-Solomon
# --------------------------------------------------------------------------
_GEN = {}


def rs_generator(n):
    """g(x) = (x - a^0)(x - a^1)...(x - a^(n-1))"""
    if n not in _GEN:
        g = [1]
        for i in range(n):
            g = p_mul(g, [1, EXP[i]])
        _GEN[n] = g
    return _GEN[n]


def rs_encode(data, n):
    """EC codewords = remainder of  data(x) * x^n  divided by  g(x)."""
    _, rem = p_div(list(data) + [0] * n, rs_generator(n))
    return rem


class RSError(Exception):
    pass


def rs_syndromes(msg, n):
    """S_i = r(a^i). All zero  <=>  msg is a valid codeword (no detectable errors)."""
    return [p_eval(msg, EXP[i]) for i in range(n)]


def rs_decode(msg, n):
    """
    Correct up to n//2 byte errors at UNKNOWN positions.
    Returns (corrected, error_positions, syndromes). Raises RSError if impossible.
    """
    msg = list(msg)
    synd_list = rs_syndromes(msg, n)
    if max(synd_list) == 0:
        return msg, [], synd_list
    synd = [0] + synd_list  # padded copy, convenient indexing

    # Berlekamp-Massey: find the error-locator polynomial
    err_loc = [1]
    old_loc = [1]
    for i in range(n):
        k = i + 1
        delta = synd[k]
        for j in range(1, len(err_loc)):
            delta ^= gmul(err_loc[-(j + 1)], synd[k - j])
        old_loc = old_loc + [0]
        if delta != 0:
            if len(old_loc) > len(err_loc):
                new_loc = p_scale(old_loc, delta)
                old_loc = p_scale(err_loc, ginv(delta))
                err_loc = new_loc
            err_loc = p_add(err_loc, p_scale(old_loc, delta))
    while err_loc and err_loc[0] == 0:
        err_loc.pop(0)
    errs = len(err_loc) - 1
    if errs * 2 > n:
        raise RSError("too many errors")

    # Chien search: the roots of the locator tell us the error POSITIONS
    L = len(msg)
    rev = err_loc[::-1]
    err_pos = []
    for i in range(L):
        if p_eval(rev, gpow(2, i)) == 0:
            err_pos.append(L - 1 - i)
    if len(err_pos) != errs:
        raise RSError("locator roots do not match")

    # Forney: compute the error MAGNITUDES (what to XOR to fix each byte)
    coef_pos = [L - 1 - p for p in err_pos]
    e_loc = [1]
    for i in coef_pos:
        e_loc = p_mul(e_loc, p_add([1], [gpow(2, i), 0]))
    _, err_eval = p_div(p_mul(synd[::-1], e_loc), [1] + [0] * len(e_loc))
    X = [gpow(2, -(255 - cp)) for cp in coef_pos]
    E = [0] * L
    for i, Xi in enumerate(X):
        Xi_inv = ginv(Xi)
        prime = 1
        for j, Xj in enumerate(X):
            if j != i:
                prime = gmul(prime, 1 ^ gmul(Xi_inv, Xj))
        if prime == 0:
            raise RSError("Forney failed")
        y = gmul(Xi, p_eval(err_eval, Xi_inv))
        E[err_pos[i]] = gdiv(y, prime)
    fixed = [a ^ b for a, b in zip(msg, E)]
    if max(rs_syndromes(fixed, n)) != 0:
        raise RSError("could not correct")
    return fixed, sorted(err_pos), synd_list


# --------------------------------------------------------------------------
# 3. QR code tables (versions 1-10)
# --------------------------------------------------------------------------
# (ec codewords per block, [(number of blocks, data codewords per block), ...])
EC_TABLE = {
    1: {'L': (7, [(1, 19)]), 'M': (10, [(1, 16)]), 'Q': (13, [(1, 13)]), 'H': (17, [(1, 9)])},
    2: {'L': (10, [(1, 34)]), 'M': (16, [(1, 28)]), 'Q': (22, [(1, 22)]), 'H': (28, [(1, 16)])},
    3: {'L': (15, [(1, 55)]), 'M': (26, [(1, 44)]), 'Q': (18, [(2, 17)]), 'H': (22, [(2, 13)])},
    4: {'L': (20, [(1, 80)]), 'M': (18, [(2, 32)]), 'Q': (26, [(2, 24)]), 'H': (16, [(4, 9)])},
    5: {'L': (26, [(1, 108)]), 'M': (24, [(2, 43)]), 'Q': (18, [(2, 15), (2, 16)]), 'H': (22, [(2, 11), (2, 12)])},
    6: {'L': (18, [(2, 68)]), 'M': (16, [(4, 27)]), 'Q': (24, [(4, 19)]), 'H': (28, [(4, 15)])},
    7: {'L': (20, [(2, 78)]), 'M': (18, [(4, 31)]), 'Q': (18, [(2, 14), (4, 15)]), 'H': (26, [(4, 13), (1, 14)])},
    8: {'L': (24, [(2, 97)]), 'M': (22, [(2, 38), (2, 39)]), 'Q': (22, [(4, 18), (2, 19)]), 'H': (26, [(4, 14), (2, 15)])},
    9: {'L': (30, [(2, 116)]), 'M': (22, [(3, 36), (2, 37)]), 'Q': (20, [(4, 16), (4, 17)]), 'H': (24, [(4, 12), (4, 13)])},
    10: {'L': (18, [(2, 68), (2, 69)]), 'M': (26, [(4, 43), (1, 44)]), 'Q': (24, [(6, 19), (2, 20)]), 'H': (28, [(6, 15), (2, 16)])},
}
ALIGN = {1: [], 2: [6, 18], 3: [6, 22], 4: [6, 26], 5: [6, 30], 6: [6, 34],
         7: [6, 22, 38], 8: [6, 24, 42], 9: [6, 26, 46], 10: [6, 28, 50]}
ECL_BITS = {'L': 1, 'M': 0, 'Q': 3, 'H': 2}
MAX_VERSION = 10

# module type codes (used by the page to colour the "anatomy" overlay)
T_DATA, T_FINDER, T_SEP, T_TIMING, T_ALIGN, T_FORMAT, T_VERSION, T_DARK = range(8)

MASKS = [
    lambda x, y: (x + y) % 2 == 0,
    lambda x, y: y % 2 == 0,
    lambda x, y: x % 3 == 0,
    lambda x, y: (x + y) % 3 == 0,
    lambda x, y: (x // 3 + y // 2) % 2 == 0,
    lambda x, y: x * y % 2 + x * y % 3 == 0,
    lambda x, y: (x * y % 2 + x * y % 3) % 2 == 0,
    lambda x, y: ((x + y) % 2 + x * y % 3) % 2 == 0,
]
MASK_FORMULAS = [
    "(x + y) mod 2 = 0", "y mod 2 = 0", "x mod 3 = 0", "(x + y) mod 3 = 0",
    "(⌊x/3⌋ + ⌊y/2⌋) mod 2 = 0", "xy mod 2 + xy mod 3 = 0",
    "(xy mod 2 + xy mod 3) mod 2 = 0", "((x+y) mod 2 + xy mod 3) mod 2 = 0",
]


def format_bits(ecl, mask):
    data = ECL_BITS[ecl] << 3 | mask
    rem = data
    for _ in range(10):
        rem = (rem << 1) ^ ((rem >> 9) * 0x537)
    return ((data << 10) | (rem & 0x3FF)) ^ 0x5412


def version_bits(ver):
    rem = ver
    for _ in range(12):
        rem = (rem << 1) ^ ((rem >> 11) * 0x1F25)
    return ver << 12 | (rem & 0xFFF)


def format_positions(size):
    """(x, y) of format bit i for both copies."""
    c1 = [(8, i) for i in range(6)] + [(8, 7), (8, 8), (7, 8)] + [(14 - i, 8) for i in range(9, 15)]
    c2 = [(size - 1 - i, 8) for i in range(8)] + [(8, size - 15 + i) for i in range(8, 15)]
    return c1, c2


VALID_FORMATS = [(e, m, format_bits(e, m)) for e in 'LMQH' for m in range(8)]


# --------------------------------------------------------------------------
# 4. Grid construction
# --------------------------------------------------------------------------
_TEMPLATES = {}


def template(ver):
    """Function patterns for a version: (modules, types, zigzag order of data modules)."""
    if ver in _TEMPLATES:
        return _TEMPLATES[ver]
    size = 17 + 4 * ver
    mod = [[0] * size for _ in range(size)]
    typ = [[T_DATA] * size for _ in range(size)]

    def setf(x, y, dark, t):
        mod[y][x] = 1 if dark else 0
        typ[y][x] = t

    for i in range(size):
        setf(6, i, i % 2 == 0, T_TIMING)
        setf(i, 6, i % 2 == 0, T_TIMING)
    for cx, cy in ((3, 3), (size - 4, 3), (3, size - 4)):
        for dy in range(-4, 5):
            for dx in range(-4, 5):
                d = max(abs(dx), abs(dy))
                x, y = cx + dx, cy + dy
                if 0 <= x < size and 0 <= y < size:
                    setf(x, y, d != 2 and d != 4, T_FINDER if d <= 3 else T_SEP)
    pos = ALIGN[ver]
    n = len(pos)
    for i in range(n):
        for j in range(n):
            if (i == 0 and j == 0) or (i == 0 and j == n - 1) or (i == n - 1 and j == 0):
                continue
            for dy in range(-2, 3):
                for dx in range(-2, 3):
                    setf(pos[i] + dx, pos[j] + dy, max(abs(dx), abs(dy)) != 1, T_ALIGN)
    c1, c2 = format_positions(size)
    for x, y in c1 + c2:
        setf(x, y, False, T_FORMAT)
    setf(8, size - 8, True, T_DARK)
    if ver >= 7:
        bits = version_bits(ver)
        for i in range(18):
            a, b = size - 11 + i % 3, i // 3
            dark = (bits >> i) & 1
            setf(a, b, dark, T_VERSION)
            setf(b, a, dark, T_VERSION)

    # zig-zag order of the data modules: pairs of columns, right to left, snaking up/down
    order = []
    right = size - 1
    while right >= 1:
        if right == 6:
            right = 5
        upward = ((right + 1) & 2) == 0
        for vert in range(size):
            y = size - 1 - vert if upward else vert
            for j in range(2):
                x = right - j
                if typ[y][x] == T_DATA:
                    order.append((x, y))
        right -= 2
    _TEMPLATES[ver] = (mod, typ, order)
    return _TEMPLATES[ver]


def block_layout(ver, ecl):
    """List of (data_len, ec_len) per block, and the interleaving origin of each codeword."""
    ec, groups = EC_TABLE[ver][ecl]
    blocks = []
    for cnt, dlen in groups:
        blocks += [(dlen, ec)] * cnt
    origin = []  # interleaved index -> (block, index within block)
    maxd = max(d for d, _ in blocks)
    for i in range(maxd):
        for b, (d, _) in enumerate(blocks):
            if i < d:
                origin.append((b, i))
    for i in range(ec):
        for b, (d, _) in enumerate(blocks):
            origin.append((b, d + i))
    return blocks, origin


def capacity_bytes(ver, ecl):
    blocks, _ = block_layout(ver, ecl)
    total_bits = sum(d for d, _ in blocks) * 8
    cc = 8 if ver < 10 else 16
    return (total_bits - 4 - cc) // 8


def penalty(mod):
    n = len(mod)
    n1 = n2 = n3 = 0
    P1 = [1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0]
    P2 = P1[::-1]
    lines = [list(r) for r in mod] + [[mod[y][x] for y in range(n)] for x in range(n)]
    for line in lines:
        run = 1
        for i in range(1, n):
            if line[i] == line[i - 1]:
                run += 1
            else:
                if run >= 5:
                    n1 += run - 2
                run = 1
        if run >= 5:
            n1 += run - 2
        s = [0] * 4 + line + [0] * 4
        for i in range(len(s) - 10):
            w = s[i:i + 11]
            if w == P1 or w == P2:
                n3 += 40
    for y in range(n - 1):
        r0, r1 = mod[y], mod[y + 1]
        for x in range(n - 1):
            c = r0[x]
            if c == r0[x + 1] == r1[x] == r1[x + 1]:
                n2 += 3
    dark = sum(map(sum, mod))
    pct = dark * 100.0 / (n * n)
    n4 = int(abs(pct - 50) // 5) * 10
    return {'n1': n1, 'n2': n2, 'n3': n3, 'n4': n4, 'total': n1 + n2 + n3 + n4, 'dark_pct': round(pct, 1)}


def place(ver, ecl, codewords, mask):
    """Return a full grid: function patterns + data bits (masked) + format info."""
    tmod, typ, order = template(ver)
    size = len(tmod)
    mod = [row[:] for row in tmod]
    nbits = len(codewords) * 8
    for i, (x, y) in enumerate(order):
        bit = (codewords[i >> 3] >> (7 - (i & 7))) & 1 if i < nbits else 0
        if mask is not None and MASKS[mask](x, y):
            bit ^= 1
        mod[y][x] = bit
    if mask is not None:
        fb = format_bits(ecl, mask)
        c1, c2 = format_positions(size)
        for i in range(15):
            b = (fb >> i) & 1
            mod[c1[i][1]][c1[i][0]] = b
            mod[c2[i][1]][c2[i][0]] = b
    return mod


def encode(text, ecl):
    data = text.encode('utf-8')
    ver = None
    for v in range(1, MAX_VERSION + 1):
        if len(data) <= capacity_bytes(v, ecl):
            ver = v
            break
    if ver is None:
        raise ValueError('Message too long: max %d bytes at level %s' % (capacity_bytes(MAX_VERSION, ecl), ecl))
    blocks, origin = block_layout(ver, ecl)
    cap = sum(d for d, _ in blocks)
    bits = []

    def put(val, n):
        for i in range(n - 1, -1, -1):
            bits.append((val >> i) & 1)
    put(0b0100, 4)                      # mode indicator: byte mode
    put(len(data), 8 if ver < 10 else 16)  # character count
    for b in data:
        put(b, 8)
    put(0, min(4, cap * 8 - len(bits)))  # terminator
    while len(bits) % 8:
        bits.append(0)
    dcw = [int(''.join(map(str, bits[i:i + 8])), 2) for i in range(0, len(bits), 8)]
    pad = [0xEC, 0x11]
    i = 0
    while len(dcw) < cap:
        dcw.append(pad[i % 2])
        i += 1
    # split into blocks and compute EC
    blk = []
    k = 0
    for d, ec in blocks:
        part = dcw[k:k + d]
        k += d
        blk.append({'data': part, 'ec': rs_encode(part, ec)})
    inter = []
    for b, idx in origin:
        d = blocks[b][0]
        inter.append(blk[b]['data'][idx] if idx < d else blk[b]['ec'][idx - d])
    return ver, blocks, origin, blk, inter, data


# --------------------------------------------------------------------------
# 5. Decoder: grid of 0/1 -> text, with a full report of every stage
# --------------------------------------------------------------------------
FINDER_LOST_AT = 8  # a finder with this many of its 49 modules wrong no longer looks like 1:1:3:1:1


def decode(mod):
    size = len(mod)
    ver = (size - 17) // 4
    tmod, typ, order = template(ver)
    rep = {'stage': 'locate', 'ok': False, 'version': ver}

    # 1. finders
    finders = []
    for cx, cy in ((3, 3), (size - 4, 3), (3, size - 4)):
        bad = 0
        for dy in range(-3, 4):
            for dx in range(-3, 4):
                if mod[cy + dy][cx + dx] != tmod[cy + dy][cx + dx]:
                    bad += 1
        finders.append({'bad': bad, 'found': bad < FINDER_LOST_AT})
    rep['finders'] = finders
    if not all(f['found'] for f in finders):
        return rep

    # 2. format info: pick the nearest valid 15-bit pattern from either copy
    rep['stage'] = 'format'
    c1, c2 = format_positions(size)
    best = None
    copies = []
    for copy in (c1, c2):
        val = 0
        for i, (x, y) in enumerate(copy):
            val |= mod[y][x] << i
        dist, e, m = min((bin(val ^ fb).count('1'), e, m) for e, m, fb in VALID_FORMATS)
        copies.append({'dist': dist})
        if best is None or dist < best[0]:
            best = (dist, e, m)
    rep['format'] = {'copies': copies, 'dist': best[0], 'ecl': best[1], 'mask': best[2]}
    if best[0] > 3:
        return rep
    _, ecl, mask = best

    # 3. unmask + read codewords in zig-zag order
    rep['stage'] = 'fix'
    blocks, origin = block_layout(ver, ecl)
    total = len(origin)
    raw = [0] * total
    for i, (x, y) in enumerate(order[:total * 8]):
        bit = mod[y][x] ^ (1 if MASKS[mask](x, y) else 0)
        raw[i >> 3] |= bit << (7 - (i & 7))
    recv = [[0] * (d + e) for d, e in blocks]
    for i, (b, idx) in enumerate(origin):
        recv[b][idx] = raw[i]

    # 4. Reed-Solomon on every block
    out_blocks = []
    data = []
    all_ok = True
    for b, (d, e) in enumerate(blocks):
        r = recv[b]
        info = {'data_len': d, 'ec_len': e, 'can_fix': e // 2, 'read': r,
                'syndromes': rs_syndromes(r, e)}
        try:
            fixed, pos, _ = rs_decode(r, e)
            info.update(status='fixed' if pos else 'clean', positions=pos, corrected=fixed)
            data += fixed[:d]
        except RSError:
            info.update(status='failed', positions=[], corrected=r)
            all_ok = False
        out_blocks.append(info)
    rep['blocks'] = out_blocks
    rep['ecl'] = ecl
    rep['mask'] = mask
    if not all_ok:
        return rep

    # 5. parse the bit stream
    rep['stage'] = 'parse'
    bits = ''.join('{:08b}'.format(c) for c in data)
    mode = int(bits[:4], 2)
    if mode != 0b0100:
        rep['error'] = 'unexpected mode %d' % mode
        return rep
    cc = 8 if ver < 10 else 16
    count = int(bits[4:4 + cc], 2)
    start = 4 + cc
    if start + count * 8 > len(bits):
        rep['error'] = 'length field is corrupt'
        return rep
    payload = bytes(int(bits[start + 8 * i:start + 8 * i + 8], 2) for i in range(count))
    rep['text'] = payload.decode('utf-8', errors='replace')
    rep['stage'] = 'done'
    rep['ok'] = True
    return rep


# --------------------------------------------------------------------------
# 6. JSON API used by the web page
# --------------------------------------------------------------------------
def _flat(mod):
    return ''.join(''.join(str(v) for v in row) for row in mod)


def api_encode(text, ecl, force_mask=-1):
    try:
        ver, blocks, origin, blk, inter, data = encode(text, ecl)
    except ValueError as e:
        return json.dumps({'error': str(e)})
    tmod, typ, order = template(ver)
    size = len(tmod)
    candidates = []
    grids = []
    for m in range(8):
        g = place(ver, ecl, inter, m)
        grids.append(g)
        candidates.append(penalty(g))
    auto = min(range(8), key=lambda m: candidates[m]['total'])
    mask = auto if force_mask < 0 else force_mask
    unmasked = place(ver, ecl, inter, None)
    # which codeword each module belongs to (-1 = function pattern / remainder bit)
    cwmap = [-1] * (size * size)
    for i, (x, y) in enumerate(order):
        if i < len(inter) * 8:
            cwmap[y * size + x] = i >> 3
    return json.dumps({
        'version': ver, 'size': size, 'ecl': ecl, 'mask': mask, 'auto_mask': auto,
        'modules': _flat(grids[mask]),
        'types': ''.join(''.join(str(v) for v in row) for row in typ),
        'cwmap': cwmap,
        'origin': origin,
        'blocks': [{'data': b['data'], 'ec': b['ec']} for b in blk],
        'penalties': candidates,
        'mask_grids': [_flat(g) for g in grids],
        'unmasked': _flat(unmasked),
        'unmasked_penalty': penalty(unmasked),
        'mask_formulas': MASK_FORMULAS,
        'payload_bytes': len(data),
        'capacity': capacity_bytes(ver, ecl),
        'generator_exps': [LOG[c] for c in rs_generator(blocks[0][1])],
        'ec_per_block': blocks[0][1],
    })


def api_decode(flat, size):
    mod = [[int(flat[y * size + x]) for x in range(size)] for y in range(size)]
    return json.dumps(decode(mod))
