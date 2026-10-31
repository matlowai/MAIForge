#!/usr/bin/env python3
"""Independent check of golden/keyed_rng_v1.golden.json.

A second implementation of vz.keyed-rng/v1 written from ADR 044 section 2 and
the crate docs, with its own pure-Python BLAKE3 (single-chunk inputs only,
<= 1024 bytes, which covers every segment and draw message). It recomputes the
seed roots, every named scope key, and every draw row (u64, uniform bits,
below, range, chance, pick, normal_q32 with the table parsed from
src/normal_table.rs). Exit 0 and a summary line on success.

Usage: python3 -I tools/verify_golden_py.py <crate dir>
"""
import json
import re
import struct
import sys

IV = [0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A, 0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19]
PERM = [2, 6, 3, 10, 7, 0, 4, 13, 1, 11, 12, 5, 9, 14, 15, 8]
CHUNK_START, CHUNK_END, ROOT, KEYED_HASH, DK_CONTEXT, DK_MATERIAL = 1, 2, 8, 16, 32, 64
M32 = 0xFFFFFFFF


def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & M32


def g(s, a, b, c, d, x, y):
    s[a] = (s[a] + s[b] + x) & M32; s[d] = rotr(s[d] ^ s[a], 16)
    s[c] = (s[c] + s[d]) & M32; s[b] = rotr(s[b] ^ s[c], 12)
    s[a] = (s[a] + s[b] + y) & M32; s[d] = rotr(s[d] ^ s[a], 8)
    s[c] = (s[c] + s[d]) & M32; s[b] = rotr(s[b] ^ s[c], 7)


def compress(cv, m, counter, blen, flags):
    s = list(cv) + IV[:4] + [counter & M32, (counter >> 32) & M32, blen, flags]
    m = list(m)
    for r in range(7):
        g(s, 0, 4, 8, 12, m[0], m[1]); g(s, 1, 5, 9, 13, m[2], m[3])
        g(s, 2, 6, 10, 14, m[4], m[5]); g(s, 3, 7, 11, 15, m[6], m[7])
        g(s, 0, 5, 10, 15, m[8], m[9]); g(s, 1, 6, 11, 12, m[10], m[11])
        g(s, 2, 7, 8, 13, m[12], m[13]); g(s, 3, 4, 9, 14, m[14], m[15])
        if r < 6:
            m = [m[p] for p in PERM]
    for i in range(8):
        s[i] ^= s[i + 8]
        s[i + 8] ^= cv[i]
    return s


def words(b):
    return list(struct.unpack('<16I', b.ljust(64, b'\0')))


def blake3_xof(data, key_words, flags, out_len):
    assert len(data) <= 1024, 'single-chunk reference only'
    blocks = [data[i:i + 64] for i in range(0, len(data), 64)] or [b'']
    cv = list(key_words)
    for i, blk in enumerate(blocks[:-1]):
        cv = compress(cv, words(blk), 0, 64, flags | (CHUNK_START if i == 0 else 0))[:8]
    last = blocks[-1]
    fl = flags | (CHUNK_START if len(blocks) == 1 else 0) | CHUNK_END | ROOT
    out = b''
    ctr = 0
    while len(out) < out_len:
        out += struct.pack('<16I', *compress(cv, words(last), ctr, len(last), fl))
        ctr += 1
    return out[:out_len]


def kw(key):
    return list(struct.unpack('<8I', key))


def keyed(key, data, n=32):
    return blake3_xof(data, kw(key), KEYED_HASH, n)


def derive_key(context, material):
    ck = blake3_xof(context.encode(), IV, DK_CONTEXT, 32)
    return blake3_xof(material, kw(ck), DK_MATERIAL, 32)


# ---- vz.keyed-rng/v1, from the ADR and crate docs ----
def lp(s):
    b = s.encode()
    return struct.pack('<I', len(b)) + b


def seg_str(tag, s):
    return bytes([tag]) + lp(s)


def root(seed):
    return derive_key('vz.keyed-rng/v1', struct.pack('<Q', seed))


def child(key, seg):
    return keyed(key, seg)


def msg(purpose, counter, attempt):
    return b'\x06' + lp(purpose) + b'\x07' + struct.pack('<Q', counter) + b'\x08' + struct.pack('<I', attempt)


class Draw:
    def __init__(self, key, purpose, counter, attempt):
        self.buf = keyed(key, msg(purpose, counter, attempt), 64 * 4)
        self.i = 0

    def word(self):
        w = struct.unpack('<Q', self.buf[self.i:self.i + 8])[0]
        self.i += 8
        return w


def below(d, n):
    m = d.word() * n
    if (m & ((1 << 64) - 1)) < n:
        t = ((1 << 64) - n) % n
        while (m & ((1 << 64) - 1)) < t:
            m = d.word() * n
    return m >> 64


def lerp40(dv, f):
    return (dv * f) >> 40  # Python >> on negative ints is floor, like i128 >>


def normal(d, coarse, tail):
    w0 = d.word()
    c = w0 >> 52
    f = (w0 >> 12) & ((1 << 40) - 1)
    if c == 0 or c == 4095:
        w1 = d.word()
        s = w1 >> 52
        f2 = (w1 >> 12) & ((1 << 40) - 1)
        z = tail[s] + lerp40(tail[s + 1] - tail[s], f2)
        return z if c == 0 else -z
    return coarse[c - 1] + lerp40(coarse[c] - coarse[c - 1], f)


def pick(d, weights):
    r = below(d, sum(weights))
    for i, w in enumerate(weights):
        if r < w:
            return i
        r -= w


def main():
    crate = sys.argv[1]
    assert blake3_xof(b'', IV, 0, 32).hex() == 'af1349b9f5f9a1a6a0404dea36dcc9499bcb25c9adc112b7cc9a93cae41f3262'
    gold = json.load(open(f'{crate}/golden/keyed_rng_v1.golden.json'))
    tab = open(f'{crate}/src/normal_table.rs').read()
    def arr(name):
        body = re.search(name + r': \[i64; \d+\] = \[(.*?)\];', tab, re.S).group(1)
        return [int(x) for x in body.replace('\n', ' ').split(',') if x.strip()]
    coarse, tail = arr('COARSE'), arr('TAIL')
    n = 0
    for row in gold['seeds']:
        assert root(int(row['seed'], 16)).hex() == row['root_key'], row; n += 1
    # Rebuild the named scopes exactly as the docs specify.
    gl = 'blake3:' + blake3_xof(b'vz.keyed-rng golden layer', IV, 0, 32).hex()
    r = root(42)
    law = child(r, seg_str(1, 'laws.golden/v1'))
    lay = child(law, seg_str(2, gl))
    ch = child(lay, b'\x03' + struct.pack('<I', 0))
    tile = lambda k, lod, ix, iy: child(k, b'\x04' + bytes([lod]) + struct.pack('<qq', ix, iy))
    t = tile(ch, 0, -12, 7)
    ad = child(t, seg_str(5, 't/3.5/12/g0'))
    leg = child(child(root(7), seg_str(0x40, 'preview-genesis')), seg_str(0x41, 'salt-1'))
    scopes = {
        's42.root': r, 's42.law': law, 's42.law.layer': lay, 's42.law.layer.chart0': ch,
        's42.law.layer.chart0.tile0.-12.7': t, 's42.law.layer.chart0.tile0.-12.7.addr': ad,
        's42.law.layer.tile0.-12.7(no chart)': tile(lay, 0, -12, 7),
        's42.law.layer.chart0.tile6.0.-1': tile(ch, 6, 0, -1),
        'legacy.7.preview-genesis.salt-1': leg,
    }
    for row in gold['scopes']:
        assert scopes[row['name']].hex() == row['key'], row['name']; n += 1
    for row in gold['draws']:
        key = scopes[row['scope']]
        a = (row['purpose'], int(row['counter'], 16), row['attempt'])
        D = lambda: Draw(key, *a)
        w = D().word()
        assert f'{w:016x}' == row['u64'], row
        u = (w >> 11) / 9007199254740992.0
        assert struct.pack('>d', u).hex() == row['uniform_f64_bits'], row
        assert below(D(), 1000) == row['below_1000'], row
        assert f'{below(D(), (1 << 63) + 1):016x}' == row['below_2p63_plus_1'], row
        assert -50 + below(D(), 100) == row['range_m50_50'], row
        assert int(below(D(), 3) < 1) == row['chance_1_3'], row
        assert pick(D(), [1, 2, 3, 4]) == row['pick_1234'], row
        assert normal(D(), coarse, tail) == row['normal_q32'], row
        n += 1
    nt = gold['normal_tail']
    for side in ('lower', 'upper'):
        for row in nt[side]:
            assert normal(Draw(scopes[nt['scope']], nt['purpose'], row['counter'], 0), coarse, tail) == row['normal_q32']; n += 1
    print(f'python reference: {n} golden rows match (independent BLAKE3 + vz.keyed-rng/v1)')


if __name__ == '__main__':
    main()
