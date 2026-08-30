"""Streaming parser for R serialization format v3 (XDR/big-endian).
Reference: R serialize.c v3 semantics (WriteItem/ReadItem/PackFlags).
Only metadata is materialized; large numeric payloads are skipped.
"""
import struct
import sys

sys.setrecursionlimit(500000)

MAX_STRVEC = 60000   # materialize string vectors up to this length
MAX_NUM = 10000      # materialize numeric vectors up to this length

PAIR_TYPES = (2, 3, 5, 6, 17)   # LISTSXP CLOSXP PROMSXP LANGSXP DOTSXP
NO_PAYLOAD = (241, 242, 250, 251, 252, 253, 254)


class RdsReader:
    def __init__(self, path):
        self.f = open(path, "rb")
        self.refs = []
        magic = self.f.read(2)
        if magic != b"X\n":
            raise ValueError("bad magic %r" % magic)
        version, writer, minr = struct.unpack(">iii", self.f.read(12))
        nelen = struct.unpack(">i", self.f.read(4))[0]
        enc = self.f.read(nelen).decode("ascii", "replace")
        self.header = {"version": version, "writer_version": writer,
                       "min_reader": minr, "encoding": enc}

    def i32(self):
        b = self.f.read(4)
        if len(b) < 4:
            raise EOFError("unexpected EOF at byte %d" % self.f.tell())
        return struct.unpack(">i", b)[0]

    def skip(self, n):
        self.f.seek(n, 1)

    def read_length(self):
        ln = self.i32()
        if ln == -1:  # long vector: two 32-bit halves
            hi = self.i32()
            lo = self.i32() & 0xFFFFFFFF
            return hi * (1 << 32) + lo
        return ln

    def read_item(self):
        flags = self.i32()
        return self.body(flags)

    def body(self, flags):
        t = flags & 0xFF
        levs = flags >> 12
        hasattr = bool(flags & 0x200)

        if t in NO_PAYLOAD:
            return ("special", t)
        if t == 255:  # REFSXP
            idx = flags >> 8
            if idx == 0:
                idx = self.i32()
            if idx < 1 or idx > len(self.refs):
                raise ValueError("ref index %d out of range (%d refs) at %d"
                                 % (idx, len(self.refs), self.f.tell()))
            return self.refs[idx - 1]
        if t == 1:  # SYMSXP
            name = self.read_item()
            tok = ("sym", name)
            self.refs.append(tok)
            return tok
        if t == 4:  # ENVSXP
            locked = self.i32()
            tok = ("env", locked)
            self.refs.append(tok)
            self.read_item(); self.read_item(); self.read_item(); self.read_item()
            return tok
        if t in PAIR_TYPES:
            cells = []
            cur = flags
            while (cur & 0xFF) in PAIR_TYPES:
                a = self.read_item() if (cur & 0x200) else None
                tg = self.read_item() if (cur & 0x400) else None
                car = self.read_item()
                cells.append((tg, car, a))
                cur = self.i32()
            tail = self.body(cur)
            return ("pairlist", cells, tail)
        if t in (247, 248, 249):  # PERSISTSXP PACKAGESXP NAMESPACESXP
            if self.i32() != 0:
                raise ValueError("persistent string names unsupported")
            ln = self.read_length()
            items = [self.read_item() for _ in range(ln)]
            tok = ("strvec", items)
            self.refs.append(tok)
            return tok
        if t == 238:  # ALTREP
            info = self.read_item()
            state = self.read_item()
            at = self.read_item()
            return ("altrep", info, state, at)
        if t == 22:  # EXTPTRSXP
            tok = ("xptr",)
            self.refs.append(tok)
            self.read_item(); self.read_item()
            if hasattr:
                self.read_item()
            return tok
        if t == 23:  # WEAKREFSXP
            tok = ("wref",)
            self.refs.append(tok)
            if hasattr:
                self.read_item()
            return tok
        if t in (7, 8):  # SPECIALSXP BUILTINSXP
            ln = self.i32()
            self.skip(ln)
            if hasattr:
                self.read_item()
            return ("builtin",)
        if t == 9:  # CHARSXP
            ln = self.i32()
            if ln < 0:
                s = None
            else:
                b = self.f.read(ln)
                if levs & 8:
                    s = b.decode("utf-8", "replace")
                elif levs & 4:
                    s = b.decode("latin-1", "replace")
                else:
                    s = b.decode("ascii", "replace")
            if hasattr:
                self.read_item()  # read and ignore (legacy)
            return s
        if t in (10, 13):  # LGLSXP INTSXP
            ln = self.read_length()
            if ln <= MAX_NUM:
                b = self.f.read(4 * ln)
                tok = ("intvec", list(struct.unpack(">%di" % ln, b)), None)
            else:
                self.skip(4 * ln)
                tok = ("intvec_skip", ln)
            if hasattr:
                at = self.read_item()
                tok = (tok[0], tok[1], at)
            return tok
        if t == 14:  # REALSXP
            ln = self.read_length()
            if ln <= MAX_NUM:
                b = self.f.read(8 * ln)
                tok = ("realvec", list(struct.unpack(">%dd" % ln, b)), None)
            else:
                self.skip(8 * ln)
                tok = ("realvec_skip", ln)
            if hasattr:
                at = self.read_item()
                tok = (tok[0], tok[1], at)
            return tok
        if t == 15:  # CPLXSXP
            ln = self.read_length()
            self.skip(16 * ln)
            tok = ("cplxvec_skip", ln)
            if hasattr:
                at = self.read_item()
                tok = (tok[0], tok[1], at)
            return tok
        if t == 16:  # STRSXP
            ln = self.read_length()
            if ln <= MAX_STRVEC:
                items = [self.read_item() for _ in range(ln)]
                tok = ("strvec", items, None)
            else:
                head = [self.read_item() for _ in range(3)]
                for _ in range(ln - 3):
                    self.read_item()
                tok = ("strvec_skip", ln, head)
            if hasattr:
                at = self.read_item()
                tok = (tok[0], tok[1], at)
            return tok
        if t in (19, 20):  # VECSXP EXPRSXP
            ln = self.read_length()
            items = [self.read_item() for _ in range(ln)]
            at = self.read_item() if hasattr else None
            return ("vec", items, at)
        if t == 24:  # RAWSXP
            ln = self.read_length()
            self.skip(ln)
            tok = ("raw_skip", ln)
            if hasattr:
                at = self.read_item()
                tok = (tok[0], tok[1], at)
            return tok
        if t == 25:  # OBJSXP / S4SXP — payload empty, slots are attributes
            at = self.read_item() if hasattr else None
            return ("obj", at)
        if t == 21:
            raise ValueError("BCODESXP not handled (byte %d)" % self.f.tell())
        raise ValueError("unknown type %d at byte %d" % (t, self.f.tell()))
