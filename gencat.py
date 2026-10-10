#!/usr/bin/env python3
"""Build-time tool: normalize the object CSVs into nightsky.dat, the compact
binary catalogue that nightsky.exe loads at runtime.

Inputs  : catalog_raw.csv  (Messier catalogue)
          catalog_ngc.csv  (curated bright non-Messier NGC/IC + Caldwell)
Output  : nightsky.dat     (little-endian packed binary; see format below)

The .dat layout mirrors the C++ CatEntry / DatHeader in nightsky_types.h /
nightsky.cpp EXACTLY, so the program can load it with a single fread. If you
ever change the record layout there, bump DAT_VERSION in both places.

  Header : magic[4]="NSKY", uint32 version, uint32 count      (12 bytes)
  Record : id[8], ngc[12], name[40], int32 type, float ra_h,
           float dec_deg, float mag, float size_arcmin         (80 bytes)
"""
import csv, re, struct, sys

DAT_VERSION = 2
DAT_PATH    = "nightsky.dat"
INPUTS      = ["catalog_raw.csv", "catalog_ngc.csv"]

REC_FMT = "<8s12s40si4f"   # must equal sizeof(CatEntry)==80
HDR_FMT = "<4sII"          # must equal sizeof(DatHeader)==12
assert struct.calcsize(REC_FMT) == 80
assert struct.calcsize(HDR_FMT) == 12

def num(s):
    m = re.search(r'[-+]?\d+(?:\.\d+)?', s)
    return float(m.group()) if m else 0.0

def parse_ra_hours(s):
    h   = re.search(r'([\d.]+)\s*h', s)
    m   = re.search(r'([\d.]+)\s*m', s)
    sec = re.search(r'([\d.]+)\s*s', s)
    val = float(h.group(1)) if h else 0.0
    if m:   val += float(m.group(1)) / 60.0
    if sec: val += float(sec.group(1)) / 3600.0
    return val

def parse_dec_deg(s):
    neg = ('-' in s) or ('−' in s)          # ascii or unicode minus
    parts = re.findall(r'[\d.]+', s)
    d = float(parts[0]) if len(parts) > 0 else 0.0
    if len(parts) > 1: d += float(parts[1]) / 60.0
    if len(parts) > 2: d += float(parts[2]) / 3600.0
    return -d if neg else d

def parse_size_arcmin(s):
    s = s.replace('×', 'x').replace('X', 'x')
    best = 0.0
    for tok in s.split('x'):
        tok = tok.strip()
        if not tok:
            continue
        v = num(tok)
        if '°' in tok:                        v *= 60.0   # degrees
        elif '″' in tok or '"' in tok:        v /= 60.0   # arcsec
        best = max(best, v)
    return best

def type_code(t):
    t = t.lower()
    if 'galaxy' in t:                                   return 0  # GALAXY
    if 'nebula' in t or 'remnant' in t or 'h ii' in t:  return 1  # NEBULA
    if 'globular' in t:                                 return 2  # GLOBULAR
    if 'open' in t:                                     return 3  # OPEN CLUSTER
    return 4                                                      # OTHER

def field(s, n):
    """Encode a string into an n-byte, null-terminated field (truncate safely)."""
    b = s.encode('ascii', 'replace')[:n-1]
    return b

def main():
    rows = []
    seen = set()
    for path in INPUTS:
        try:
            f = open(path, encoding='utf-8')
        except FileNotFoundError:
            print("  (skipping missing %s)" % path)
            continue
        with f:
            for r in csv.DictReader(f):
                mid = r['Messier'].strip()
                if not mid or mid in seen:
                    continue
                seen.add(mid)
                rows.append((
                    mid,
                    r['NGC/IC'].strip().replace(' ', ''),
                    r['Common Name'].strip(),
                    type_code(r['Object Type']),
                    parse_ra_hours(r['RA']),
                    parse_dec_deg(r['Dec']),
                    num(r['Mag']),
                    parse_size_arcmin(r['Size']),
                ))

    with open(DAT_PATH, 'wb') as o:
        o.write(struct.pack(HDR_FMT, b'NSKY', DAT_VERSION, len(rows)))
        for (mid, ngc, name, tc, ra, dec, mag, size) in rows:
            o.write(struct.pack(REC_FMT,
                                field(mid, 8), field(ngc, 12), field(name, 40),
                                tc, ra, dec, mag, size))
    print("Wrote %s: %d objects, %d bytes." %
          (DAT_PATH, len(rows), 12 + 80 * len(rows)))

if __name__ == '__main__':
    main()
