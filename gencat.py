#!/usr/bin/env python3
"""Build-time tool: merge the object sources into nightsky.dat, the compact
binary catalogue that nightsky.exe loads at runtime.

Inputs
  catalog_raw.csv        Messier catalogue                    (override names)
  catalog_ngc.csv        curated bright NGC/IC/Caldwell list  (override names)
  catalog_groups.csv     named galaxy groups; position, magnitude and size are
                         computed from the member galaxies in OpenNGC
  openngc_NGC.csv        OpenNGC NGC/IC database      } https://github.com/
  openngc_addendum.csv   OpenNGC addendum objects     } mattiaverga/OpenNGC
                         (CC BY-SA 4.0, (c) Mattia Verga)

Objects in the two override CSVs win over their OpenNGC duplicates; their
positions are cross-checked against OpenNGC and a warning is printed when they
disagree by more than POS_TOLERANCE_DEG or a quarter of the object's size,
whichever is larger (the "centre" of a large, diffuse object is ill-defined).

Output: nightsky.dat (little-endian). Must match nightsky_types.h exactly; if
the layout changes, bump DAT_VERSION here and kDatVersion in nightsky.cpp.
  Header : magic[4]="NSKY", u32 version, u32 count, u32 strBytes   (16 bytes)
  Record : f32 ra_hours, f32 dec_deg, f32 mag, f32 size_arcmin,
           u32 id, u32 name (byte offsets into the string pool), u32 type
                                                                   (28 bytes)
  Pool   : strBytes of NUL-terminated ASCII strings; offset 0 is "".
"""
import csv, math, re, struct, unicodedata

DAT_VERSION = 3
DAT_PATH    = "nightsky.dat"
HDR_FMT     = "<4sIII"
REC_FMT     = "<4f3I"
assert struct.calcsize(HDR_FMT) == 16
assert struct.calcsize(REC_FMT) == 28

OVERRIDE_CSVS = ["catalog_raw.csv", "catalog_ngc.csv"]
GROUPS_CSV    = "catalog_groups.csv"
OPENNGC_CSVS  = ["openngc_NGC.csv", "openngc_addendum.csv"]

# Selection: galaxies need a known magnitude at or below this limit (roughly
# what a 1000-2000mm suburban rig can pick up with long stacked exposures).
# Nebulae and clusters are kept regardless, since many have no catalogued
# magnitude at all.
GALAXY_MAG_LIMIT  = 15.0
MAG_UNKNOWN       = 99.0          # sentinel; nightsky prints "--"
POS_TOLERANCE_DEG = 0.25

# Type codes; must match kTypeName[] in nightsky.cpp.
T_GALAXY, T_NEBULA, T_GLOBULAR, T_OPEN, T_OTHER, T_GROUP = range(6)
OPENNGC_TYPES = {
    'G': T_GALAXY, 'GPair': T_GALAXY, 'GTrpl': T_GALAXY, 'GGroup': T_GALAXY,
    'PN': T_NEBULA, 'HII': T_NEBULA, 'EmN': T_NEBULA, 'Neb': T_NEBULA,
    'RfN': T_NEBULA, 'SNR': T_NEBULA, 'Cl+N': T_NEBULA, 'DrkN': T_NEBULA,
    'GCl': T_GLOBULAR, 'OCl': T_OPEN, '*Ass': T_OPEN,
}   # anything else (Dup, *, **, Other, NonEx, Nova) is skipped

# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------
def num(s, default=MAG_UNKNOWN):
    m = re.search(r'[-+]?\d+(?:\.\d+)?', s or '')
    return float(m.group()) if m else default

def parse_ra_hours(s):            # "05h 34m 31.9s" / "17h 40.1m" / "05:34:31.9"
    if ':' in s:
        p = [float(x) for x in s.split(':')]
        return p[0] + p[1] / 60.0 + (p[2] if len(p) > 2 else 0) / 3600.0
    h   = re.search(r'([\d.]+)\s*h', s)
    m   = re.search(r'([\d.]+)\s*m', s)
    sec = re.search(r'([\d.]+)\s*s', s)
    val = float(h.group(1)) if h else 0.0
    if m:   val += float(m.group(1)) / 60.0
    if sec: val += float(sec.group(1)) / 3600.0
    return val

def parse_dec_deg(s):             # "+22° 00′ 52.2″" / "−32° 13′" / "-32:13:00"
    neg = ('-' in s) or ('−' in s)
    parts = re.findall(r'[\d.]+', s)
    d = float(parts[0]) if len(parts) > 0 else 0.0
    if len(parts) > 1: d += float(parts[1]) / 60.0
    if len(parts) > 2: d += float(parts[2]) / 3600.0
    return -d if neg else d

def parse_size_arcmin(s):         # largest axis of "8.0′ × 5.6′", "420″", "3.17°"
    s = s.replace('×', 'x').replace('X', 'x')
    best = 0.0
    for tok in s.split('x'):
        tok = tok.strip()
        if not tok:
            continue
        v = num(tok, 0.0)
        if '°' in tok:                 v *= 60.0
        elif '″' in tok or '"' in tok: v /= 60.0
        best = max(best, v)
    return best

def override_type(t):
    t = t.lower()
    if 'galaxy' in t:                                   return T_GALAXY
    if 'nebula' in t or 'remnant' in t or 'h ii' in t:  return T_NEBULA
    if 'globular' in t:                                 return T_GLOBULAR
    if 'open' in t:                                     return T_OPEN
    return T_OTHER

def key(name):
    """Lookup key with zero padding and spaces removed: 'NGC0224'->'NGC224'."""
    name = name.replace(' ', '').upper()
    m = re.match(r'([A-Z]+)0*(\d.*)$', name)
    return m.group(1) + m.group(2) if m else name

def display_id(name):
    """Conventional short form: drop zero padding for NGC/IC/C/B/Mel ids."""
    m = re.match(r'(NGC|IC|C|B|Mel)0*(\d.*)$', name)
    return m.group(1) + m.group(2) if m else name

def ascii_text(s):
    return unicodedata.normalize('NFKD', s).encode('ascii', 'ignore').decode()

def sep_deg(ra1, dec1, ra2, dec2):
    a1, a2 = math.radians(ra1 * 15), math.radians(ra2 * 15)
    d1, d2 = math.radians(dec1), math.radians(dec2)
    c = (math.sin(d1) * math.sin(d2) +
         math.cos(d1) * math.cos(d2) * math.cos(a1 - a2))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))

def combined_mag(mags):
    known = [m for m in mags if m < MAG_UNKNOWN]
    if not known:
        return MAG_UNKNOWN
    return -2.5 * math.log10(sum(10 ** (-0.4 * m) for m in known))

# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------
def load_openngc():
    """All OpenNGC rows (any type) as dicts, plus a key->row lookup."""
    rows, lookup = [], {}
    for path in OPENNGC_CSVS:
        with open(path, encoding='utf-8') as f:
            for r in csv.DictReader(f, delimiter=';'):
                if not r['RA'] or not r['Dec']:
                    continue
                v, b = num(r['V-Mag']), num(r['B-Mag'])
                o = {
                    'name':  r['Name'],
                    'type':  r['Type'],
                    'ra':    parse_ra_hours(r['RA']),
                    'dec':   parse_dec_deg(r['Dec']),
                    'mag':   v if v < MAG_UNKNOWN else b,   # V preferred, else B
                    'size':  num(r['MajAx'], 0.0),
                    'common': r['Common names'].split(',')[0].strip(),
                    'idents': r['Identifiers'],
                    'xref':  [key(r['Name'])],
                }
                if r['M']:
                    o['xref'].append('M' + str(int(r['M'])))
                for col, pre in (('NGC', 'NGC'), ('IC', 'IC')):
                    for x in filter(None, r[col].split(',')):
                        o['xref'].append(key(pre + x.strip()))
                rows.append(o)
                if o['type'] != 'Dup':
                    for k in o['xref']:
                        lookup.setdefault(k, o)
    return rows, lookup

def alt_identifier(idents):
    """Most useful secondary id for an unnamed object: its Caldwell number,
    else UGC / PGC / MCG, with zero padding dropped ('UGC 06439'->'UGC 6439')."""
    ids = [s.strip() for s in idents.split(',') if s.strip()]
    for pre, label in (('C ', 'Caldwell '), ('UGC ', 'UGC '), ('PGC ', 'PGC '),
                       ('MCG ', 'MCG ')):
        for s in ids:
            if s.startswith(pre):
                rest = s[len(pre):]
                return label + ((rest.lstrip('0') or '0') if rest.isdigit() else rest)
    return ''

def load_overrides():
    out = []
    for path in OVERRIDE_CSVS:
        with open(path, encoding='utf-8') as f:
            for r in csv.DictReader(f):
                mid = r['Messier'].strip()
                if not mid:
                    continue
                ngc = r['NGC/IC'].strip().replace(' ', '')
                out.append({
                    'id':   mid,
                    'name': r['Common Name'].strip() or ngc,
                    'type': override_type(r['Object Type']),
                    'ra':   parse_ra_hours(r['RA']),
                    'dec':  parse_dec_deg(r['Dec']),
                    'mag':  num(r['Mag']),
                    'size': parse_size_arcmin(r['Size']),
                    'keys': [key(mid)] + ([key(ngc)] if ngc else []),
                })
    return out

def load_groups(lookup):
    groups, member_of = [], {}
    with open(GROUPS_CSV, encoding='utf-8') as f:
        for r in csv.DictReader(f):
            names = r['Members'].split()
            members = [lookup.get(key(n)) for n in names]
            missing = [n for n, m in zip(names, members) if m is None]
            if missing:
                print("  WARNING: group %s skipped, unknown members: %s"
                      % (r['ID'], ' '.join(missing)))
                continue
            # Centroid as the normalized mean of unit vectors (safe across 0h).
            x = y = z = 0.0
            for m in members:
                a, d = math.radians(m['ra'] * 15), math.radians(m['dec'])
                x += math.cos(d) * math.cos(a)
                y += math.cos(d) * math.sin(a)
                z += math.sin(d)
            ra  = (math.degrees(math.atan2(y, x)) / 15.0) % 24.0
            dec = math.degrees(math.atan2(z, math.hypot(x, y)))
            # Size: diameter of the circle about the centroid enclosing members.
            radius = max(sep_deg(ra, dec, m['ra'], m['dec']) * 60.0 + m['size'] / 2
                         for m in members)
            groups.append({
                'id': r['ID'], 'name': r['Name'], 'type': T_GROUP,
                'ra': ra, 'dec': dec,
                'mag': combined_mag([m['mag'] for m in members]),
                'size': 2 * radius,
            })
            for m in members:
                member_of.setdefault(m['name'], r['Name'])
    return groups, member_of

# ---------------------------------------------------------------------------
def main():
    ngc_rows, lookup = load_openngc()
    overrides = load_overrides()
    groups, member_of = load_groups(lookup)

    # Cross-check hand-entered override positions against OpenNGC.
    for o in overrides:
        ref = next((lookup[k] for k in o['keys'] if k in lookup), None)
        if ref:
            s = sep_deg(o['ra'], o['dec'], ref['ra'], ref['dec'])
            if s > max(POS_TOLERANCE_DEG, o['size'] / 60.0 / 4):
                print("  WARNING: %s position differs from OpenNGC %s by %.2f deg"
                      % (o['id'], ref['name'], s))

    taken_keys  = {k for o in overrides for k in o['keys']}
    taken_names = {o['name'].lower() for o in overrides} | \
                  {g['name'].lower() for g in groups}

    objects = [dict(o) for o in overrides] + groups
    skipped_dim = 0
    for r in ngc_rows:
        t = OPENNGC_TYPES.get(r['type'])
        if t is None:
            continue
        if any(k in taken_keys for k in r['xref']):
            continue                                  # covered by an override
        if r['common'] and r['common'].lower() in taken_names:
            continue                                  # same object, other id
        if t == T_GALAXY and r['mag'] > GALAXY_MAG_LIMIT:
            skipped_dim += 1
            continue
        group = member_of.get(r['name'])
        name = (r['common'] or (group and 'in ' + group) or
                alt_identifier(r['idents']))
        objects.append({
            'id': display_id(r['name']), 'name': name, 'type': t,
            'ra': r['ra'], 'dec': r['dec'], 'mag': r['mag'], 'size': r['size'],
        })

    # String pool (deduplicated) + fixed-size records.
    pool, offsets = bytearray(b'\0'), {'': 0}
    def intern(s):
        s = ascii_text(s)
        if s not in offsets:
            offsets[s] = len(pool)
            pool.extend(s.encode('ascii') + b'\0')
        return offsets[s]

    recs = bytearray()
    for o in objects:
        recs += struct.pack(REC_FMT, o['ra'], o['dec'], o['mag'], o['size'],
                            intern(o['id']), intern(o['name']), o['type'])

    with open(DAT_PATH, 'wb') as f:
        f.write(struct.pack(HDR_FMT, b'NSKY', DAT_VERSION, len(objects), len(pool)))
        f.write(recs)
        f.write(pool)

    counts = [0] * 6
    for o in objects:
        counts[o['type']] += 1
    print("Wrote %s: %d objects (%d galaxies, %d nebulae, %d globulars, "
          "%d open clusters, %d other, %d groups), %d bytes."
          % (DAT_PATH, len(objects), *counts, 16 + len(recs) + len(pool)))
    print("  (%d galaxies fainter than mag %.1f left out)"
          % (skipped_dim, GALAXY_MAG_LIMIT))

if __name__ == '__main__':
    main()
