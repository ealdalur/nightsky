// nightsky_types.h -- on-disk layout of nightsky.dat, written by gencat.py.
//
//   DatHeader | CatRec[count] | string pool (strBytes, NUL-terminated strings)
//
// The whole file is read into memory in one go and these structs are used in
// place, so the layout must match gencat.py's HDR_FMT / REC_FMT exactly.
#pragma once
#include <cstdint>

#pragma pack(push,1)
struct DatHeader {
    char     magic[4];    // "NSKY"
    uint32_t version;     // kDatVersion
    uint32_t count;       // number of CatRec records
    uint32_t strBytes;    // size of the trailing string pool
};
#pragma pack(pop)

struct CatRec {
    float    ra_hours;    // right ascension, decimal hours (J2000)
    float    dec_deg;     // declination, decimal degrees (J2000)
    float    mag;         // apparent magnitude (V, else B); >= 90 means unknown
    float    size_arcmin; // largest angular dimension, arcmin; 0 means unknown
    uint32_t id;          // string-pool offset of the designation, e.g. "NGC7331"
    uint32_t name;        // string-pool offset of the display name ("" if none)
    uint32_t type;        // 0 Galaxy, 1 Nebula, 2 Globular, 3 OpenCluster,
                          // 4 Other, 5 Galaxy group
};

static_assert(sizeof(DatHeader) == 16, "DatHeader must match gencat.py HDR_FMT");
static_assert(sizeof(CatRec)    == 28, "CatRec must match gencat.py REC_FMT");
