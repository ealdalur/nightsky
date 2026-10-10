// nightsky_types.h -- POD record shared between the catalog generator's output
// (catalog_data.h) and the main program. Kept tiny and fixed-size so the whole
// catalog can be blitted to/from the binary cache file in a single read/write.
#pragma once

// Fixed-width record. Sizes chosen to hold every Messier designation/name with
// room to spare; short strings are null-terminated within the buffer.
struct CatEntry {
    char  id[8];      // primary designation, e.g. "M31"
    char  ngc[12];    // NGC/IC designation, e.g. "NGC224" ("" if none)
    char  name[40];   // common name ("" if none)
    int   type;       // 0 Galaxy, 1 Nebula, 2 Globular, 3 OpenCluster, 4 Other
    float ra_hours;   // right ascension, decimal hours (J2000)
    float dec_deg;    // declination, decimal degrees (J2000)
    float mag;        // apparent visual magnitude
    float size_arcmin;// largest angular dimension, arcminutes
};
