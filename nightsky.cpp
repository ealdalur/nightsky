// nightsky.cpp -- standalone terminal star/deep-sky observing planner.
//
// Ranks the most rewarding celestial objects for a given night and location:
// the Moon, the planets, and ~8,600 deep-sky objects (Messier, NGC/IC,
// Caldwell and named galaxy groups). Ranking blends how long the object clears
// a useful altitude during the dark hours, its apparent brightness, and how
// favourably it transits the meridian relative to the middle of the night.
//
// No third-party libraries. The deep-sky catalogue lives only in nightsky.dat,
// a compact binary file generated at build time by gencat.py and loaded with a
// single read; planet and Moon positions are computed.
//
// Build: see build.bat (MSVC) or build.sh.  Usage: see printUsage() below.

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <chrono>
#include <string>
#include <vector>
#include <algorithm>

// Device location uses the OS-native Windows Location Platform via C++/WinRT
// (headers + windowsapp.lib ship with the Windows SDK -- no third-party libs).
#ifdef _WIN32
  #ifndef NOMINMAX
    #define NOMINMAX
  #endif
  #ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN
  #endif
  #include <winrt/Windows.Foundation.h>
  #include <winrt/Windows.Devices.Geolocation.h>
#endif

#include "nightsky_types.h"

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------
static const double PI   = 3.14159265358979323846;
static const double DEG  = PI / 180.0;
static const double RAD  = 180.0 / PI;

static inline double d2r(double d){ return d * DEG; }
static inline double r2d(double r){ return r * RAD; }

// Reduce an angle in degrees to the range [0,360).
static double rev(double x){
    x = fmod(x, 360.0);
    if (x < 0) x += 360.0;
    return x;
}
static double sind(double d){ return sin(d2r(d)); }
static double cosd(double d){ return cos(d2r(d)); }
static double tand(double d){ return tan(d2r(d)); }
static double asind(double x){ return r2d(asin(x)); }
static double atan2d(double y,double x){ return r2d(atan2(y,x)); }

// Type labels shared with the CatRec.type codes emitted by gencat.py.
static const char* kTypeName[] = { "Galaxy", "Nebula", "Globular", "OpenClus", "Other", "GalGroup" };
static const char* typeName(int t){ return (t>=0 && t<6) ? kTypeName[t] : "?"; }

// ---------------------------------------------------------------------------
// Julian date / calendar conversions (all internal times are UT in JD)
// ---------------------------------------------------------------------------

// Gregorian calendar date + UT hours -> Julian Date.
static double calToJD(int y, int m, int d, double utHours){
    if (m <= 2){ y -= 1; m += 12; }
    int A = y / 100;
    int B = 2 - A + A / 4;
    double jd = floor(365.25 * (y + 4716))
              + floor(30.6001 * (m + 1))
              + d + B - 1524.5;
    return jd + utHours / 24.0;
}

// JD (already shifted to the desired zone by the caller) -> broken-down civil
// date/time. Returns whole H:M via out params and the calendar Y/M/D.
static void jdToCal(double jd, int& Y,int& M,int& D,int& hh,int& mm){
    double z = floor(jd + 0.5);
    double f = (jd + 0.5) - z;
    double A = z;
    if (z >= 2299161){
        double alpha = floor((z - 1867216.25) / 36524.25);
        A = z + 1 + alpha - floor(alpha / 4);
    }
    double B = A + 1524;
    double C = floor((B - 122.1) / 365.25);
    double Dd = floor(365.25 * C);
    double E = floor((B - Dd) / 30.6001);
    double day = B - Dd - floor(30.6001 * E) + f;
    D = (int)floor(day);
    M = (E < 14) ? (int)E - 1 : (int)E - 13;
    Y = (M > 2) ? (int)C - 4716 : (int)C - 4715;
    double dayFrac = day - D;
    double hours = dayFrac * 24.0;
    hh = (int)floor(hours);
    mm = (int)floor((hours - hh) * 60.0 + 0.5);
    if (mm == 60){ mm = 0; hh++; }
    if (hh == 24){ hh = 0; D++; }   // rollover is cosmetic for our display use
}

// Greenwich Mean Sidereal Time in degrees for a given UT JD.
static double gmstDeg(double jd){
    double T = (jd - 2451545.0) / 36525.0;
    double g = 280.46061837
             + 360.98564736629 * (jd - 2451545.0)
             + 0.000387933 * T * T
             - T * T * T / 38710000.0;
    return rev(g);
}

// ---------------------------------------------------------------------------
// Equatorial coordinates -> local altitude (degrees)
// ---------------------------------------------------------------------------
static double altitude(double raHours, double decDeg,
                       double jd, double latDeg, double lonDeg){
    double lstDeg = rev(gmstDeg(jd) + lonDeg);
    double haDeg  = rev(lstDeg - raHours * 15.0);
    double sinAlt = sind(latDeg) * sind(decDeg)
                  + cosd(latDeg) * cosd(decDeg) * cosd(haDeg);
    if (sinAlt >  1) sinAlt =  1;
    if (sinAlt < -1) sinAlt = -1;
    return asind(sinAlt);
}

// ---------------------------------------------------------------------------
// Solar-system positions
//
// Sun: Meeus low-precision. Moon & planets: Schlyter's method with the main
// periodic perturbations for the Moon, Jupiter, Saturn and Uranus. Accuracy is
// ~1-2 arcmin (Moon) to <0.05 deg (planets) -- far better than needed to rank
// rise/set/transit to the minute and estimate brightness.
// ---------------------------------------------------------------------------
struct EqPos { double ra; double dec; };          // ra in hours, dec in deg

// Solve Kepler's equation (all in degrees) -> eccentric anomaly.
static double kepler(double M, double e){
    M = rev(M);
    double E = M + r2d(e) * sind(M) * (1.0 + e * cosd(M));
    for (int i = 0; i < 8; ++i){
        double dE = (E - r2d(e) * sind(E) - M) / (1.0 - e * cosd(E));
        E -= dE;
        if (fabs(dE) < 1e-8) break;
    }
    return E;
}

// Obliquity of the ecliptic for the given day-number d (Schlyter).
static double oblEcl(double d){ return 23.4393 - 3.563e-7 * d; }

// Sun geocentric: returns ecliptic rectangular (xs,ys in AU) plus RA/Dec and
// its ecliptic longitude / mean anomaly (needed by the planet perturbations).
struct SunResult { double xs, ys; double ra, dec; double lon; double M; double dist; };
static SunResult computeSun(double d){
    double w = 282.9404 + 4.70935e-5 * d;
    double e = 0.016709 - 1.151e-9 * d;
    double M = rev(356.0470 + 0.9856002585 * d);
    double E = kepler(M, e);
    double xv = cosd(E) - e;
    double yv = sqrt(1 - e*e) * sind(E);
    double v  = atan2d(yv, xv);
    double r  = sqrt(xv*xv + yv*yv);
    double lon = rev(v + w);
    double xs = r * cosd(lon);
    double ys = r * sind(lon);
    double ecl = oblEcl(d);
    double xe = xs;
    double ye = ys * cosd(ecl);
    double ze = ys * sind(ecl);
    SunResult s;
    s.xs = xs; s.ys = ys; s.lon = lon; s.M = M; s.dist = r;
    s.ra  = rev(atan2d(ye, xe)) / 15.0;
    s.dec = atan2d(ze, sqrt(xe*xe + ye*ye));
    return s;
}
static EqPos sunEq(double jd){
    SunResult s = computeSun(jd - 2451543.5);
    return { s.ra, s.dec };
}

// Moon geocentric equatorial. Also outputs distance (Earth radii) and the
// Sun-Moon elongation (deg) for phase/brightness.
static EqPos moonEq(double jd, double* distER = nullptr, double* elong = nullptr){
    double d = jd - 2451543.5;
    double N = rev(125.1228 - 0.0529538083 * d);
    double i = 5.1454;
    double w = rev(318.0634 + 0.1643573223 * d);
    double a = 60.2666;
    double e = 0.054900;
    double M = rev(115.3654 + 13.0649929509 * d);
    double E = kepler(M, e);
    double xv = a * (cosd(E) - e);
    double yv = a * (sqrt(1 - e*e) * sind(E));
    double v  = atan2d(yv, xv);
    double r  = sqrt(xv*xv + yv*yv);
    // position in ecliptic coordinates
    double xh = r * (cosd(N)*cosd(v+w) - sind(N)*sind(v+w)*cosd(i));
    double yh = r * (sind(N)*cosd(v+w) + cosd(N)*sind(v+w)*cosd(i));
    double zh = r * (sind(v+w) * sind(i));
    double lon = atan2d(yh, xh);
    double lat = atan2d(zh, sqrt(xh*xh + yh*yh));
    // perturbations
    SunResult s = computeSun(d);
    double Ms = s.M;
    double Ls = rev(s.M + 282.9404 + 4.70935e-5 * d);   // sun mean longitude
    double Lm = rev(N + w + M);                          // moon mean longitude
    double D  = rev(Lm - Ls);                            // mean elongation
    double F  = rev(Lm - N);                             // argument of latitude
    double Mm = M;
    lon += -1.274 * sind(Mm - 2*D)
         +  0.658 * sind(2*D)
         -  0.186 * sind(Ms)
         -  0.059 * sind(2*Mm - 2*D)
         -  0.057 * sind(Mm - 2*D + Ms)
         +  0.053 * sind(Mm + 2*D)
         +  0.046 * sind(2*D - Ms)
         +  0.041 * sind(Mm - Ms)
         -  0.035 * sind(D)
         -  0.031 * sind(Mm + Ms)
         -  0.015 * sind(2*F - 2*D)
         +  0.011 * sind(Mm - 4*D);
    lat += -0.173 * sind(F - 2*D)
         -  0.055 * sind(Mm - F - 2*D)
         -  0.046 * sind(Mm + F - 2*D)
         +  0.033 * sind(F + 2*D)
         +  0.017 * sind(2*Mm + F);
    r   += -0.58 * cosd(Mm - 2*D)
         -  0.46 * cosd(2*D);
    // ecliptic -> equatorial
    double ecl = oblEcl(d);
    double xg = r * cosd(lon) * cosd(lat);
    double yg = r * sind(lon) * cosd(lat);
    double zg = r * sind(lat);
    double xe = xg;
    double ye = yg * cosd(ecl) - zg * sind(ecl);
    double ze = yg * sind(ecl) + zg * cosd(ecl);
    EqPos p;
    p.ra  = rev(atan2d(ye, xe)) / 15.0;
    p.dec = atan2d(ze, sqrt(xe*xe + ye*ye));
    if (distER) *distER = r;
    if (elong){
        // geocentric elongation from the Sun
        double cosElong = cosd(s.lon - lon) * cosd(lat);
        if (cosElong >  1) cosElong = 1;
        if (cosElong < -1) cosElong = -1;
        *elong = r2d(acos(cosElong));
    }
    return p;
}

// Planet indices
enum { P_MERCURY=0, P_VENUS, P_MARS, P_JUPITER, P_SATURN, P_URANUS, P_NEPTUNE, P_COUNT };

struct PlanetElem { double N0,Nd, i0,id_, w0,wd, a0,ad, e0,ed, M0,Md; };
// Orbital elements (Schlyter): value = base + rate*d
static const PlanetElem kPlanet[P_COUNT] = {
  // Mercury
  { 48.3313, 3.24587e-5, 7.0047, 5.00e-8, 29.1241, 1.01444e-5, 0.387098, 0.0,      0.205635, 5.59e-10, 168.6562, 4.0923344368 },
  // Venus
  { 76.6799, 2.46590e-5, 3.3946, 2.75e-8, 54.8910, 1.38374e-5, 0.723330, 0.0,      0.006773,-1.302e-9, 48.0052,  1.6021302244 },
  // Mars
  { 49.5574, 2.11081e-5, 1.8497,-1.78e-8, 286.5016,2.92961e-5, 1.523688, 0.0,      0.093405, 2.516e-9, 18.6021,  0.5240207766 },
  // Jupiter
  { 100.4542,2.76854e-5, 1.3030,-1.557e-7,273.8777,1.64505e-5, 5.20256,  0.0,      0.048498, 4.469e-9, 19.8950,  0.0830853001 },
  // Saturn
  { 113.6634,2.38980e-5, 2.4886,-1.081e-7,339.3939,2.97661e-5, 9.55475,  0.0,      0.055546,-9.499e-9, 316.9670, 0.0334442282 },
  // Uranus
  { 74.0005, 1.3978e-5,  0.7733, 1.9e-8,  96.6612, 3.0565e-5,  19.18171,-1.55e-8,  0.047318, 7.45e-9,  142.5905, 0.011725806 },
  // Neptune
  { 131.7806,3.0173e-5,  1.7700,-2.55e-7, 272.8461,-6.027e-6,  30.05826, 3.313e-8, 0.008606, 2.15e-9,  260.2471, 0.005995147 }
};

// Full geocentric equatorial position of a planet, plus heliocentric (r) and
// geocentric (R) distances for the magnitude/size model.
static EqPos planetEq(int idx, double jd, double* rOut=nullptr, double* ROut=nullptr,
                      double* fvOut=nullptr){
    double d = jd - 2451543.5;
    const PlanetElem& p = kPlanet[idx];
    double N = rev(p.N0 + p.Nd * d);
    double i = p.i0 + p.id_ * d;
    double w = rev(p.w0 + p.wd * d);
    double a = p.a0 + p.ad * d;
    double e = p.e0 + p.ed * d;
    double M = rev(p.M0 + p.Md * d);
    double E = kepler(M, e);
    double xv = a * (cosd(E) - e);
    double yv = a * (sqrt(1 - e*e) * sind(E));
    double v  = atan2d(yv, xv);
    double r  = sqrt(xv*xv + yv*yv);
    double xh = r * (cosd(N)*cosd(v+w) - sind(N)*sind(v+w)*cosd(i));
    double yh = r * (sind(N)*cosd(v+w) + cosd(N)*sind(v+w)*cosd(i));
    double zh = r * (sind(v+w) * sind(i));
    double lon = atan2d(yh, xh);
    double lat = atan2d(zh, sqrt(xh*xh + yh*yh));
    // major mutual perturbations
    double Mj = rev(kPlanet[P_JUPITER].M0 + kPlanet[P_JUPITER].Md * d);
    double Msa= rev(kPlanet[P_SATURN].M0  + kPlanet[P_SATURN].Md  * d);
    double Mu = rev(kPlanet[P_URANUS].M0  + kPlanet[P_URANUS].Md  * d);
    if (idx == P_JUPITER){
        lon += -0.332 * sind(2*Mj - 5*Msa - 67.6)
             -  0.056 * sind(2*Mj - 2*Msa + 21)
             +  0.042 * sind(3*Mj - 5*Msa + 21)
             -  0.036 * sind(Mj - 2*Msa)
             +  0.022 * cosd(Mj - Msa)
             +  0.023 * sind(2*Mj - 3*Msa + 52)
             -  0.016 * sind(Mj - 5*Msa - 69);
    } else if (idx == P_SATURN){
        lon += +0.812 * sind(2*Mj - 5*Msa - 67.6)
             -  0.229 * cosd(2*Mj - 4*Msa - 2)
             +  0.119 * sind(Mj - 2*Msa - 3)
             +  0.046 * sind(2*Mj - 6*Msa - 69)
             +  0.014 * sind(Mj - 3*Msa + 32);
        lat += -0.020 * cosd(2*Mj - 4*Msa - 2)
             +  0.018 * sind(2*Mj - 6*Msa - 49);
    } else if (idx == P_URANUS){
        lon += +0.040 * sind(Msa - 2*Mu + 6)
             +  0.035 * sind(Msa - 3*Mu + 33)
             -  0.015 * sind(Mj - Mu + 20);
    }
    // heliocentric ecliptic rectangular (perturbed)
    double xhp = r * cosd(lon) * cosd(lat);
    double yhp = r * sind(lon) * cosd(lat);
    double zhp = r * sind(lat);
    // add Sun's geocentric rectangular to get planet geocentric ecliptic
    SunResult s = computeSun(d);
    double xg = xhp + s.xs;
    double yg = yhp + s.ys;
    double zg = zhp;
    double ecl = oblEcl(d);
    double xe = xg;
    double ye = yg * cosd(ecl) - zg * sind(ecl);
    double ze = yg * sind(ecl) + zg * cosd(ecl);
    double R = sqrt(xg*xg + yg*yg + zg*zg);   // geocentric distance
    EqPos out;
    out.ra  = rev(atan2d(ye, xe)) / 15.0;
    out.dec = atan2d(ze, sqrt(xe*xe + ye*ye));
    if (rOut) *rOut = r;
    if (ROut) *ROut = R;
    if (fvOut){
        // phase angle FV: cos = (r^2 + R^2 - s^2)/(2rR), s = sun distance
        double sdist = s.dist;
        double cf = (r*r + R*R - sdist*sdist) / (2*r*R);
        if (cf >  1) cf = 1;
        if (cf < -1) cf = -1;
        *fvOut = r2d(acos(cf));
    }
    return out;
}

// ---------------------------------------------------------------------------
// Target model -- unifies fixed catalogue objects with the moving bodies.
// ---------------------------------------------------------------------------
enum Kind { K_FIXED, K_MOON, K_PLANET };

struct Target {
    const char* id;     // string literal or pointer into the loaded catalogue
    const char* name;
    int    type;        // CatRec type code
    Kind   kind;
    int    planet;      // valid when kind==K_PLANET
    // fixed-object data (also filled for moving bodies at reference time)
    double ra, dec;     // hours / deg
    double mag;
    double size;        // arcmin
};

static EqPos targetEq(const Target& t, double jd){
    switch (t.kind){
        case K_MOON:   return moonEq(jd);
        case K_PLANET: return planetEq(t.planet, jd);
        default:       return { t.ra, t.dec };
    }
}

// ---------------------------------------------------------------------------
// Event finding: transit within the noon-to-noon day of the night, the rise
// before and set after it, plus in-window visibility statistics for scoring.
//
// Every object is sampled on the same grid (2-minute steps), so the sidereal
// angle at each sample is computed once per night. The grid runs 12 h either
// side of the noon-to-noon day: an object is never more than 12 h from its
// transit when it rises or sets, so every rise/set is found, and "no rise/set"
// genuinely means the object never crosses the horizon (circumpolar).
//
// Objects are evaluated in sin(altitude) space: for a fixed object that is
//   sin(alt) = sinLat*sinDec + cosLat*cosDec*cos(LST - RA)
// with cos(LST - RA) expanded so each sample costs a few multiply-adds and no
// trig. asin() is only taken for the handful of values that are reported.
// ---------------------------------------------------------------------------
struct SkyGrid {
    double startJD, step;
    int    N;
    int    dayK0, dayK1;     // sample range of the noon-to-noon day
    double sinLat, cosLat;
    std::vector<double> jd, cosLst, sinLst;
};

static SkyGrid makeGrid(double noonJD, double lat, double lon){
    SkyGrid g;
    g.step    = 2.0 / (24.0 * 60.0);
    int perDay = (int)(1.0 / g.step + 0.5);
    g.startJD = noonJD - 0.5;
    g.N       = 2 * perDay + 1;
    g.dayK0   = perDay / 2;
    g.dayK1   = g.dayK0 + perDay;
    g.sinLat  = sind(lat);
    g.cosLat  = cosd(lat);
    g.jd.resize(g.N); g.cosLst.resize(g.N); g.sinLst.resize(g.N);
    for (int k = 0; k < g.N; ++k){
        double jd  = g.startJD + k * g.step;
        double lst = gmstDeg(jd) + lon;
        g.jd[k] = jd; g.cosLst[k] = cosd(lst); g.sinLst[k] = sind(lst);
    }
    return g;
}

// sin(altitude) of the target at every grid sample.
static void sinAltSeries(const Target& t, const SkyGrid& g, double* s){
    if (t.kind == K_FIXED){
        double a  = g.sinLat * sind(t.dec), b = g.cosLat * cosd(t.dec);
        double cr = b * cosd(t.ra * 15.0), sr = b * sind(t.ra * 15.0);
        for (int k = 0; k < g.N; ++k)
            s[k] = a + g.cosLst[k] * cr + g.sinLst[k] * sr;
    } else {
        for (int k = 0; k < g.N; ++k){
            EqPos p = targetEq(t, g.jd[k]);
            double b = g.cosLat * cosd(p.dec);
            s[k] = g.sinLat * sind(p.dec)
                 + g.cosLst[k] * b * cosd(p.ra * 15.0)
                 + g.sinLst[k] * b * sind(p.ra * 15.0);
        }
    }
}

static double altFromSin(double s){ return asind(s > 1 ? 1 : (s < -1 ? -1 : s)); }

struct Events {
    bool   hasRise, hasSet;
    double riseJD, transitJD, setJD;
    double transitAlt;       // altitude at transit (deg)
    double visHours;         // hours above visThreshold within [winStart,winEnd]
    double peakJD;           // time of max altitude within the window
    double peakAlt;          // that max altitude
    bool   everUp;           // rose above visThreshold at all during the window
};

static Events findEvents(const SkyGrid& g, const double* s,
                         double winStart, double winEnd, double visThresh){
    const int N = g.N;
    Events ev{};
    ev.peakAlt = -999; ev.peakJD = winStart;

    // Transit = highest sample of the day, refined with a parabola through
    // the neighbouring altitudes.
    int bestI = g.dayK0;
    for (int k = g.dayK0 + 1; k <= g.dayK1; ++k) if (s[k] > s[bestI]) bestI = k;
    ev.transitJD  = g.jd[bestI];
    ev.transitAlt = altFromSin(s[bestI]);
    if (bestI > 0 && bestI < N-1){
        double y0 = altFromSin(s[bestI-1]), y1 = ev.transitAlt, y2 = altFromSin(s[bestI+1]);
        double denom = (y0 - 2*y1 + y2);
        if (fabs(denom) > 1e-9){
            double dx = 0.5 * (y0 - y2) / denom;    // in samples
            ev.transitJD  = g.startJD + (bestI + dx) * g.step;
            ev.transitAlt = y1 - 0.25 * (y0 - y2) * dx;
        }
    }
    // Rise: last upward crossing of the horizon before transit.
    for (int k = bestI; k > 0; --k){
        if (s[k] >= 0 && s[k-1] < 0){
            double f = -s[k-1] / (s[k] - s[k-1]);
            ev.riseJD = g.startJD + (k-1 + f) * g.step;
            ev.hasRise = true;
            break;
        }
    }
    // Set: first downward crossing of the horizon after transit.
    for (int k = bestI; k < N-1; ++k){
        if (s[k] >= 0 && s[k+1] < 0){
            double f = s[k] / (s[k] - s[k+1]);
            ev.setJD = g.startJD + (k + f) * g.step;
            ev.hasSet = true;
            break;
        }
    }
    // In-window visibility above the threshold altitude.
    int k0 = (int)ceil ((winStart - g.startJD) / g.step);
    int k1 = (int)floor((winEnd   - g.startJD) / g.step);
    if (k0 < 0) k0 = 0;
    if (k1 > N-1) k1 = N-1;
    double sThresh = sind(visThresh);
    int peakK = -1, upSamples = 0;
    for (int k = k0; k <= k1; ++k){
        if (s[k] >= sThresh) ++upSamples;
        if (peakK < 0 || s[k] > s[peakK]) peakK = k;
    }
    ev.everUp   = upSamples > 0;
    ev.visHours = upSamples * g.step * 24.0;
    if (peakK >= 0){ ev.peakJD = g.jd[peakK]; ev.peakAlt = altFromSin(s[peakK]); }
    return ev;
}

// Find the Sun crossing of altitude h0 (deg) between t0 and t1, going in the
// given direction (-1 = setting, +1 = rising). Returns JD or -1 if none.
static double sunCrossing(double t0, double t1, double lat, double lon,
                          double h0, int dir){
    const double step = 1.0 / (24.0 * 60.0);   // 1-min resolution
    double prevJD = t0;
    EqPos sp = sunEq(t0);
    double prev = altitude(sp.ra, sp.dec, t0, lat, lon) - h0;
    for (double jd = t0 + step; jd <= t1; jd += step){
        EqPos p = sunEq(jd);
        double cur = altitude(p.ra, p.dec, jd, lat, lon) - h0;
        bool cross = (dir < 0) ? (prev >= 0 && cur < 0) : (prev < 0 && cur >= 0);
        if (cross){
            double f = (0 - prev) / (cur - prev);
            return prevJD + f * (jd - prevJD);
        }
        prev = cur; prevJD = jd;
    }
    return -1;
}

// ---------------------------------------------------------------------------
// Brightness / size models for the moving bodies
// ---------------------------------------------------------------------------
static void moonMagSize(double jd, double& mag, double& sizeArcmin, double& illum){
    double distER = 60.0, elong = 90.0;
    moonEq(jd, &distER, &elong);
    double phaseAngle = 180.0 - elong;              // 0 = full
    mag = -12.7 + 0.026 * fabs(phaseAngle)
        + 4.0e-9 * pow(fabs(phaseAngle), 4.0);
    double sizeArcsec = 112489.0 / distER;          // 2*Rmoon/d in arcsec
    sizeArcmin = sizeArcsec / 60.0;
    illum = (1.0 + cosd(phaseAngle)) / 2.0;
}

static void planetMagSize(int idx, double jd, double& mag, double& sizeArcmin){
    double r=1, R=1, fv=0;
    planetEq(idx, jd, &r, &R, &fv);
    double logrR = 5.0 * log10(r * R);
    switch (idx){
        case P_MERCURY: mag = -0.36 + logrR + 0.027*fv + 2.2e-13*pow(fv,6); break;
        case P_VENUS:   mag = -4.34 + logrR + 0.013*fv + 4.2e-7*pow(fv,3);  break;
        case P_MARS:    mag = -1.51 + logrR + 0.016*fv;                     break;
        case P_JUPITER: mag = -9.25 + logrR + 0.014*fv;                     break;
        case P_SATURN:  mag = -9.00 + logrR + 0.044*fv;                     break;
        case P_URANUS:  mag = -7.15 + logrR + 0.001*fv;                     break;
        default:        mag = -6.90 + logrR + 0.001*fv;                     break;
    }
    // apparent equatorial diameter in arcsec at 1 AU / R(AU)
    static const double d0[P_COUNT] = {6.728,16.688,9.36,196.74,165.5,65.8,62.2};
    double arcsec = d0[idx] / R;
    sizeArcmin = arcsec / 60.0;
}

// ---------------------------------------------------------------------------
// Catalogue loading: read the binary data file produced at build time by
// gencat.py. The catalogue lives ONLY in this file (the program embeds no
// copy), so a missing/incompatible file is a hard error with a clear hint.
// The file is read with a single fread and its records and strings are used
// in place -- no per-record parsing or copying.
// ---------------------------------------------------------------------------
static const uint32_t kDatVersion = 3;
static const char*    kDatPath    = "nightsky.dat";
static const int      kTypeCount  = 6;

struct Catalog {
    std::vector<char> buf;                 // the whole of nightsky.dat
    const CatRec*     rec   = nullptr;
    uint32_t          count = 0;
    const char*       str   = nullptr;     // string pool
    const char* s(uint32_t off) const { return str + off; }
};

// Returns false on any failure, with a reason in `err`.
static bool loadCatalog(Catalog& cat, std::string& err){
    FILE* f = fopen(kDatPath, "rb");
    if (!f){ err = "cannot open " + std::string(kDatPath); return false; }
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (size < (long)sizeof(DatHeader)){
        fclose(f); err = std::string(kDatPath) + " is not a valid catalogue file"; return false;
    }
    cat.buf.resize((size_t)size);
    size_t got = fread(cat.buf.data(), 1, cat.buf.size(), f);
    fclose(f);
    if (got != cat.buf.size()){ err = std::string(kDatPath) + " could not be read"; return false; }

    DatHeader h;
    memcpy(&h, cat.buf.data(), sizeof(h));
    if (memcmp(h.magic, "NSKY", 4) != 0){
        err = std::string(kDatPath) + " is not a valid catalogue file"; return false;
    }
    if (h.version != kDatVersion){
        char b[96]; snprintf(b,sizeof(b),"%s version %u, expected %u",
                             kDatPath, h.version, kDatVersion);
        err = b; return false;
    }
    uint64_t expect = sizeof(DatHeader) + (uint64_t)h.count * sizeof(CatRec) + h.strBytes;
    if (h.count == 0 || h.strBytes == 0 || expect != (uint64_t)size ||
        cat.buf.back() != '\0'){
        err = std::string(kDatPath) + " is corrupt or truncated"; return false;
    }
    cat.rec   = reinterpret_cast<const CatRec*>(cat.buf.data() + sizeof(DatHeader));
    cat.count = h.count;
    cat.str   = cat.buf.data() + sizeof(DatHeader) + (size_t)h.count * sizeof(CatRec);
    for (uint32_t i = 0; i < cat.count; ++i){
        const CatRec& r = cat.rec[i];
        if (r.id >= h.strBytes || r.name >= h.strBytes || r.type >= (uint32_t)kTypeCount){
            err = std::string(kDatPath) + " is corrupt (bad record)"; return false;
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// Formatting helpers
// ---------------------------------------------------------------------------
// Format a UT JD as local HH:MM (zone offset already added by caller path).
static std::string fmtLocalHM(double jdUT, double tzHours){
    int Y,M,D,hh,mm;
    jdToCal(jdUT + tzHours/24.0, Y,M,D,hh,mm);
    char b[8]; snprintf(b, sizeof(b), "%02d:%02d", hh, mm);
    return std::string(b);
}
static std::string fmtMag(double mag){
    if (mag >= 90.0) return "--";                                           // not catalogued
    char b[16]; snprintf(b, sizeof(b), "%.1f", mag);
    return std::string(b);
}
static std::string fmtSize(double arcmin){
    char b[16];
    if (arcmin <= 0.0) return "--";                                         // not catalogued
    if (arcmin < 1.0) snprintf(b, sizeof(b), "%.0f\"", arcmin * 60.0);       // arcsec
    else              snprintf(b, sizeof(b), "%.1f'", arcmin);               // arcmin
    return std::string(b);
}

static void printUsage(){
    printf(
    "nightsky - top celestial objects to observe/photograph tonight\n\n"
    "Usage:\n"
    "  nightsky [count] [YYYYMMDD] [lat lon]\n\n"
    "  count      number of objects to list (default 10)\n"
    "  YYYYMMDD   evening date of the night to evaluate (the date just before\n"
    "             the midnight of that night). Also accepts YYYY/MM/DD.\n"
    "  lat lon    observer latitude and longitude in decimal degrees\n"
    "  --no-locate  skip device-location lookup (use nightsky.cfg / default)\n\n"
    "Examples:\n"
    "  nightsky                     tonight, current location, top 10\n"
    "  nightsky 20                  tonight, current location, top 20\n"
    "  nightsky 10 20270304         night of 2027-03-04, current location\n"
    "  nightsky 10 20270304 34.8633 -111.7510   at the given coordinates\n\n"
    "With no coordinates, nightsky uses the device's location via the Windows\n"
    "Location Platform (GPS/Wi-Fi). If Location Services are off or denied, it\n"
    "falls back to nightsky.cfg (a single line \"lat lon\") then a built-in default.\n");
}

// ---------------------------------------------------------------------------
// Location resolution
// ---------------------------------------------------------------------------
// Fallback when neither CLI coordinates nor device location nor a config file
// are available.
static const double kDefaultLat =  33.8493053;   // La Palma Central Park,
static const double kDefaultLon = -118.0379991;   // La Palma, CA
static const char*  kCfgPath    = "nightsky.cfg";

// The device's real position, resolved from the OS Windows Location Platform
// (fused GPS / Wi-Fi / cellular). Returns false -- with a human-readable reason
// in `note` -- if location services are off, permission is denied, or no fix is
// obtained within the timeout, so callers can fall back to cfg/default.
#ifdef _WIN32
static bool tryDeviceLocation(double& lat, double& lon, std::string& note){
    using namespace winrt::Windows::Devices::Geolocation;
    try {
        winrt::init_apartment();
        auto access = Geolocator::RequestAccessAsync().get();
        if (access != GeolocationAccessStatus::Allowed){
            note = "location access not allowed "
                   "(Settings > Privacy & security > Location)";
            return false;
        }
        Geolocator locator;
        locator.DesiredAccuracyInMeters(500u);   // city-scale is plenty; quicker fix
        // Accept a cached fix up to 5 min old; give up after 10 s so we never hang.
        auto pos = locator.GetGeopositionAsync(std::chrono::minutes(5),
                                               std::chrono::seconds(10)).get();
        auto p = pos.Coordinate().Point().Position();
        lat = p.Latitude;
        lon = p.Longitude;
        return true;
    } catch (winrt::hresult_error const& e){
        note = winrt::to_string(e.message());
        return false;
    } catch (...) {
        note = "location query failed";
        return false;
    }
}
#else
static bool tryDeviceLocation(double&, double&, std::string& note){
    note = "device location only supported on Windows";
    return false;
}
#endif

// A manually configured location: a single "lat lon" line in nightsky.cfg.
static bool readCfgLocation(double& lat, double& lon){
    FILE* f = fopen(kCfgPath, "r");
    if (!f) return false;
    double a, b;
    bool ok = (fscanf(f, "%lf %lf", &a, &b) == 2);
    fclose(f);
    if (ok){ lat = a; lon = b; }
    return ok;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
// Ranking weights (sum to 1). Brightness carries the most weight because, with
// thousands of catalogue objects, it is what separates showpieces from filler.
// Brightness scores 1 at mag -2 and falls linearly to 0 at kMagFaint.
static const double kWDur     = 0.25;   // share of the dark window above threshold
static const double kWBright  = 0.50;   // apparent magnitude
static const double kWMerid   = 0.25;   // peak altitude, timed near mid-night
static const double kMagFaint = 14.0;
struct Ranked { Target t; Events ev; double score; double dispMag, dispSize; };

int main(int argc, char** argv){
    int    count = 10;
    bool   dateGiven = false;
    int    yy=0, mo=0, dd=0;
    bool   locGiven = false;
    double lat = kDefaultLat, lon = kDefaultLon;
    const char* locSource = "built-in default";

    // ---- parse arguments -------------------------------------------------
    bool autoLocate = true;
    std::vector<std::string> args;
    for (int i = 1; i < argc; ++i){
        std::string a = argv[i];
        if (a == "-h" || a == "--help" || a == "-help" || a == "/?"){ printUsage(); return 0; }
        if (a == "--no-locate"){ autoLocate = false; continue; }
        args.push_back(a);
    }
    if (args.size() >= 1){
        count = atoi(args[0].c_str());
        if (count <= 0 || count > 100000){ fprintf(stderr, "Invalid count '%s'.\n\n", args[0].c_str()); printUsage(); return 1; }
    }
    if (args.size() >= 2){
        std::string digits;
        for (char c : args[1]) if (c >= '0' && c <= '9') digits.push_back(c);
        if (digits.size() != 8){ fprintf(stderr, "Date must be YYYYMMDD (got '%s').\n\n", args[1].c_str()); printUsage(); return 1; }
        yy = atoi(digits.substr(0,4).c_str());
        mo = atoi(digits.substr(4,2).c_str());
        dd = atoi(digits.substr(6,2).c_str());
        if (mo < 1 || mo > 12 || dd < 1 || dd > 31){ fprintf(stderr, "Invalid date '%s'.\n\n", args[1].c_str()); printUsage(); return 1; }
        dateGiven = true;
    }
    if (args.size() >= 4){
        lat = atof(args[2].c_str());
        lon = atof(args[3].c_str());
        if (lat < -90 || lat > 90 || lon < -180 || lon > 180){ fprintf(stderr, "Latitude must be -90..90 and longitude -180..180.\n\n"); printUsage(); return 1; }
        locGiven = true;
        locSource = "command line";
    } else if (args.size() == 3){
        fprintf(stderr, "Provide both latitude AND longitude, or neither.\n\n"); printUsage(); return 1;
    }

    // ---- resolve location ------------------------------------------------
    // Priority: command line > device location (Windows) > nightsky.cfg > default.
    if (!locGiven){
        double dlat, dlon; std::string note = "auto-locate disabled";
        if (autoLocate && tryDeviceLocation(dlat, dlon, note)){
            lat = dlat; lon = dlon; locSource = "device location";
        } else {
            double clat, clon;
            if (readCfgLocation(clat, clon)){ lat = clat; lon = clon; locSource = "nightsky.cfg"; }
            else fprintf(stderr, "nightsky: device location unavailable (%s); using %s.\n",
                         note.c_str(), locSource);
        }
    }

    // ---- resolve timezone (display offset) -------------------------------
    // The observing site's local clock zone is estimated from its longitude
    // (15 deg per hour). This keeps the reported times consistent with the
    // location being evaluated rather than with whatever zone this machine's
    // clock happens to be set to. Political zones / DST are intentionally
    // ignored; the offset is stated plainly in the header.
    time_t nowT = time(nullptr);
    double tzHours = floor(lon / 15.0 + 0.5);

    // ---- resolve the evening date ----------------------------------------
    if (!dateGiven){
        double nowJD = 2440587.5 + (double)nowT / 86400.0;
        int Y,M,D,hh,mm;
        jdToCal(nowJD + tzHours/24.0, Y,M,D,hh,mm);   // local calendar date now
        yy = Y; mo = M; dd = D;
    }

    // ---- load catalogue --------------------------------------------------
    std::string catErr;
    Catalog cat;
    if (!loadCatalog(cat, catErr)){
        fprintf(stderr,
            "nightsky: %s.\n"
            "The object catalogue file '%s' is required. Run build.bat (which\n"
            "runs gencat.py) to generate it, and keep it next to nightsky.exe.\n",
            catErr.c_str(), kDatPath);
        return 1;
    }

    // ---- night window: sunset(date) -> sunrise(date+1) -------------------
    double noonJD = calToJD(yy, mo, dd, 12.0 - tzHours);   // local noon of the evening, in UT
    double sunAlt = -0.833;   // standard sunset/sunrise altitude
    double sunset  = sunCrossing(noonJD, noonJD + 0.75, lat, lon, sunAlt, -1);
    double sunrise = sunCrossing(sunset > 0 ? sunset : noonJD + 0.25, noonJD + 1.0, lat, lon, sunAlt, +1);

    // Polar day/night fallbacks: if the sun never crosses, use civil clock night.
    bool noSunset = (sunset < 0), noSunrise = (sunrise < 0);
    if (noSunset)  sunset  = calToJD(yy, mo, dd, 20.0 - tzHours);         // 20:00 local
    if (noSunrise) sunrise = calToJD(yy, mo, dd, 6.0  - tzHours) + 1.0;   // 06:00 next day
    if (sunrise <= sunset) sunrise = sunset + 0.35;                       // safety

    double winStart = sunset, winEnd = sunrise;

    // For "tonight" (no explicit date), if we're already past sunset, start the
    // observing window at the current moment.
    double nowJD = 2440587.5 + (double)nowT / 86400.0;
    bool clampedNow = false;
    if (!dateGiven && nowJD > sunset && nowJD < sunrise){
        winStart = nowJD; clampedNow = true;
    }
    double winMid  = 0.5 * (winStart + winEnd);
    double winHrs  = (winEnd - winStart) * 24.0;
    double halfWin = 0.5 * (winEnd - winStart);

    // ---- build the target pool -------------------------------------------
    std::vector<Target> pool;
    pool.reserve(cat.count + 1 + P_COUNT);
    // Moon
    { Target t{}; t.id="Moon"; t.name="The Moon"; t.type=4; t.kind=K_MOON; t.planet=-1; pool.push_back(t); }
    // Planets
    { static const char* pid[P_COUNT]={"Mercury","Venus","Mars","Jupiter","Saturn","Uranus","Neptune"};
      for (int i=0;i<P_COUNT;++i){ Target t{}; t.id=pid[i]; t.name=pid[i]; t.type=4; t.kind=K_PLANET; t.planet=i; pool.push_back(t); } }
    // Catalogue objects
    for (uint32_t i = 0; i < cat.count; ++i){
        const CatRec& e = cat.rec[i];
        Target t{}; t.id=cat.s(e.id); t.name=cat.s(e.name); t.type=(int)e.type; t.kind=K_FIXED; t.planet=-1;
        t.ra=e.ra_hours; t.dec=e.dec_deg; t.mag=e.mag; t.size=e.size_arcmin;
        pool.push_back(t);
    }

    // ---- evaluate + score ------------------------------------------------
    SkyGrid grid = makeGrid(noonJD, lat, lon);
    std::vector<double> sinAlt(grid.N);
    std::vector<Ranked> ranked;
    for (const Target& t : pool){
        double visThresh = (t.kind==K_FIXED) ? 20.0 : 10.0;   // photographic vs solar-system
        // A fixed object culminates at 90 - |lat - dec|; skip any that can
        // never clear the threshold without sampling them.
        if (t.kind==K_FIXED && 90.0 - fabs(lat - t.dec) < visThresh) continue;
        sinAltSeries(t, grid, sinAlt.data());
        Events ev = findEvents(grid, sinAlt.data(), winStart, winEnd, visThresh);
        if (!ev.everUp) continue;                             // not observable this night

        // magnitude / size at mid-window for moving bodies
        double mag = t.mag, size = t.size, illum = 0;
        if (t.kind==K_MOON)   moonMagSize(winMid, mag, size, illum);
        if (t.kind==K_PLANET) planetMagSize(t.planet, winMid, mag, size);

        // scoring components (0..1)
        double durScore = ev.visHours / (winHrs > 0 ? winHrs : 1.0);
        if (durScore > 1) durScore = 1;
        double brightScore = (kMagFaint - mag) / (kMagFaint + 2.0);   // -2 ->1, kMagFaint ->0
        if (brightScore < 0) brightScore = 0; if (brightScore > 1) brightScore = 1;
        double altScore = ev.peakAlt > 0 ? sind(ev.peakAlt) : 0.0;
        double timeScore = 1.0 - fabs(ev.peakJD - winMid) / (halfWin > 0 ? halfWin : 1.0);
        if (timeScore < 0) timeScore = 0;
        double meridScore = altScore * (0.5 + 0.5 * timeScore);

        double score = kWDur*durScore + kWBright*brightScore + kWMerid*meridScore;

        Ranked r; r.t=t; r.ev=ev; r.score=score; r.dispMag=mag; r.dispSize=size;
        ranked.push_back(r);
    }
    std::sort(ranked.begin(), ranked.end(),
              [](const Ranked&a,const Ranked&b){ return a.score > b.score; });
    if ((int)ranked.size() > count) ranked.resize(count);

    // ---- moon phase for the header ---------------------------------------
    double mMag, mSize, mIllum; moonMagSize(winMid, mMag, mSize, mIllum);

    // ---- build the report (to both stdout and the text file) -------------
    static const char* kRuleEq =
        "==============================================================================================\n";
    static const char* kRuleDash =
        "----------------------------------------------------------------------------------------------\n";
    std::string out;
    char line[512];
    auto emit = [&](const char* s){ out += s; };

    emit(kRuleEq);
    snprintf(line,sizeof(line),
        " NIGHTSKY  -  top %d objects for the night of %04d-%02d-%02d (evening)\n",
        count, yy, mo, dd); emit(line);
    emit(kRuleEq);
    snprintf(line,sizeof(line),
        " Location : %+.4f, %+.4f  (%s)\n", lat, lon, locSource); emit(line);
    snprintf(line,sizeof(line),
        " Timezone : UTC%+05.1f    Night: %s -> %s   (%.1f h dark)\n",
        tzHours,
        fmtLocalHM(sunset, tzHours).c_str(), fmtLocalHM(sunrise, tzHours).c_str(),
        (sunrise - sunset) * 24.0); emit(line);
    if (clampedNow){
        snprintf(line,sizeof(line),
        " Window   : ranking from now (%s) until dawn\n",
        fmtLocalHM(winStart, tzHours).c_str()); emit(line);
    }
    snprintf(line,sizeof(line),
        " Moon     : %.0f%% illuminated, mag %.1f\n", mIllum*100.0, mMag); emit(line);
    if (noSunset || noSunrise){ emit(" Note     : sun does not cross horizon (polar) - used civil-night hours\n"); }
    snprintf(line,sizeof(line),
        " Times are local (UTC%+05.1f).  Rise/Set at 0deg horizon, transit=meridian.\n", tzHours); emit(line);
    emit(kRuleDash);
    snprintf(line,sizeof(line),
        " %-4s %-10s %-22s %-9s %6s %8s  %5s %5s %5s %6s\n",
        "#","ID","Name","Type","Mag","Size","Rise","Trans","Set","MaxAlt"); emit(line);
    emit(kRuleDash);

    int rank = 1;
    for (const Ranked& r : ranked){
        std::string nm = r.t.name; if (nm.size() > 22) nm = nm.substr(0,22);
        std::string rise = r.ev.hasRise ? fmtLocalHM(r.ev.riseJD, tzHours) : " up  ";
        std::string set  = r.ev.hasSet  ? fmtLocalHM(r.ev.setJD,  tzHours) : " up  ";
        std::string tr   = fmtLocalHM(r.ev.transitJD, tzHours);
        snprintf(line,sizeof(line),
            " %-4d %-10s %-22s %-9s %6s %8s  %5s %5s %5s %5.0f deg\n",
            rank++, r.t.id, nm.c_str(), typeName(r.t.type),
            fmtMag(r.dispMag).c_str(), fmtSize(r.dispSize).c_str(),
            rise.c_str(), tr.c_str(), set.c_str(),
            r.ev.peakAlt);
        emit(line);
    }
    emit(kRuleDash);
    snprintf(line,sizeof(line),
        " Ranking = %.2f*visibility-duration + %.2f*brightness + %.2f*meridian-favourability.\n",
        kWDur, kWBright, kWMerid); emit(line);
    emit(" Deep-sky counted as visible above 20deg; solar-system above 10deg.\n");
    snprintf(line,sizeof(line),
        " Objects: Moon, planets, and %u deep-sky objects (Messier, NGC/IC, Caldwell, groups).\n"
        " Deep-sky data from OpenNGC (CC BY-SA 4.0). Mag is V where known, else B. '--' = not catalogued.\n",
        cat.count); emit(line);
    emit(kRuleEq);

    // ---- output ----------------------------------------------------------
    fputs(out.c_str(), stdout);

    char fname[64]; snprintf(fname,sizeof(fname),"nightsky_%04d%02d%02d.txt", yy, mo, dd);
    FILE* fo = fopen(fname, "w");
    if (fo){ fputs(out.c_str(), fo); fclose(fo); printf("Wrote %s\n", fname); }
    else    { fprintf(stderr, "Warning: could not write %s\n", fname); }
    return 0;
}
