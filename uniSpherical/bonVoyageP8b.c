/* bonVoyageP8b: produce the itinerary report assuming first command line
   argument is a path/file name of a .ptb (Us8 coordinate format) binary
   file array of coordinates sorted in the itinerary order. The report is
   printed to the stdout fille.

   For instance:

   bonVoyageP8b w1904711.p8b
 */

#define PGM_DSCR "Report itinerary of .p8b (Us8 binary format) file"
#define PGM_LAST_EDIT_DATE "2025.085"
#include <stdio.h>

#include <nemo.h>
#include "../scullions/scullions.h" /* include after nemo.h has been included */

#define LINE_MAX      256
#define METERS2NM       0.0005399568
double arcLegLength(nemoPtUs8, nemoPtUs8);
double geodesicLegLength(nemoPtUs8, nemoPtUs8);

static const char *progName;    /* for error logging by this source file only */
/* ========================================================================== */
int main (int argc,
          const char *argv[],
          const char *envr[]) {
   int n;
   FILE *inFpCoords;               /* input binary file, location coordinates */
   nemoPtUs8 startUs8, locUs8, xLocUs8;    /* first, current, previous as Us8 */
   double arc, arcTotal, arcMin, arcMax, arcStartEnd;   /* on spherical Earth */
   double gds, gdsTotal, gdsMin, gdsMax, gdsStartEnd; /* on ellipsoidal Earth */
/* -------------------------------------------------------------------------- */
   progName = strrchr(argv[0], '/');                                 /* POSIX */
   if (progName == NULL) progName = strrchr(argv[0], '\\');         /* MS Win */
   if (progName == NULL) progName = argv[0];                      /* neither? */
   else progName += 1;                        /* strip leading path separator */
   fprintf(stderr, "\n[%s]: %s\nsource as of: %s, Nemo Library: %.3f\n",
                   progName, PGM_DSCR, PGM_LAST_EDIT_DATE, NEMO_LIBRARY_DATE);

   if (argc < 2) errorExit(progName, __LINE__,
      "command-line arguments: w1904711.ptb\n");

   inFpCoords = fopen(argv[1], "rb");                      /* Open input file */
   if (inFpCoords == NULL) errorExit(progName, __LINE__,
                                "Can't open [%s] for reading\n", argv[1]);

   arcTotal = arcMax = 0.0;
   gdsTotal = gdsMax = 0.0;
   arcMin = gdsMin = NEMO_DOUBLE_HUGE;
   n = 0;
/* first leg starting coordinate */
   fread(&startUs8, sizeof(nemoPtUs8), 1, inFpCoords);
   xLocUs8 = startUs8;
   while (fread(&locUs8, sizeof(nemoPtUs8), 1, inFpCoords)) { /* leg ending.. */
      arc = arcLegLength(xLocUs8, locUs8);
      if (arc < arcMin) arcMin = arc;
      if (arc > arcMax) arcMax = arc;
      arcTotal += arc;
      gds = geodesicLegLength(xLocUs8, locUs8);
      if (gds < gdsMin) gdsMin = gds;
      if (gds > gdsMax) gdsMax = gds;
      gdsTotal += gds;
      xLocUs8 = locUs8;
      n++;
      }
   fclose(inFpCoords);

/* let's hope the peddler does not end up exactly at the antipodes... */
   arcStartEnd = arcLegLength(startUs8, xLocUs8);
   gdsStartEnd = geodesicLegLength(startUs8, xLocUs8);

   printf("Itinerary from: %s, legs: %d\n", argv[1], n);
   printf("\"Open\" itinerary distances (meters, nautical miles):\n");
   printf("Spherical Earth (radius %10.3f meters):\n", NEMO_EARTH_RADIUS);
   printf("   minumum leg: %18.3f, %12.3f\n", arcMin, METERS2NM * arcMin);
   printf("   maximum leg: %18.3f, %12.3f\n", arcMax, METERS2NM * arcMax);
   printf("   total:       %18.3f, %12.3f\n", arcTotal, METERS2NM * arcTotal);
   printf("  (return home: %18.3f, %12.3f)\n", arcStartEnd, METERS2NM * arcStartEnd);
   printf("WGS84 Ellipsoid Earth:\n");
   printf("   minumum leg: %18.3f, %12.3f\n", gdsMin, METERS2NM * gdsMin);
   printf("   maximum leg: %18.3f, %12.3f\n", gdsMax, METERS2NM * gdsMax);
   printf("   total:       %18.3f, %12.3f\n", gdsTotal, METERS2NM * gdsTotal);
   printf("  (return home: %18.3f, %12.3f)\n", gdsStartEnd, METERS2NM * gdsStartEnd);
   printf("(For \"circular\" itinerary, add return leg to total. Et Bon Voyage!)\n");

   return(0);
   }
/* ========================================================================== */
double arcLegLength(nemoPtUs8 locA, nemoPtUs8 locB) {
   nemoPtNcs ptNcsA, ptNcsB;
   nemo_Us8ToNcs(locA, &ptNcsA);
   nemo_Us8ToNcs(locB, &ptNcsB);
   return(NEMO_EARTH_RADIUS * nemo_ArcV3(ptNcsA.dc, ptNcsB.dc));
   }
/* ========================================================================== */
double geodesicLegLength(nemoPtUs8 locA, nemoPtUs8 locB) {
   int nIter;
   double gds;
   nemoPtNcs ptNcsA, ptNcsB;
   nemoPtEnr ptEnrA, ptEnrB;
   nemo_Us8ToNcs(locA, &ptNcsA);
   nemo_Us8ToNcs(locB, &ptNcsB);
   nemo_NcsToEnr(nemo_ElrWgs84(), &ptNcsA, &ptEnrA);
   nemo_NcsToEnr(nemo_ElrWgs84(), &ptNcsB, &ptEnrB);
   gds = nemo_GeodesicSzpila(nemo_ElrWgs84(), &ptEnrA, &ptEnrB, &nIter);
   return(gds);
   }
/* ========================================================================== */
#include "../scullions/errorExit.c"
#include "../scullions/nemoStrings.c"
/* ========================================================================== */
