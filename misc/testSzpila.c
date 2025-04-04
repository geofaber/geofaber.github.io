/* A simple "smoke test" program, confirming the linking (or source inclusion)
   and execution of a Nemo library function. The Library function used here
   is Szpila's calculation (Vincenty formulae) of the lengh of geodesic with
   two Australian triangulation stations used as the example in the original
   publication of the formulae. The expected result is 54972.271 meters.

   Programmer: Hrvoje Lukatela, 2022.
 */

#define PGM_DSCR "Vincenty geodesics Nemo Librarty \"smoke-test\""
#define PGM_LAST_EDIT_DATE "2024.199"         /* format as from 'date +%Y.%j' */

#include <stdio.h>
#include <string.h>

#include <nemo.h>

/* ========================================================================== */
int main (int argc,
          const char *argv[],
          const char *envr[]) {

   int iterCount;
   nemoPtEll ptIn;                                           /* given φ and λ */
   nemoPtEnr ptA, ptB;                       /* two given points on ellipsoid */
   double geodesicLength;                  /* length of geodesic between them */
   const char *progName;        /* for error logging by this source file only */
/* -------------------------------------------------------------------------- */
   progName = strrchr(argv[0], '/');                                 /* POSIX */
   if (progName == NULL) progName = strrchr(argv[0], '\\');         /* MS Win */
   if (progName) progName += 1;             /* skip found last path separator */
   if (progName == NULL) progName = argv[0];     /* neither? Just use argv[0] */
   fprintf(stderr, "\n[%s]: %s\nsource as of: %s, Nemo Library: %.3f\n",
                   progName, PGM_DSCR, PGM_LAST_EDIT_DATE, NEMO_LIBRARY_DATE);

/* Coordinates of Flinder's Peak: */
   ptIn.a[NEMO_LAT] = NEMO_DEG2RAD * -37.951033417;
   ptIn.a[NEMO_LNG] = NEMO_DEG2RAD * 144.424867889;
   nemo_LatLongToDcos3(ptIn.a, ptA.dc);
   printf("Flinder's Peak: %14.9f, %15.9f\n", NEMO_RAD2DEG * ptIn.a[NEMO_LAT],
                                              NEMO_RAD2DEG * ptIn.a[NEMO_LNG]);
/* Coordinates of Buninyong: */
   ptIn.a[NEMO_LAT] = NEMO_DEG2RAD * -37.652821139;
   ptIn.a[NEMO_LNG] = NEMO_DEG2RAD * 143.926495528;
   nemo_LatLongToDcos3(ptIn.a, ptB.dc);
   printf("Buninyong:      %14.9f, %15.9f\n", NEMO_RAD2DEG * ptIn.a[NEMO_LAT],
                                              NEMO_RAD2DEG * ptIn.a[NEMO_LNG]);

   geodesicLength = nemo_GeodesicSzpila(nemo_ElrWgs84(), &ptA, &ptB, &iterCount);

   if (geodesicLength == NEMO_DOUBLE_UNDEF) printf(
      "Unexpected error; nemo_GeodesicSzpila() failed to converge?\n");
   else printf(
      "geodesic length: %.3f, %d iterations\n", geodesicLength, iterCount);

   return(0);
   }
/* ========================================================================== */
