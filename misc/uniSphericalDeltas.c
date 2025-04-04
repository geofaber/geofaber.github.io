/* Monte-Carlo test of UniSpherical coordinate encoding Δ distribution.
   The program requires no files and takes one command line argument, the
   number of random location tests.

   (run uniSphericalDeltas -h for usage summary).
 */

#define PGM_DSCR "UniSpherical coordinate encoding Δ\'s"
#define PGM_LAST_EDIT_DATE "2025.090"         /* format as from 'date +%Y.%j' */

#include <stdio.h>
#include <time.h>
#include <nemo.h>
#include "../scullions/scullions.h" /* include after nemo.h has been included */

#define TEST_NUMBER  10000000
void usage(const char *, const char *);
static const char *progName;                             /* messaging/logging */
/* -------------------------------------------------------------------------- */
int main (int argc,
          const char *argv[],
          const char *envr[]) {

   const char *optKey;
   const char *optVal;


   int n, testNum;
   nemoPtUs8 u64globLoc;  /* UniSpherical global location 8-byte unsigned int */
   nemoPtUs4 u32globLoc;                        /* ...and 4-byte unsigned int */

   nemoPtNcs nsFull, nsTruncated; /* NCS locations to be moved to and from it */

   double cdcDelta;                                         /* iu granularity */
   double cdcDeltaMax;
   double stDev;                                        /* standard deviation */
/* -------------------------------------------------------------------------- */
   progName = strrchr(argv[0], '/');                                 /* POSIX */
   if (progName == NULL) progName = strrchr(argv[0], '\\');         /* MS Win */
   if (progName) progName += 1;             /* skip found last path separator */
   if (progName == NULL) progName = argv[0];     /* neither? Just use argv[0] */
   fprintf(stderr, "\n[%s]: %s\nsource as of: %s, Nemo Library: %.3f\n",
                   progName, PGM_DSCR, PGM_LAST_EDIT_DATE, NEMO_LIBRARY_DATE);

   testNum = TEST_NUMBER;                    /* default random location tests */

   while ((optKey = clOption(argc, argv, &optVal))) {
      if (*optKey == 'h') usage(NULL, NULL);
      else if (*optKey == 'r') testNum = atoi(optVal);
      else usage("unrecognized option", optKey);
      }

   printf("Test with %.1f M random locations\n", (double)testNum / 1000000.0);
   printf("Direct/inverse 8-byte UniSpherical transformations:\n");
   stDev = 0.0;
   cdcDeltaMax = 0.0;
   srand(time(NULL));
   for (n = 0; n < testNum; n++) {
      if (n%1000000 == 0) fprintf(stderr, "%d M\r", n / 1000000);
      nsFull = *nemo_SphereRandomPointGlobal(NULL);
      u64globLoc = nemo_NcsToUs8(&nsFull);
      nemo_Us8ToNcs(u64globLoc, &nsTruncated);
      cdcDelta = sqrt(NEMO_ChordSq3(nsFull.dc, nsTruncated.dc)) * NEMO_EARTH_RADIUS;
      if (cdcDelta > cdcDeltaMax) cdcDeltaMax = cdcDelta;
      stDev += (cdcDelta * cdcDelta);
      }
   stDev = sqrt(stDev  / (double)(n - 1));
   printf("Δ max: %2d mm\n", (int)(cdcDeltaMax * 1000.0));
   printf("σ    : %2d mm\n", (int)(stDev * 1000.0));

   printf("Direct/inverse 4-byte UniSpherical transformations:\n");
   stDev = 0.0;
   cdcDeltaMax = 0.0;
   srand(time(NULL));
   for (n = 0; n < testNum; n++) {
      if (n%1000000 == 0) fprintf(stderr, "%d M\r", n / 1000000);
      nsFull = *nemo_SphereRandomPointGlobal(NULL);
      u32globLoc = nemo_NcsToUs4(&nsFull);
      nemo_Us4ToNcs(u32globLoc, &nsTruncated);
      cdcDelta = sqrt(NEMO_ChordSq3(nsFull.dc, nsTruncated.dc)) * NEMO_EARTH_RADIUS;
      if (cdcDelta > cdcDeltaMax) cdcDeltaMax = cdcDelta;
      stDev += (cdcDelta * cdcDelta);
      }
   stDev = sqrt(stDev  / (double)(n - 1));
   printf("Δ max: %3d m\n", (int)(cdcDeltaMax));
   printf("σ    : %3d m\n", (int)(stDev));
   }
/* ========================================================================== */
void usage(const char *mA,                  /* first message string (or NULL) */
           const char *mB) {               /* second message string (or NULL) */
   if (mA || mB) fprintf (stderr, "Error: %s %s\n", mA ? mA : "\0", mB ? mB : "\0");
   fprintf (stderr, "Usage: %s [options]\n", progName);
   fprintf (stderr, "Options:\n");
   fprintf (stderr, " -h(elp)      to print this usage help and exit\n");
   fprintf (stderr, " -r(andlocs)=nnnn: random locations to test (default:%d)\n", TEST_NUMBER);
   exit(1);
   }
/* ========================================================================== */
#include "../scullions/clFileOpt.c"
/* ========================================================================== */
