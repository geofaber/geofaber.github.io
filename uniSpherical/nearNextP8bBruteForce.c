/* nearNextP8bBruteForce: "Brute force" nearest neighbour TSP solution,
   starting with the first location in the file (a binary Us8 coordinate array
   read from input file given as the first command-line argument. The second
   command line argument is the name of the output file; the same binary Us8
   coordinates, re-ordered to form the the proposed  itinerary.

   Note the extreme simplicity of the algorithm!

   The program is invoked as:

   nearNextP8bBruteForce w1904711.p8b w1904711_nn.p8b
 */

#define PGM_DSCR "Brute-force nearest-next itinerary for .p8b input file"
#define PGM_LAST_EDIT_DATE "2025.085"
#define METERS2NM       0.0005399568
#include <stdio.h>
#include <time.h>

#include <nemo.h>
#include "../scullions/scullions.h" /* include after nemo.h has been included */

static const char *progName;    /* for error logging by this source file only */
/* ========================================================================== */
int main (int argc,
          const char *argv[],
          const char *envr[]) {

   int n, nn, nx;
   FILE *inFp;
   int inFileSize;
   int lcnCnt;              /* count of locations in input binary (.ptb) file */
   FILE *outFp;
   nemoPtUs8 location;
   nemoPtUs8 *locations;
   nemoPtUs8 holdLoc;
   nemoPtNcs ptNcs;                                     /* leg start location */
   nemoPtNcs cndNcs;                                           /* "candidate" */
   double cld;                               /* "candidate" location distance */
   double nld;                                   /* nearest location distance */
   nemoPtEnr ptEnr, ptEnrPrev;      /* ellipsoid leg start/end coordinates... */
   double totalLength, returnLegLength;    /* ...used for itinerary reporting */
   double clockSeconds, clockHours;                              /* timing... */
   time_t clockStart;                                      /* ...paraphenalia */
/* -------------------------------------------------------------------------- */
   progName = strrchr(argv[0], '/');                                 /* POSIX */
   if (progName == NULL) progName = strrchr(argv[0], '\\');         /* MS Win */
   if (progName == NULL) progName = argv[0];                      /* neither? */
   else progName += 1;                        /* strip leading path separator */
   fprintf(stderr, "\n[%s]: %s\nsource as of: %s, Nemo Library: %.3f\n",
                   progName, PGM_DSCR, PGM_LAST_EDIT_DATE, NEMO_LIBRARY_DATE);

   if (argc < 3) errorExit(progName, __LINE__,
      "command-line argumens: input.p8b output.p8b\n");

   inFp = fopen(argv[1], "rb");                            /* Open input file */
   if (inFp == NULL) errorExit(progName, __LINE__,
                                "Can't open [%s] for reading\n", argv[1]);
/* Input file is a flat array of Us8 coordinates */
/* Load the input file as memory-resident array, first get the file size:     */
   if (fseek(inFp, 0, SEEK_END)) errorExit(progName, __LINE__,
                                "Can't read [%s]?\n", argv[1]);
   inFileSize = ftell(inFp);
   rewind(inFp);
   lcnCnt = (int)(inFileSize / sizeof(nemoPtUs8));
   if (inFileSize%lcnCnt) errorExit(progName, __LINE__,
    "Input file size (%d) not multiple of %d\n", inFileSize, sizeof(nemoPtUs8));

   fprintf(stderr, "Input file has: %d records\n", (int)(inFileSize / sizeof(nemoPtUs8)));
   locations = malloc(lcnCnt * sizeof(location));
   if (locations == NULL) errorExit(progName, __LINE__, "No memory?\n");
   n = fread(locations, sizeof(nemoPtUs8), lcnCnt, inFp);
   fclose(inFp);
   fprintf(stderr, "Locations loaded: %d\n", n);

   clockStart = clock();
   for (n = 0; n < (lcnCnt - 2); n++) {
      if (n%100 == 0) fprintf(stderr, "Itinerary head at %d\r", n);
      nemo_Us8ToNcs(locations[n], &ptNcs); /* got leg start, find closest end */
      nld = NEMO_DOUBLE_HUGE;
      nx = -1;
      for (nn = n + 1; nn < lcnCnt; nn++) {
         nemo_Us8ToNcs(locations[nn], &cndNcs);
         cld = NEMO_ChordSq3(ptNcs.dc, cndNcs.dc);
         if (cld < nld) {
            nld = cld;
            nx = nn;
            }
         }
      if (nx == -1) errorExit(progName, __LINE__,                     /* WtF? */
               "Unexpected error while searching for next location at %d\n", n);
/*    found next location to visit. Swap it with n + 1, move to the next n */
      holdLoc = locations[n + 1];
      locations[n + 1] = locations[nx];
      locations[nx] = holdLoc;
      n++;
      }
   clockSeconds = (double)(clock() - clockStart) / (double)CLOCKS_PER_SEC;
   fprintf(stderr, "TSP itinerary sort of %d locations completed, duration: ", lcnCnt);
   clockHours = clockSeconds / (60.0 * 60.0);
   if (clockHours < 1) fprintf(stderr, "%.3f seconds\n", clockSeconds);
   else fprintf(stderr, "%.3f hours (%.3f seconds)\n", clockHours, clockSeconds);

/* report total itinerary length along geodesics */
   totalLength = 0.0;
   nemo_Us8ToNcs(locations[0], &ptNcs);
   nemo_NcsToEnr(nemo_ElrWgs84(), &ptNcs, &ptEnrPrev);
   for (n = 1; n < lcnCnt; n++) {
      nemo_Us8ToNcs(locations[n], &ptNcs);
      nemo_NcsToEnr(nemo_ElrWgs84(), &ptNcs, &ptEnr);
      totalLength += nemo_GeodesicSzpila(nemo_ElrWgs84(), &ptEnrPrev, &ptEnr, NULL);
      ptEnrPrev = ptEnr;
      }
   fprintf(stderr, "Open itinerary total: %12.3f\n", METERS2NM * totalLength);
   nemo_Us8ToNcs(locations[0], &ptNcs);
   nemo_NcsToEnr(nemo_ElrWgs84(), &ptNcs, &ptEnrPrev);
   returnLegLength = nemo_GeodesicSzpila(nemo_ElrWgs84(), &ptEnrPrev, &ptEnr, NULL);
   fprintf(stderr, "Return leg length: %12.3f\n", METERS2NM * returnLegLength);

/* write output file */
   outFp = fopen(argv[2], "wb");                /* Open output locations file */
   if (outFp == NULL) errorExit(progName, __LINE__,
                     "Can't open [%s] for writing itinerary sorted locations\n", argv[2]);

   n = fwrite(locations, sizeof(nemoPtUs8), lcnCnt, outFp);
   if (n != lcnCnt) errorExit(progName, __LINE__,
               "Error in writing itinerary sorted locations (record:%d)\n", n);
   fclose(outFp);
   free(locations);

   return(0);
   }
/* ========================================================================== */
#include "../scullions/errorExit.c"
/* ========================================================================== */
