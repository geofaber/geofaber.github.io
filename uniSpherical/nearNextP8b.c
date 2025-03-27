/* nearCs8.c Create pseudo-nearest-neighbor itinerary throughout locations
   organized in a sorted binary .ptb file - an array of Cs8 point locations.
   The program takes advantage of the fact that a Cs8 point coordinate array
   sorted on the coordinate numeric value keeps - to the maximum extent
   possible - the locations close to each other on spherical (or ellipsoidal)
   surface close to each other in the numerically ordered coorsinate array.
   This is achieved by a "two-level-search": first within a restricted segment
   of the coordinate array, and when no "un-visited" locations ate found in it,
   searching for the nearest location outside the segment.

   First command line argument is the input (one to be sorted in TSP-like
   itinerary order) .ptb file, the second is the (itinerary order sorted)
   output. It is assumed the input file is either ordered by natural Cs8
   coordinate numeric order, or in some itinerary order that can be further
   improved. The third argument is an integer, the size of "hopping window"
   that is used limit the initial search for the nearest neighbour.

   For instance:

   nearCs8 w1904711.ptb w1904711Itin_A.ptb 1000
 */

#define PGM_DSCR "Itinerary from Cs8 sorted binary (.ptb) file"
#define PGM_LAST_EDIT_DATE "2025.042"
#include <stdio.h>
#include <time.h>

#include <nemo.h>
#include "../scullions/scullions.h" /* include after nemo.h has been included */

#define METERS2NM       0.0005399568
#define MIN_WIN        16         /* the low value only for testing/debbuging */
#define MAX_WIN     32000

static int closeInWin(int);
static int closeOutWin(int);
static int compIntsIx0(const void *, const void *);

struct loc {                                 /* locations to be "TSP ordered" */
   int iOrd;             /* sort order - fist the cluster, then the itinerary */
   int unused;
   nemoPtCs8 ptCs8;                                   /* location coordinates */
   };

static int loCnt;                                      /* number of locations */
static struct loc *locs;                                         /* locations */
static int iWin;       /* first-pass nerest neighbour search window (half of) */
static const char *progName;    /* for error logging by this source file only */
static int nInsideWin, nOutsideWin;          /* to report algorithm behaviour */
/* ========================================================================== */
int main (int argc,
          const char *argv[],
          const char *envr[]) {

   int k, n;
   int nPrev, nNext;
   int inFileSize;
   FILE *inFp;                       /* input binary file, locations to visit */
   FILE *outFp;                                    /* itinerary-sorted output */
   nemoPtCs8 locCs8;                              /* location Cs8 coordinates */
   nemoPtNcs ptNcs;                     /* used only in itinerary report pass */
   nemoPtEnr ptEnr, ptEnrPrev;                                       /* ditto */
   double totalLength, returnLegLength;                              /* ditto */
   double clockSeconds;                                          /* timing... */
   time_t clockStart;                                      /* ...paraphenalia */
/* -------------------------------------------------------------------------- */
   progName = strrchr(argv[0], '/');                                 /* POSIX */
   if (progName == NULL) progName = strrchr(argv[0], '\\');         /* MS Win */
   if (progName == NULL) progName = argv[0];                      /* neither? */
   else progName += 1;                        /* strip leading path separator */
   fprintf(stderr, "\n[%s]: %s\nsource as of: %s, Nemo Library: %.3f\n",
                   progName, PGM_DSCR, PGM_LAST_EDIT_DATE, NEMO_LIBRARY_DATE);

   if (argc < 4) errorExit(progName, __LINE__,
      "command-line arguments: w1904711.ptb w1904711_nn.itin\n");

   fprintf(stderr, "Binary input from: %s\n", argv[1]);
   inFp = fopen(argv[1], "rb");                  /* Open input locations file */
   if (inFp == NULL) errorExit(progName, __LINE__,
                                "Can't open [%s] for reading locations\n", argv[1]);

   fprintf(stderr, "Binary output to: %s\n", argv[2]);
   outFp = fopen(argv[2], "wb");                /* Open output locations file */
   if (outFp == NULL) errorExit(progName, __LINE__,
                                "Can't open [%s] for writing locations\n", argv[2]);
   fclose(outFp);              /* we'll open it again when it's time to writa */

   n = atoi(argv[3]);
   if ((n < MIN_WIN) || (n > MAX_WIN)) errorExit(progName, __LINE__,
       "invalid window size (%d < n < %)\n", MIN_WIN, MAX_WIN);
   fprintf(stderr, "Search window :%d\n", n);
   iWin = n / 2;   /* half "below" and half "above" the last visited location */

/* Load locations into memory-resident array of location structs */
   if (fseek(inFp, 0, SEEK_END)) errorExit(progName, __LINE__,
                                "Can't read [%s]?\n", argv[1]);
   inFileSize = ftell(inFp);
   rewind(inFp);
   loCnt = (int)(inFileSize / sizeof(nemoPtCs8));
   if (inFileSize%loCnt) errorExit(progName, __LINE__,
    "Input file size (%d) not multiple of %d\n", inFileSize, sizeof(nemoPtCs8));

   fprintf(stderr, "Input file has: %d records\n",
                    (int)(inFileSize / sizeof(nemoPtCs8)));

   locs = malloc(loCnt * sizeof(struct loc));
   if (locs == NULL) errorExit(progName, __LINE__, "No memory for locations?\n");

   for (n = 0; n < loCnt; n++) {
      fread(&locCs8, sizeof(nemoPtCs8), 1, inFp);
      locs[n].iOrd = -1;                            /* set all to "unvisited" */
      locs[n].ptCs8 = locCs8;
/*    fprintf(stderr, "%d %s\n", n, nemo_StrCs8Coords(locs[n].ptCs8)); */
      }
   fclose(inFp);
   fprintf(stderr, "Locations loaded: %d\n", n);

   nInsideWin = nOutsideWin = 0;
   clockStart = clock();                                     /* time TSP sort */
   locs[0].iOrd = 0;
   nPrev = 0;                               /* index of last visited location */
   k = 1;
   while (k < loCnt) {
     if (k%1000 == 0) fprintf(stderr, "Itinerary stations: %dK\r", k);
      nNext = closeInWin(nPrev); /* Search for closest inside "search window" */
/*    fprintf(stderr, "next: %2d prev: %2d, in Win next: %2d\n", k, nPrev, nNext); */
      if (nNext == -1) nNext = closeOutWin(nPrev);  /* none found, go outside */
      if (nNext == -1) errorExit(progName, __LINE__, "Program assertion?\n");
      locs[nNext].iOrd = k++;  /* assign itinerary visitation order to loc... */
      nPrev = nNext;                              /* ...and resume the search */
      }
   qsort(locs, loCnt, sizeof(struct loc), compIntsIx0);    /* sort for output */
   clockSeconds = (double)(clock() - clockStart) / (double)CLOCKS_PER_SEC;
   fprintf(stderr, "Cc8 coordinates itinerary ordering  %6.3f seconds\n", clockSeconds);

   fprintf(stderr, "found inWin: %d, found outWin %d\n", nInsideWin, nOutsideWin);

   fprintf(stderr, "Locations sorted: %d\n", k);

/* report total itinerary length along geodesics */
   totalLength = 0.0;
    nemo_Cs8ToNcs(locs[0].ptCs8, &ptNcs);
    nemo_NcsToEnr(nemo_ElrWgs84(), &ptNcs, &ptEnrPrev);
    for (n = 1; n < loCnt; n++) {
       nemo_Cs8ToNcs(locs[n].ptCs8, &ptNcs);
       nemo_NcsToEnr(nemo_ElrWgs84(), &ptNcs, &ptEnr);
       totalLength += nemo_GeodesicSzpila(nemo_ElrWgs84(), &ptEnrPrev, &ptEnr, NULL);
       ptEnrPrev = ptEnr;
       }
   fprintf(stderr, "Open itinerary total: %12.3f\n", METERS2NM * totalLength);
   nemo_Cs8ToNcs(locs[0].ptCs8, &ptNcs);
   nemo_NcsToEnr(nemo_ElrWgs84(), &ptNcs, &ptEnrPrev);
   returnLegLength = nemo_GeodesicSzpila(nemo_ElrWgs84(), &ptEnrPrev, &ptEnr, NULL);
   fprintf(stderr, "Return leg length:    %12.3f\n", METERS2NM * returnLegLength);

/* write output file */
   outFp = fopen(argv[2], "wb");                /* Open output locations file */
   if (outFp == NULL) errorExit(progName, __LINE__,
                     "Can't open [%s] for writing itinerary sorted locations\n", argv[2]);
   for (n = 0; n < loCnt; n++) {
      locCs8 = locs[n].ptCs8;
      fwrite(&locCs8, sizeof(nemoPtCs8), 1, outFp);
      }
   fclose(outFp);
   free(locs);

   printf("Itinerary total, nautical miles: %.3f; Sort duration: %.3f\n",
           METERS2NM * (totalLength + returnLegLength), clockSeconds);
   return(0);
   }
/* ========================================================================== */
/* Search for the closest location inside the window +/- slots from nLast
 */
static int closeInWin(int nLast) {
   int n, nStart, nEnd, nMin;
   double csqMin, csq;                       /* square chord search distances */
   nemoPtNcs ptNcsLast, ptNcs;  /* ncs coords of the "stable", "moving" point */
/* -------------------------------------------------------------------------- */
   nStart = nLast - iWin;
   if (nStart < 0) nStart = 0;
   nEnd = nLast + iWin;
   if (nEnd > loCnt) nEnd = loCnt;
/* fprintf(stderr, "closeInWin, from %d to %d\n", nStart, nEnd); */
   nemo_Cs8ToNcs(locs[nLast].ptCs8, &ptNcsLast);
   csqMin = NEMO_DOUBLE_HUGE;
   nMin = -1;
   for (n = nStart; n < nEnd; n++) {
      if (locs[n].iOrd != -1) continue;          /* this location was visited */
/*    fprintf(stderr, "candidate next %d\n", n); */
      nemo_Cs8ToNcs(locs[n].ptCs8, &ptNcs);
      csq = NEMO_ChordSq3(ptNcsLast.dc, ptNcs.dc);
      if (csq < csqMin) {            /* new closest in-window point was found */
         csqMin = csq;
         nMin = n;                               /* location index of closest */
/*       fprintf(stderr, "best so far %d\n", nMin); */
         }
      }
   nInsideWin++;
   return(nMin);
   }
/* ========================================================================== */
/* Search for the next point to resume itinerary outside the "close search"
   window, alternating between high and low side of last visited location.
   Not finding a point would be an obvious "two-level-search" algorithm error.
 */
static int closeOutWin(int nLast) {
   int nLow, nHigh;
/* -------------------------------------------------------------------------- */
   nLow = nLast - iWin;
   if (nLow < 0) nLow = 0;
   nHigh = nLast + iWin;
   if (nHigh > loCnt) nHigh = loCnt;
   nOutsideWin++;
   while ((nLow > 0) || (nHigh < loCnt)) {
      if (nLow > 0) {
         if (locs[nLow].iOrd == -1) return(nLow);
         nLow--;
         }
      if (nHigh < loCnt) {
         if (locs[nHigh].iOrd == -1) return(nHigh);
         nHigh++;
         }
      }
   return(-1);                                     /* this better not happen! */
   }
/* ========================================================================== */
/* Natural integer compare function. Two structures are assumed to start with
   an array of natural integers. Compare them (in a manner required by C
   standard library qsort() and bsearch() functions. This implementation
   assumes first (ix:0) element of the array of int's 'determines the order. */
#define COMP_INDEX 0
int compIntsIx0(const void *p1, const void *p2) {
   return(((int *)p1)[COMP_INDEX] - ((int *)p2)[COMP_INDEX]);
   }
/* ========================================================================== */
#include "../scullions/errorExit.c"
#include "../scullions/nemoStrings.c"
/* ========================================================================== */
