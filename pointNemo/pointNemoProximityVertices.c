/* pointNemoProximityVertices.c: Find proximity vertices, with calculations
   using only spherical chord squared for distance criterion. This program can
   not be depended upon to find the three proximity vertices in all cases,
   this task requires some form of spatial sort and search capabilities
   implemented in the applications that may use the Nemo library functions as
   "building blocks" or, indeed a completely different fundamental approach to
   ellipsoidal geometry computing. This program only provides a heuristical
   solution that works well enough with physical coastline configuration in
   several tested geographical locations, most notably in the South Pacific
   region of planet Earth ("Point Nemo" or "The longest swim" georaphical
   trivia problem).

   In general, the input consists of:

   1) Binary file name containing coastline vertices' coordinates in the
      UniSpherical Us8 form.
   2) Coordinates of the centre of the search region.
   3) Radius of the search region as meters on the planetary surface.

   See usage for program input, output and options - for instance:

   ./pointNemoProxVertices osmPacific.p8b -center="-49.0, -123.4" -r=1e6

   At the end of the run, the three vertices' coordinates are written to the
   standard output, which can be re-directed to a file used as the input to
   the pointNemoIterate program by appending the above line, like so:

   ./pointNemoProximityVertices ... > proximityVertices.pts
*/

#include <time.h>

#define PGM_DSCR "Find three Nemo Proximity Vertices"
#define PGM_LAST_EDIT_DATE "2025.092"         /* format as from 'date +%Y.%j' */

#include <nemo.h>
#include "../scullions/scullions.h" /* include after nemo.h has been included */

#define MAX_COORD_STR                       64
#define TEST_COUNT                     2000000           /* what's a million? */
#define GLOBAL_LOCAL_RANDOM_CUTOFF     1500000            /* 1.5 K kilometers */
#define PROX_VRTX_SEPARATION              5000             /* five kilometers */
void usage(const char *, const char *);               /* program command-line */

static nemoPtNcs *cvx;                  /* a large array of coastlin vertices */
static const char *progName;    /* for error logging by this source file only */
/* ========================================================================== */
int main (int argc,
          const char *argv[],
          const char *envr[]) {

   const char *optKey;
   const char *optVal;
   const char *optValCenter;
   const char *optValRadius;
   const char *fnIn;            /* input binary file, OSM coastline (extract) */
   FILE *fpIn;

   int j, n;
   int testCount;
   int iPlate;
   int nCstVtx;                   /* number of points (vertices) on the coast */
   double globalLocalCutoff;                   /* for random point generation */
   double parms[NEMO_GNOMONIC_PCNT];
   double proxVrtxSeparation;                 /* vertex coincidence criterion */
   char coordStr[MAX_COORD_STR + 2];    /* text parsing, as simple as it gets */
   const char delimiters[] = ", ";
   char *token;
   nemoPtEll ptEll;           /* command linee input angular φ, λ coordinates */
   double srgnGround;         /* Search radius, as given on planeraty surface */
   double srgnArc;                      /* Search radius, as NCS arc (approx) */
   double srgnChSq;                    /* Search radius, as NCS chord squared */
   nemoPtNcs srgnCntr;                        /* search region centre, on NCS */

   size_t inFileSize;
   nemoPtUs8 ptUs8;
   int nTotalTests, moreTests;
   int nIn, nOut;
   double chSq, chSq0, chSq1;                  /* transient use square chords */
   nemoPtNcs ptRand;        /* random location, is it approximate Point Nemo? */
   nemoPtNcs ptNemo;
   double nemoDist;                /* ptNemo to nearest coast vertex, minimum */
   double randDist;          /* random point to nearest coast vertex, minimum */
   int i, ii, iii;
   nemoPtNcs proxVrtx[3];                         /* three proximity vertices */
   double proxVrtxChSq[3];             /* and chord squared distances to them */
   nemoPtNcs ptVrtx;                                      /* coastline vertex */
   double clockSeconds;                                /* timing paraphenalia */
   time_t clockStart;
/* -------------------------------------------------------------------------- */
   progName = strrchr(argv[0], '/');                                 /* POSIX */
   if (progName == NULL) progName = strrchr(argv[0], '\\');         /* MS Win */
   if (progName) progName += 1;             /* skip found last path separator */
   if (progName == NULL) progName = argv[0];     /* neither? Just use argv[0] */
   fprintf(stderr, "\n[%s]: %s\nsource as of: %s, Nemo Library: %.3f\n",
                   progName, PGM_DSCR, PGM_LAST_EDIT_DATE, NEMO_LIBRARY_DATE);

   if (argc < 2) usage("Missing command line argument(s)", NULL);

/* =============
   Preliminaries
   =============
 */
/* Retrieve command-line options and input filename */

   testCount = TEST_COUNT;                                         /* default */
   optValCenter = optValRadius = NULL;                    /* required options */
   while ((optKey = clOption(argc, argv, &optVal))) {
      if (*optKey == 'h') usage(NULL, NULL);
      else if (*optKey == 'c') optValCenter = optVal;
      else if (*optKey == 'r') optValRadius = optVal;
      else if (*optKey == 't') testCount=atoi(optVal);
      else usage("unrecognized option", optKey);
      }

   fnIn = clFileName(argc, argv);                      /* required input file */
   if (fnIn == NULL) usage("Missing input file name", NULL);

/* Search region centre coordinates */
   if (optValCenter == NULL) usage("Missing parameter - center coordinates", NULL);
   strncpy(coordStr, optValCenter, MAX_COORD_STR);
   token = strtok(coordStr, delimiters);
   ptEll.a[NEMO_LAT] = NEMO_DEG2RAD * strtod(token, NULL);
   token = strtok(NULL, delimiters);
   ptEll.a[NEMO_LNG] = NEMO_DEG2RAD * strtod(token, NULL);
/* Convert search region centre to NCS */
   nemo_EllToNcs(nemo_ElrWgs84(), &ptEll, &srgnCntr);

/* Search region radius */
   if (optValRadius == NULL) usage("Missing parameter - search radius", NULL);
   srgnGround = strtod(optValRadius, NULL);
/* fprintf(stderr, "Search region radis, arc on planet surface: %f\n", srgnGround); */
   srgnArc = srgnGround / NEMO_EARTH_RADIUS;         /* as arc on unit sphere */
   srgnChSq = nemo_ArcToChordApprox(srgnArc);
   srgnChSq *= srgnChSq;                                  /* as chord squared */

/* initialize local random point generation cutoff */
/* First, as arc on planet surface */
   globalLocalCutoff = GLOBAL_LOCAL_RANDOM_CUTOFF / NEMO_EARTH_RADIUS;
   globalLocalCutoff = nemo_ArcToChordApprox(globalLocalCutoff);  /* chord on unit sphere */
   globalLocalCutoff = globalLocalCutoff * globalLocalCutoff;  /* as chord squared */
   parms[0] = NEMO_DOUBLE_UNDEF;
   proxVrtxSeparation = PROX_VRTX_SEPARATION / NEMO_EARTH_RADIUS; /* arc on NCS */
   proxVrtxSeparation = nemo_ArcToChordApprox(proxVrtxSeparation); /* chord */
   proxVrtxSeparation = proxVrtxSeparation * proxVrtxSeparation;

/* Report all optios - as specified on the command line or defaults: */
   fprintf(stderr, "Input file: %s\n", fnIn);
   fprintf(stderr, "Search centre: %s\n", nemo_StrNcsCoords(&srgnCntr));
   fprintf(stderr, "Search region radius: %.0f\n", srgnGround);
   fprintf(stderr, "Vertex to vertex saparation criterion: %s\n",
                     nemo_StrChSqDist(proxVrtxSeparation));
/*
   fprintf(stderr, "Global/local randGen cutoff, given on planet: %.1f km\n",
                    GLOBAL_LOCAL_RANDOM_CUTOFF/1000.0);
   fprintf(stderr, "Global/local randGen cutoff, arc on NCS:      %.12f\n",
                    globalLocalCutoff);
   fprintf(stderr, "Global/local randGen cutoff cgSq:             %.12f\n", globalLocalCutoff);
   fprintf(stderr, "Global/local randGen cutoff: %.s\n", nemo_StrChSqDist(globalLocalCutoff));
 */

/* Open input file */
   fpIn = fopen(fnIn, "rb");
   if (fpIn == NULL) errorExit(progName, __LINE__,
                                "Can't open [%s] for reading\n", fnIn);

/* Load the input file as memory-resident array; first get the file size:     */
   if (fseek(fpIn, 0, SEEK_END)) errorExit(progName, __LINE__,
                                "Can't read [%s]?\n", fnIn);
   inFileSize = ftell(fpIn);
   rewind(fpIn);
   fprintf(stderr, "Input file has: %d records\n", (int)(inFileSize / sizeof(nemoPtUs8)));
   cvx = malloc((inFileSize / sizeof(nemoPtUs8)) * sizeof(nemoPtNcs));
   n = fread(&ptUs8, sizeof(nemoPtUs8), 1, fpIn);
   nCstVtx = 0;                                     /* count number of points */
   while (n) {
      iPlate = NEMO_Us8Plate(ptUs8);
      if (iPlate) nemo_Us8ToNcs(ptUs8, cvx + nCstVtx++);
      n = fread(&ptUs8, sizeof(nemoPtUs8), 1, fpIn);
      }
   fprintf(stderr, "Loaded search array of %d coastline vertices\n", nCstVtx);
   fclose(fpIn);

/* ==============================
   Phase 1: testing random points
   ==============================
 */
   fprintf(stderr, "Testing %d random points\n", testCount);
   moreTests = testCount;
   nIn = nOut = nTotalTests = n = 0;
   nemoDist = -(NEMO_DOUBLE_HUGE);
   clockStart = clock();
   while (moreTests) {              /* more random points remain to be tested */
      if (moreTests%1000 == 0) fprintf(stderr, "tests remaining: %d K      \r", moreTests/1000);
/*    If the search region is large, random points are generated as "global"
      and rejected if outside of it. Otherwise, random points will be
      generated as "local". (but the test will still be done). */
      if (srgnChSq > globalLocalCutoff) {
          ptRand = *nemo_SphereRandomPointGlobal(NULL);
//          fprintf(stderr, "Random Point global: %s\n", nemo_StrNcsCoords(&ptRand));
          }
      else {
         ptRand = *nemo_SphereRandomPointLocal(&srgnCntr, srgnArc, parms, NULL);
//         fprintf(stderr, "Random Point local: %s\n", nemo_StrNcsCoords(&ptRand));
         }
/*    its distance to search region centre: */
      chSq = NEMO_ChordSq3(ptRand.dc, srgnCntr.dc);
//    fprintf(stderr, "Distance from search region centre: %s\n", nemo_StrChSqDist(chSq));
/*    Reject generated random point if it is out of the search region */
      if (chSq > srgnChSq) nOut++;
      else nIn++;
      if (chSq > srgnChSq) continue;       /* ptRand is outside search region */
/*    This random point shiold be tested against all verices, and the closest
      one (distance and its index, ii recorded. We can abandon the traverse if
      the distance to a vertex is closer than the "best Nemo" found so far.   */
//      fprintf(stderr, "Random Point: %s\n", nemo_StrNcsCoords(&ptRand));
      nTotalTests++;                          /* used only to report progress */
      moreTests--;                                   /* next Monte Carlo test */
      randDist = NEMO_DOUBLE_HUGE;
      for (i = 0; i < nCstVtx; i++) {              /* traverse coast vertices */
         chSq = NEMO_ChordSq3(ptRand.dc, cvx[i].dc);      /* random to vertex */
//         fprintf(stderr, "Distance to vertex %3d: %s\n", i, nemo_StrChSqDist(chSq));
         if (chSq < nemoDist) {                       /* this random point... */
            randDist = 0.0;        /* ...can't possibly be next point nemo... */
            break;                    /* ...so break out of random point loop */
            }
         if (chSq < randDist) {
            randDist = chSq;                        /* nearest coastal so far */
            ii = i;                /* its vertex index */
            }
         }

      if (randDist > nemoDist) {                       /* better "Nemo" found */
         ptNemo = ptRand;
         nemoDist = randDist;
         iii = ii;
//         fprintf(stderr, "New Nemo: %3d %s\n", iii, nemo_StrNcsCoords(&ptRand));
//         fprintf(stderr, "New dstance: %s\n", nemo_StrChSqDist(randDist));
         }
      }

   if (nemoDist == -(NEMO_DOUBLE_HUGE)) errorExit(progName, __LINE__,
                        "Unexpected condition: no far point found?\n");
   fprintf(stderr, "Approximate Point Nemo  %s\n", nemo_StrNcsCoords(&ptNemo));
   fprintf(stderr, "Near coast vertex index: %d\n", iii);

   ptVrtx = cvx[iii];                      /* furthest vertex Ncs coordinates */
   fprintf(stderr, "Near coast vertex:      %s\n", nemo_StrNcsCoords(&ptVrtx));
   chSq = NEMO_ChordSq3(ptVrtx.dc, ptNemo.dc);            /* vertex to random */
   fprintf(stderr, "Distance to it: %s\n", nemo_StrChSqDist(chSq));
   proxVrtx[0] = ptVrtx;
   proxVrtxChSq[0] = chSq;

/* verification pass: is it really the closest one? */
   for (i = 0; i < nCstVtx; i++) {                 /* traverse coast vertices */
      chSq = NEMO_ChordSq3(cvx[i].dc, ptNemo.dc);     /* vertex to point Nemo */
      if (chSq < proxVrtxChSq[0]) {
         errorExit(progName, __LINE__, "Unexpected distance: vertex %d, %s\n",
         i, nemo_StrChSqDist(chSq));
         }
      }
   clockSeconds = (double)(clock() - clockStart) / (double)CLOCKS_PER_SEC;
   fprintf(stderr, "Phase 1 duration: %6.3f seconds, verification passed\n\n", clockSeconds);

/* ===============================================
   Phase 3: Find three closest points on the coast
   ===============================================
 */
/* proxVrtx[0] (and its associated distance) has been identified and set.
   Search for [1] and [2] will be done in two separate passes of coastline.
 */
   proxVrtxChSq[1] = NEMO_DOUBLE_HUGE;           /* search for second closest */
   for (i = 0; i < nCstVtx; i++) {             /* traverse all coast vertices */
      chSq = NEMO_ChordSq3(ptNemo.dc, cvx[i].dc);     /* point Nemo to vertex */
      if (chSq > proxVrtxChSq[1]) continue;   /* we already have a better one */
      chSq0 = NEMO_ChordSq3(cvx[i].dc, proxVrtx[0].dc);
      if (chSq0 < proxVrtxSeparation) continue;    /* close to [0], ignore it */
      proxVrtx[1] = cvx[i];                        /* record proximity vertex */
      proxVrtxChSq[1] = chSq;                           /* and distance to it */
      }

   proxVrtxChSq[2] = NEMO_DOUBLE_HUGE;            /* search for third closest */
   for (i = 0; i < nCstVtx; i++) {             /* traverse all coast vertices */
      chSq = NEMO_ChordSq3(ptNemo.dc, cvx[i].dc);     /* point Nemo to vertex */
      if (chSq > proxVrtxChSq[2]) continue;   /* we already have a better one */
      chSq0 = NEMO_ChordSq3(cvx[i].dc, proxVrtx[0].dc);
      chSq1 = NEMO_ChordSq3(cvx[i].dc, proxVrtx[1].dc);
      if ((chSq0 < proxVrtxSeparation) ||               /* close to [0] or... */
          (chSq1 < proxVrtxSeparation)) continue;        /* ...[1], ignore it */
      proxVrtx[2] = cvx[i];                        /* record proximity vertex */
      proxVrtxChSq[2] = chSq;                           /* and distance to it */
      }

   fprintf(stderr, "Phase 2 done\n");

/* all that's left to do, is to report the results: */
   nemo_NcsToEll(nemo_ElrWgs84(), &ptNemo, &ptEll);
   fprintf(stdout, "# approximate point Nemo φ, λ: %13.9f,%14.9f\n",
           NEMO_RAD2DEG * ptEll.a[0], NEMO_RAD2DEG * ptEll.a[1]);
   fprintf(stdout, "# three proximity vertices and distances to them:\n");
   for (j = 0; j < 3; j++) {                 /* traverse 3 proximity vertices */
      nemo_NcsToEll(nemo_ElrWgs84(), proxVrtx + j, &ptEll);
      fprintf(stdout, "%11.7f,%12.7f, %s\n",
              NEMO_RAD2DEG * ptEll.a[NEMO_LAT], NEMO_RAD2DEG * ptEll.a[NEMO_LNG],
              nemo_StrChSqDist(proxVrtxChSq[j]));
      }

   free(cvx);
   return(0);
   }
/* ========================================================================== */
void usage(const char *mA,                  /* first message string (or NULL) */
           const char *mB) {               /* second message string (or NULL) */

   if (mA || mB) fprintf (stderr, "Error: %s %s\n", mA ? mA : "\0", mB ? mB : "\0");
   fprintf (stderr, "Usage: %s [options] inFile\n", progName);
   fprintf (stderr, "  inFile: Binary .ptb input file of coastline vertices\n");
   fprintf (stderr, "Required options:\n");
   fprintf (stderr, " -c(enter)=\"φ,λ\": coordinate string, center of search area\n");
   fprintf (stderr, " -r(adius)=rrrr: meters, search radius\n");
   fprintf (stderr, "Other options:\n");
   fprintf (stderr, " -t(estcount)=nnn: integer, random test count, (default:%d)\n", TEST_COUNT);
   fprintf (stderr, " -h(elp): to print this usage help and exit\n");
   exit(1);
   }
/* ========================================================================== */
#include "../scullions/clFileOpt.c"
#include "../scullions/errorExit.c"
#include "../scullions/nemoStrings.c"
/* ========================================================================== */
