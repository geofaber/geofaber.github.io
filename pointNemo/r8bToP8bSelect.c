/* r8bToP8bSelect.c: Program to select/extract points from a large
   agglomerate set of point coordinates in UniSpherical (Us8) binary format.
   Line segment/ring break markers, if present, are ignored.

   The extraction criterion is an "ellipsoid circle"; i.e. extracted points
   must be closer to the given extraction center than a given ~geodesic~
   distance. Note that this level of precision is not usually required; the
   purpose of the code presented here is to demonstrate that even such high
   accuracy computations can be carried out very efficiently by a purposeful
   combination of spherical and ellipsoidal geometry productions).

   The command-line program is invoked with two (ordered) filenames and two
   (mandatory) program execution options - for instance:

   ./extractFromR8b input.r8b output.p8b -center="41.5,18.1" -radius=150000

   where:

      input.r8b
         input file in UniSpherical binary coordinate format. File format is
         that used by GGW (Galileo Geodetic Workbench) Note that the file can
         be either a point (.ptb) a line segment (.lnb) or a region (.rgb) file.

      output.p8b
         output point agglomerate file with binary UnSpherical coordinates.
         (Note: current version of this program writes no line segment or
         ring terminator markers - thus it is only capable of selecting
         0-dimensional seb-sets (.p8b files, see GGW documentation).

      -center
         extraction center quote-string of comma-separated φ, λ coordinates
         in angular decimal degrees - for instance "-50.0000003, -125.0000007"
         for a point in South Pacific far from any land.

      -radius
         extraction distance as a length of ellipsoid geodesic measured in
         meters - for instance, 3000000.003 (in this example three thousand
         kilometres and a little bit more) on the planetary surface).
 */
#include <time.h>

#include <nemo.h>
#include "../scullions/scullions.h" /* include after nemo.h has been included */

#define PGM_DSCR "Extraction of point records from .ptb/.lnb file"
#define PGM_LAST_EDIT_DATE "2025.091"         /* format as from 'date +%Y.%j' */

#define BLOCK_POINTS   1024                   /* reading/writing is in blocks */
#define MAX_COORD_STR   128

/* pending inclusion to nemo.h */
//#define NEMO_Us8Plate(u8) ((int)((u8 & 0xf000000000000000) >> 60))

int proxChordTest(nemoPtNcs *, nemoPtNcs *, double, double);
int proxGeodesicTest(nemoPtEnr *, nemoPtEnr *, double);

static const char *progName;    /* for error logging by this source file only */
void usage(const char *, const char *);
/* ========================================================================== */
int main (int argc,
          const char *argv[],
          const char *envr[]) {

   int i, iPlate, n, nb;
   int isClose;                         /* 1:is close, -1:is far, 0:uncertain */
   int nBin, nBout;             /* number of points in input and output block */
   int nMarksIn;                /* number of segment/ring marks in input file */
   int nGeodTests;       /* number of point classified by geodesic evaluation */
   int nPtIn, nPtOut, nPtFar;       /* number of input, output and far points */
   const char *optKey, *optVal;            /* options, in -keyword=value form */
   const char *fnIn, *fnOut;                  /* as given on the command line */
   FILE *fpIn;                              /* input: .ptb, .lnb or .rgb file */
   nemoPtUs8 ptUs8in[BLOCK_POINTS];        /* input block of U64 (CDC) points */
   FILE *fpOut;                                           /* output .ptb file */
   nemoPtUs8 ptUs8out[BLOCK_POINTS];                          /* output block */
   char coordStr[MAX_COORD_STR + 2];    /* text parsing, as simple as it gets */
   const char delimiters[] = ", \r\n";
   const char *strCenter, *strRadius;
   char *token;
   nemoPtEll ptEll;            /* command line input angular φ, λ coordinates */
   nemoPtEnr rtcEnr;                 /* retrieval center, as ellipsoid normal */
   nemoPtNcs rtcNcs;             /* retrieval center as near-conformal sphere */

   nemoPtNcs ptNcs;              /* input file point on near-conformal sphere */
   nemoPtEnr ptEnr;                   /* input file point as ellipsoid normal */

   double exRadGeodesic;                       /* extraction radius, geodesic */
   double chSqNear, chSqFar;             /* chord squared inclusion/exclusion */
   double a, c;                                  /* ...to calculate the above */

   double clockSeconds;                               /* timing paraphernalia */
   time_t clockStart;
/* -------------------------------------------------------------------------- */
   progName = strrchr(argv[0], '/');                                 /* POSIX */
   if (progName == NULL) progName = strrchr(argv[0], '\\');         /* MS Win */
   if (progName) progName += 1;             /* skip found last path separator */
   if (progName == NULL) progName = argv[0];     /* neither? Just use argv[0] */
   fprintf(stderr, "\n[%s]: %s\nsource as of: %s, Nemo Library: %.3f\n",
                   progName, PGM_DSCR, PGM_LAST_EDIT_DATE, NEMO_LIBRARY_DATE);

   strCenter = strRadius = NULL;
   while ((optKey = clOption(argc, argv, &optVal))) {
      if (*optKey == 'h') usage(NULL, NULL);
      else if (*optKey == 'c') strCenter = optVal;
      else if (*optKey == 'r') strRadius = optVal;
      else usage("unrecognized option", optKey);
      }

/* Extract retrieval center φ, λ coordinates */
   if (strCenter == NULL) usage("missing argument:", "extraction center");
   strncpy(coordStr, strCenter, MAX_COORD_STR);
   token = strtok(coordStr, delimiters);
   ptEll.a[NEMO_LAT] = NEMO_DEG2RAD * strtod(token, NULL);
   token = strtok(NULL, delimiters);
   ptEll.a[NEMO_LNG] = NEMO_DEG2RAD * strtod(token, NULL);
   fprintf(stderr, "Retrieval center, φ,λ: %.7f, %.7f\n",
           NEMO_RAD2DEG * ptEll.a[NEMO_LAT], NEMO_RAD2DEG * ptEll.a[NEMO_LNG]);

   nemo_LatLongToDcos3(ptEll.a, rtcEnr.dc); /* Convert to ellipsoid normal... */
   nemo_EnrToNcs(nemo_ElrWgs84(), &rtcEnr, &rtcNcs);      /* ...and NC sphere */
/* fprintf(stderr, "Retrieval center NCS: %f, %f %f\n",
           rtcNcs.dc[0], rtcNcs.dc[1], rtcNcs.dc[2]); */

/* Extract retrieval radius as geodesic, meters on the surface */
   if (strRadius == NULL) usage("missing argument:", "extraction radius");
   exRadGeodesic = strtod(strRadius, NULL);
   fprintf(stderr, "Retrieval radius: %.0f meters\n", exRadGeodesic);

/* First file argument: input file path/name */
   fnIn = clFileName(argc, argv);                               /* input file */
   if (fnIn == NULL) usage("Missing input file name", NULL);
   fpIn = fopen(fnIn, "rb");
   if (fpIn == NULL) errorExit(progName, __LINE__,
                                "Can't open [%s] for reading\n", fnIn);
   fprintf(stderr, "Input from: [%s]\n", fnIn);

/* Second file argument: output file path/name */
   fnOut = clFileName(argc, argv);                             /* output file */
   if (fnOut == NULL) usage("Missing output file name", NULL);
   fpOut = fopen(fnOut, "wb");
   if (fpOut == NULL) errorExit(progName, __LINE__,
                                "Can't open [%s] for writing\n", fpOut);
   fprintf(stderr, "Output to: [%s]\n", fnOut);

/* Find squared chord magnitude below which the point is included, and the
   one above which it is rejected. For points with squared chord distances
   between the two values, the more expensive geodesic length computation
   will be required in order to decide whether to include or reject. */
   a = NEMO_GEOARC_MIN * (exRadGeodesic / NEMO_EARTH_RADIUS);
   c = nemo_ArcToChordApprox(a);
   chSqNear = c * c;
   a = NEMO_GEOARC_MAX * (exRadGeodesic / NEMO_EARTH_RADIUS);
   c = nemo_ArcToChordApprox(a);
   chSqFar = c * c;
   fprintf(stderr, "geodesic, (NCS limits): %f, (%f, %f)\n", exRadGeodesic, chSqNear, chSqFar);

   nPtIn = nPtOut = nPtFar = nBin = nBout = nGeodTests = 0;
   nBin = fread(ptUs8in, sizeof(nemoPtUs8), BLOCK_POINTS, fpIn);
   nPtIn += nBin;
   nb = 0;

   clockStart = clock();
   while (nBin > 0) {
    if (nb++%1000 == 0) fprintf(stderr, "%d M\r", nb / 1000);
    for (i = 0; i < nBin; i++) {            /* traverse points in input block */
         iPlate = NEMO_Us8Plate(ptUs8in[i]);
         if (iPlate == 0) {
            nMarksIn++;
            continue;
            }
/*       Determine if the point is within the given proximity; if it is,
         transfer the coordinates to the next free slot in the output block */
         isClose = 0;                           /* unknown, we must determine */
/*       Point coordinates on near-conformal sphere - fast transformation */
         nemo_Us8ToNcs(ptUs8in[i], &ptNcs);
/*       For curios cats: comment out the following statement,
         recompile and observe the change in reported duration */
         isClose = proxChordTest(&rtcNcs, &ptNcs, chSqNear, chSqFar);
         if (isClose == 0) {   /* chord proximity test did not provide answer */
/*          Somewhat more expensive transformation of point coordinates to
            the ellipsoid, followed by a much more expensive geodesic test... */
            nemo_NcsToEnr(nemo_ElrWgs84(), &ptNcs, &ptEnr);   /* to ellipsoid */
            isClose = proxGeodesicTest(&rtcEnr, &ptEnr, exRadGeodesic);
            }
         if (isClose > 0) {      /* point is close, must be written to output */
            ptUs8out[nBout++] = ptUs8in[i];
            if (nBout == BLOCK_POINTS) {        /* tile to write output block */
               n = fwrite(ptUs8out, sizeof(nemoPtUs8), BLOCK_POINTS, fpOut);
               if (n != BLOCK_POINTS) errorExit(progName, __LINE__,
                   "write error at %d input, %d ouput record\n", nPtIn, nPtOut);
               nBout = 0;                        /* output block is now empty */
               }
            nPtOut++;
            }
         else nPtFar++;                     /* point is far, was not included */
         }
      nBin = fread(ptUs8in, sizeof(nemoPtUs8), BLOCK_POINTS, fpIn);
      nPtIn += nBin;
      }

   clockSeconds = (double)(clock() - clockStart) / (double)CLOCKS_PER_SEC;
   printf("duration: %6.3f seconds\n", clockSeconds);

   if (nBout > 0) {                    /* write partially filled block and... */
      n = fwrite(ptUs8out, sizeof(nemoPtUs8), nBout, fpOut);
      if (n != nBout) errorExit(progName, __LINE__,
          "write error at %d input, %d ouput record\n", nPtIn, nPtOut);
      nBout = 0;            /* not that it matters, but output block is empty */
      }

   fclose(fpIn);
   fclose(fpOut);

   fprintf(stderr, "Input points (records):     %8d\n", nPtIn);
   fprintf(stderr, "Input segments or rings:    %8d\n", nMarksIn);
   fprintf(stderr, "Points included:            %8d\n", nPtOut);
   fprintf(stderr, "Points excluded:            %8d\n", nPtFar);
   fprintf(stderr, "Geodesic tests required:    %8d\n", nGeodTests);
   return(0);
   }
/* ========================================================================== */
/* Determine point proximity based on quick (but possibly inconclusive)
   spherical square proximity. Return 1 for close, -1 for far, 0 for
   undetermined.
 */
int proxChordTest(nemoPtNcs *ncsA, nemoPtNcs *ncsB,       /* two given points */
                  double chSqNear, double chSqFar) {     /* near/far criteria */
   double chSq;
   chSq = NEMO_ChordSq3(ncsA->dc, ncsB->dc);
   if (chSq < chSqNear) return(1);                                  /* inside */
   if (chSq > chSqFar) return(-1);                                 /* outside */
   return(0);
   }
/* ========================================================================== */
/* Determine point proximity based on rigorous geodesic evaluation. Return
   1 for close, -1 for far. (ε is so minuscule we can - somewhat arbitrary -
   consider equal length to be "in").
 */
int proxGeodesicTest(nemoPtEnr *enrA, nemoPtEnr *enrB,    /* two given points */
                     double geodesic) {      /* proximity criterion, geodesic */
   double g;
   g = nemo_GeodesicSzpila(nemo_ElrWgs84(), enrA,  enrB, NULL);
   if (g == NEMO_DOUBLE_UNDEF) errorExit(progName, __LINE__,
                            "Unexpected Vincenty failure\n");
   if (g > geodesic) return(-1);
   return(1);
   }
/* ========================================================================== */
void usage(const char *mA,                  /* first message string (or NULL) */
           const char *mB) {               /* second message string (or NULL) */
   if (mA || mB) fprintf (stderr, "Error: %s %s\n", mA ? mA : "\0", mB ? mB : "\0");
   fprintf (stderr, "Usage: %s [options] inFile outFile\n", progName);
   fprintf (stderr, "  inFile:  .r8b (or.lnb, .p8b) binary coordinate input file\n");
   fprintf (stderr, "  outFile: .p8b binary coordinate output file\n");
   fprintf (stderr, "Options:\n");
   fprintf (stderr, " -h[elp|  to print this usage help and exit\n");
   fprintf (stderr, " -c[enter]=\"φ,λ\" extraction center, in decimal degrees\n");
   fprintf (stderr, " -r[adius]=nnn extraction radius, meters on planetary surface\n");
   exit(1);
   }
/* ========================================================================== */
#include "../scullions/errorExit.c"
#include "../scullions/clFileOpt.c"
/* ========================================================================== */
