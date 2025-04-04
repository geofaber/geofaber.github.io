/* r8bToAscii.c: Convert binary UniSpherical coordinate file into a text
   coordinate file. The text coordinates are by default simply the
   (single-number) hexadecimal representation of the input UniSpherical
   coordinate; optionally (see command-line options) it can be in the
   decimal degree φ and λ.

   Since the UniSpherical coordinates in all Nemo Library demo/example
   programs are on the near-conformal sphere, the transformation to φ and λ
   requires the knowledge of the ellipsoid parameters. (Important note: This
   is the case for ANY computation performed with φ, λ numbers!) This programs
   assumes the output coordinates (if φ, λ output was requuested) are on
   Wgs84 ellipsoid).

   The program assumes binary coordinate array might include "markers" that
   terminate the line segment or ring. Default format of marker records in
   binary UniSpherical file are 8-byte records (just like a binary form of
   UniSpherical coordinate number) with most significant nibble (half-byte)
   is 0 - this makes them "undefined" coordinate. If the rest of the record
   is not 0, it is assumed two items are packed as two unsigned 32-bit
   integers: the high order one being a line segment or ring unsigned integer
   id-number, the low order one the count of the vertices in the preceding
   segment.

   Note that binary input coordinate record is assumed to be of little-endian
   variety; as are all binary UniSpherical coordinates files in the Nemo
   Library API example prograsm collection.

Command line program invocation example:

   ./r8bToAscii xyz.r8b -n=100 -f=0 > xyz.pts

   First (and only) argument is the input .r8b (or .l8b, p8b) file of 8-byte
   UniSpherical binary records. Command line options (see "usage" or type
   "r8bToAscii -h") specify the limit of the input file top records to process
   and the coordinate output format integer specifier (0: hexadecimal
   UniSpherical - the default, n:decimal degree φ, λ where n is number of
   decimal fraction digits. Somewhat odd-valued decimal digits correspond
   to more commonly useful (Earth!) ground resolutions: 10 meters for 4
   and 1 millimeter for 8 decimal digits. Both are an order of magnitude
   higher than 4- and 8-byte (respectively) UniSpherical coordinates).

   Output can be redirected (see example above) for further processing.
 */

#define PGM_DSCR "Convert binary UniSpherical coordinates to text"
#define PGM_LAST_EDIT_DATE "2025.093"
#include <stdio.h>
#include <nemo.h>
#include "../scullions/scullions.h" /* include after nemo.h has been included */

static const char *progName;    /* for error logging by this source file only */
void usage(const char *, const char *);
/* ========================================================================== */
int main (int argc,
          const char *argv[],
          const char *envr[]) {
   int n, m;
   int iPlate, iFormat, nRecs, nCoords, nMarkers, idSeg, nSegPts;
   const char *optKey, *optVal;            /* options, in -keyword=value form */
   const char *inFn;                                       /* input file name */
   FILE *inFp;                  /* input binary file, coordinate to trabsform */
   nemoPtUs8 ptUs8;
   nemoPtEll locEll;
   nemoPtNcs locNcs;
/* -------------------------------------------------------------------------- */
   progName = strrchr(argv[0], '/');                                 /* POSIX */
   if (progName == NULL) progName = strrchr(argv[0], '\\');         /* MS Win */
   if (progName == NULL) progName = argv[0];                      /* neither? */
   else progName += 1;                        /* strip leading path separator */
   fprintf(stderr, "\n[%s]: %s\nsource as of: %s, Nemo Library: %.3f\n",
                   progName, PGM_DSCR, PGM_LAST_EDIT_DATE, NEMO_LIBRARY_DATE);

   iFormat = nRecs = 0; /* assume whole file, hexadecimal UniSpherical output */
   while ((optKey = clOption(argc, argv, &optVal))) {
      if (*optKey == 'h') usage(NULL, NULL);
      else if (*optKey == 'f') iFormat = atoi(optVal);
      else if (*optKey == 'n') nRecs = atoi(optVal);
      else usage("unrecognized option", optKey);
      }
   if (iFormat > 4) iFormat = 8;             /* about 1 mm along the meridian */
   else if (iFormat > 0) iFormat = 4;   /* about 10 meters along the meridian */
   else iFormat = 0;                       /* just to nix input of a negative */
/* fprintf(stderr, "limit: %d, format: %d\n", nRecs, iFormat); */

   inFn = clFileName(argc, argv);                          /* input file name */
   if (inFn == NULL) errorExit(progName, __LINE__,
                 "missing command line argument (input file name)\n");
   inFp = fopen(inFn, "rb");                            /* Open input file */
   if (inFp == NULL) errorExit(progName, __LINE__,
                                "Can't open [%s] for reading\n", inFn);
   n = nCoords = nMarkers = 0;
   m = fread(&ptUs8, sizeof(nemoPtUs8), 1, inFp);
/* fprintf(stderr, "first read: %d\n", m); */
   while (m) {
      if ((nRecs) && (n >= nRecs)) break;                     /* want no more */
      iPlate = NEMO_Us8Plate(ptUs8);
      if (iPlate == 0) {             /* line segment/ring end "marker" record */
         nMarkers++;
         idSeg = (int)(ptUs8 >> 32);
         nSegPts = (int)(ptUs8 & 0x00000000ffffffff);
         putchar('*');
         if ((idSeg) || (nSegPts)) printf(" %d %d", idSeg, nSegPts);
         putchar('\n');
         }
      else {                                /* UniSpherical coordinate record */
         nCoords++;
         if (iFormat) {                                    /* convert to φ, λ */
            nemo_Us8ToNcs(ptUs8, &locNcs);
            nemo_NcsToEll(nemo_ElrWgs84(), &locNcs, &locEll);
            if (iFormat > 4) printf("%12.8f %13.8f\n", NEMO_RAD2DEG * locEll.a[0],
                                                       NEMO_RAD2DEG * locEll.a[1]);
            else printf("%8.4f %9.4f\n", NEMO_RAD2DEG * locEll.a[0],
                                         NEMO_RAD2DEG * locEll.a[1]);
            }
         else printf("%016lx\n", ptUs8);
         }
      n++;
      m = fread(&ptUs8, sizeof(nemoPtUs8), 1, inFp);
      }
   fclose(inFp);

   fprintf(stderr, "%s done, coordinates: %d markers: %d\n",
                    progName, nCoords, nMarkers);
   return(0);
   }
/* ========================================================================== */
void usage(const char *mA,                  /* first message string (or NULL) */
           const char *mB) {               /* second message string (or NULL) */
   if (mA || mB) fprintf (stderr, "Error: %s %s\n", mA ? mA : "\0", mB ? mB : "\0");
   fprintf (stderr, "Usage: %s [options] inFile\n", progName);
   fprintf (stderr, "  inFile:  .r8b (or.l8b, .p8b) binary coordinate input file\n");
   fprintf (stderr, "Options:\n");
   fprintf (stderr, " -h[elp|  to print this usage help and exit\n");
   fprintf (stderr, " -f[ormat]=[0|4|8]: 0:hexUniS, n:φ,λ decimal° fraction digits\n");
   fprintf (stderr, " -n[umber]=nn restrict processing to first nn input records\n");
   exit(1);
   }
/* ========================================================================== */
#include "../scullions/errorExit.c"
#include "../scullions/clFileOpt.c"
/* ========================================================================== */
