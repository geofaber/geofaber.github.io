/* pointNemoDisqualify.c: Traverse a file of coast vertex coordinates and
   attemp to disqualify the proposed solution by finding that there are not
   exactly three points at equal geodesic distance from it, or that there
   is a point or points closer to it than the proposed solution distance.
   ===========================================================================
   This test does not prove that the proposed Point Nemo geometry solution
   is correct, but - if a closer point or points are found - it does prove
   that the proposed solution IS NOT correct. The feasibility of opposite
   proof, i.e., that a proposed solution has no alternative with shorter
   Nemo distance that would pass this test) depends on the coastline
   configuration: from simple and obvious to computationally infeasible.
   ===========================================================================
   The program is somewhat similar to the one that extracts the coastal
   points from large OSM .p8b file; it is however assumed this program
   will traverse a smaller set of points, in a smaller file, and that it
   is quite important to keep the code as simple as possible: the
   disqualification verdict must be simple in order to be convincing.

   Input file is assumed to be the same as the one that was used to compute
   the solution to the "longest swim problem": coordinates of the Point Nemo,
   its distance to the nearest point on land and three "proximity vertices".

   The program takes three command-line arguments: the file name of coastline
   points UniSpherical coordinates used to find Nmeo Proximity vertices), the
   coordinates of the proposed Point Nemo and the proposed Nemo Distance;
   for example:

   ./pointNemoDisqualify osmPacific.p8b \
            -pointNemo="-49.002579500,-123.391860387", -distance=2701065.845

   Coordinates of the Point Nemo (φ, λ) are in decimal degrees and Nemo
   distance is measured in meters, as the length of geodesic on ellipsoid.

   Programmer: Hrvoje Lukatela, 2023.
 */

#define PGM_DSCR "Point Nemo Disqualification"
#define PGM_LAST_EDIT_DATE "2025.093"         /* format as from 'date +%Y.%j' */

#include <nemo.h>
#include "../scullions/scullions.h" /* include after nemo.h has been included */

#define MAX_COORD_STR   64
#define DIST_EPSILON     0.025                              /* 25 millimetres */

static const char *progName;    /* for error logging by this source file only */
void usage(const char *, const char *);
/* ========================================================================== */
int main (int argc,
          const char *argv[],
          const char *envr[]) {

   int n, nIn, nOut;
   const char *fnIn;
   FILE *fpIn;
   nemoPtEll ptEll;           /* command linee input angular φ, λ coordinates */
   nemoPtEnr ptNemo;                                    /* claimet Point Nemo */
   double nemoDist;                        /* claimed Nemo distance, geodesic */
   const char delimiters[] = ", ";
   char coordStr[MAX_COORD_STR + 2];
   char *token;
   const char *optKey, *optVal;            /* options, in -keyword=value form */
   const char *strPtNemo, *strDistance;
   nemoPtUs8 ptUs8;                     /* Coastal point from file as CDC/U64 */
   nemoPtNcs ptNcs;                     /* as above, on near-conformal sphere */
   nemoPtEnr ptCoast;                        /* as above, as ellipsoid normal */
   nemoPtEll llCoast;                         /* as above, latitude/longitude */
   double g;
/* -------------------------------------------------------------------------- */
   progName = strrchr(argv[0], '/');                                 /* POSIX */
   if (progName == NULL) progName = strrchr(argv[0], '\\');         /* MS Win */
   if (progName) progName += 1;             /* skip found last path separator */
   if (progName == NULL) progName = argv[0];     /* neither? Just use argv[0] */
   fprintf(stderr, "\n[%s]: %s\nsource as of: %s, Nemo Library: %.3f\n",
                   progName, PGM_DSCR, PGM_LAST_EDIT_DATE, NEMO_LIBRARY_DATE);

   strPtNemo = strDistance = NULL;                        /* mamdatory values */
   while ((optKey = clOption(argc, argv, &optVal))) {
      if (*optKey == 'h') usage(NULL, NULL);
      else if (*optKey == 'p') strPtNemo = optVal;
      else if (*optKey == 'd') strDistance = optVal;
      else usage("unrecognized option", optKey);
      }

/* Extract claimed Point Nemo coordinates  (ptEll is kept, see below) */
   if (strPtNemo == NULL) usage("missing argument:", "Point Nemo coordinates");
   strncpy(coordStr, strPtNemo, MAX_COORD_STR);
   token = strtok(coordStr, delimiters);
   ptEll.a[NEMO_LAT] = NEMO_DEG2RAD * strtod(token, NULL);
   token = strtok(NULL, delimiters);
   ptEll.a[NEMO_LNG] = NEMO_DEG2RAD * strtod(token, NULL);
/* Convert Point Nemo coordinates to ellipsoid normal */
   nemo_LatLongToDcos3(ptEll.a, ptNemo.dc);

/* Get Nemo distance, geodesic meters on the surface */
   if (strDistance == NULL) usage("missing argument:", "Nemo distance");
   nemoDist = strtod(strDistance, NULL);

/* First and only file argument: input file path/name */
   fnIn = clFileName(argc, argv);                               /* input file */
   if (fnIn == NULL) usage("Missing input file name", NULL);
   fpIn = fopen(fnIn, "rb");                                 /* Open the file */
   if (fpIn == NULL) errorExit(progName, __LINE__,
                                "Can't open [%s] for reading\n", fnIn);

/* Report what was specified on the command line. This is quite necessary when
   the error checking is minimal or non-existant (in this case done to keep
   Library API demo code as succinct as possible) */
   fprintf(stderr, "Input file: %s\n", fnIn);
   fprintf(stderr, "Claimed Point Nemo:      %13.9f,%14.9f\n",
            NEMO_RAD2DEG * ptEll.a[NEMO_LAT], NEMO_RAD2DEG * ptEll.a[NEMO_LNG]);
   fprintf(stderr, "Claimed Nemo Distance: %13.3f\n", nemoDist);
   nIn = nOut = 0;
   n = fread(&ptUs8, sizeof(nemoPtUs8), 1, fpIn);
   while (n) {
      if (NEMO_Us8Plate(ptUs8) == 0) {    /* ignore possible ring-end markers */
         n = fread(&ptUs8, sizeof(nemoPtUs8), 1, fpIn);
         continue;
         }
      nIn++;
      nemo_Us8ToNcs(ptUs8, &ptNcs);
      nemo_NcsToEnr(nemo_ElrWgs84(), &ptNcs, &ptCoast);
      g = nemo_GeodesicSzpila(nemo_ElrWgs84(), &ptNemo,  &ptCoast, NULL);
      if (g == NEMO_DOUBLE_UNDEF) errorExit(progName, __LINE__,
                                            "Unexpected Vincenty failure\n");
      if (g < (nemoDist + DIST_EPSILON)) {   /* point is within Nemo distance */
         nemo_Dcos3ToLatLong(ptCoast.dc, llCoast.a);

         printf("%13.9f,%14.9f %6.3f\n", NEMO_RAD2DEG * llCoast.a[NEMO_LAT],
                                         NEMO_RAD2DEG * llCoast.a[NEMO_LNG],
                                         g - nemoDist);
         nOut++;
         }
      n = fread(&ptUs8, sizeof(nemoPtUs8), 1, fpIn);
      }
   fclose(fpIn);
   fprintf(stderr, "Points read: %8d, written: %8d\n", nIn, nOut);
   if (nOut > 3) return(1);
   else if (nOut < 3) return(-1);
   else return(0);
   }
/* ========================================================================== */
void usage(const char *mA,                  /* first message string (or NULL) */
           const char *mB) {               /* second message string (or NULL) */
   if (mA || mB) fprintf (stderr, "Error: %s %s\n", mA ? mA : "\0", mB ? mB : "\0");
   fprintf (stderr, "Usage: %s [options] inFile\n", progName);
   fprintf (stderr, "  inFile:  .p8b binary coordinate input file\n");
   fprintf (stderr, "Options:\n");
   fprintf (stderr, " -h[elp|  to print this usage help and exit\n");
   fprintf (stderr, " -p[ointNemo]=\"φ,λ\" Point Nemo coordinates, in decimal degrees\n");
   fprintf (stderr, " -d[istance]=nnn Nemo distance, meters on planetary surface\n");
   exit(1);
   }
/* ========================================================================== */
#include "../scullions/errorExit.c"
#include "../scullions/nemoStrings.c"
#include "../scullions/clFileOpt.c"
/* ========================================================================== */
