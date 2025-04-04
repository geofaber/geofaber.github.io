/* Programs calculating elements of a large, easy to visualize ellipsoid
   triangle, using chord direct and inverse geodetic problem functions.
   (This program is tipically used to confirm compile-time inclusion of
   functions in "scullions"" collection).
 */

#define PGM_DSCR "Calculating ellipsoid triangle geometry"
#define PGM_LAST_EDIT_DATE "2025.090"         /* format as from 'date +%Y.%j' */

#include <stdio.h>
#include <nemo.h>
#include "../scullions/scullions.h" /* include after nemo.h has been included */

#define ICY3(i) (((i)+3)%3)             /* index of an array element modulo 3 */
static const char *progName;    /* for error logging by this source file only */
void usage(const char *, const char *);
/* ========================================================================== */
int main (int argc,
          const char *argv[],
          const char *envr[]) {
   int i;                                         /* index of triangle vertex */
   const char *optKey;
   const char *optVal;
   nemoPtEll ptEll;
   nemoPtEnr vrtxEnr[3], vtxRet;                           /* triangle vertex */
   char *name[3] = {"Zagreb", "Dublin", "Timbak"};       /* names of vertices */
   double chSq, chord[3], geodVin;                   /* triangle side lengths */
   nemoDxPln dxAB[3], dxPrev, dxNext;                           /* directions */
   double aPrev, aNext, aSum, intAng[3], eps; /* spherical excess calculation */
/* -------------------------------------------------------------------------- */
   progName = strrchr(argv[0], '/');                                 /* POSIX */
   if (progName == NULL) progName = strrchr(argv[0], '\\');         /* MS Win */
   if (progName == NULL) progName = argv[0];                      /* neither? */
   else progName += 1;                        /* strip leading path separator */
   fprintf(stderr, "\n[%s]: %s\nsource as of: %s, Nemo Library: %.3f\n",
                   progName, PGM_DSCR, PGM_LAST_EDIT_DATE, NEMO_LIBRARY_DATE);

   while ((optKey = clOption(argc, argv, &optVal))) {
      if (*optKey == 'h') usage(NULL, NULL);
      else usage("unrecognized option", optKey);
      }

   fprintf(stderr, "                φ°              λ°\n");
   ptEll.a[NEMO_LAT] = NEMO_DEG2RAD * 45.814565201;
   ptEll.a[NEMO_LNG] = NEMO_DEG2RAD * 15.979425507;
   nemo_LatLongToDcos3(ptEll.a, vrtxEnr[0].dc);
   fprintf(stderr, "%s:  %s\n", name[0], nemo_StrEnrCoords(vrtxEnr + 0));

   ptEll.a[NEMO_LAT] = NEMO_DEG2RAD * 53.339754879;
   ptEll.a[NEMO_LNG] = NEMO_DEG2RAD * -6.272038955;
   nemo_LatLongToDcos3(ptEll.a, vrtxEnr[1].dc);
   fprintf(stderr, "%s:  %s\n", name[1], nemo_StrEnrCoords(vrtxEnr + 1));

   ptEll.a[NEMO_LAT] = NEMO_DEG2RAD * 16.775833333;
   ptEll.a[NEMO_LNG] = NEMO_DEG2RAD * -3.009444444;

   nemo_LatLongToDcos3(ptEll.a, vrtxEnr[2].dc);
   fprintf(stderr, "%s:  %s\n", name[2], nemo_StrEnrCoords(vrtxEnr + 2));

   fprintf(stderr, "\n          ");
   for (i = 0; i < 3; i++) fprintf(stderr, "%s-%s  ", name[i], name[ICY3(i+1)]);

   fprintf(stderr, "\nchord:    ");
   for (i = 0; i < 3; i++) {
      chSq = nemo_EllipsoidChordInverse(nemo_ElrWgs84(),
                                        vrtxEnr + i, vrtxEnr + (ICY3(i+1)),
                                        NULL, NULL);
      chord[i] = sqrt(chSq);
      fprintf(stderr, "%12.3fm  ", chord[i]);
      }

   fprintf(stderr, "\ngeodesic: ");
   for (i = 0; i < 3; i++) {
      geodVin = nemo_GeodesicSzpila(nemo_ElrWgs84(),
                                       vrtxEnr + i, vrtxEnr + (ICY3(i+1)), NULL);
      fprintf(stderr, "%12.3fm  ", geodVin);
      }

   fprintf(stderr, "\nazimuth:    ");
   for (i = 0; i < 3; i++) {
      chSq = nemo_EllipsoidChordInverse(nemo_ElrWgs84(),
                                        vrtxEnr + i, vrtxEnr + (ICY3(i+1)),
                                        dxAB + i, NULL);
      fprintf(stderr, "%10.6f°    ", NEMO_RAD2DEG * nemo_DirectionToAzimuth(dxAB[i].dc));
      }

   fprintf(stderr, "\n\nDirect problem chord, end points:\n");
   for (i = 0; i < 3; i++) {
      nemo_EllipsoidChordDirect(nemo_ElrWgs84(),
                                vrtxEnr + (ICY3(i-1)), dxAB + (ICY3(i-1)),
                                chord[ICY3(i-1)], &vtxRet, 0.0001, NULL);
      fprintf(stderr, "         %s\n", nemo_StrEnrCoords(&vtxRet));
      }

/* Spherical excess: */
   aSum = 0.0;
   for (i = 0; i < 3; i++) {
      nemo_EllipsoidChordInverse(nemo_ElrWgs84(),
                                 vrtxEnr + i,
                                 vrtxEnr + (ICY3(i-1)),
                                 &dxPrev, NULL);
      nemo_EllipsoidChordInverse(nemo_ElrWgs84(),
                                 vrtxEnr + i,
                                 vrtxEnr + (ICY3(i+1)),
                                 &dxNext, NULL);
      aPrev = nemo_DirectionToAzimuth(dxPrev.dc);
      aNext = nemo_DirectionToAzimuth(dxNext.dc);
      intAng[i] = aNext - aPrev;
      if (intAng[i] < 0.0) intAng[i] += NEMO_TWOPI;
/*    fprintf(stderr, "at %s, to %s, a:%10.6f\n", name[i], name[ICY3(i-1)], aPrev);
      fprintf(stderr, "at %s, to %s, a:%10.6f\n", name[i], name[ICY3(i+1)], aNext);
      fprintf(stderr, "%13.9f\n", NEMO_RAD2DEG * intAng[i]); */
      aSum += intAng[i];
      }
   eps = aSum - NEMO_PI;
   fprintf(stderr, "\nε: %13.9f°\n", NEMO_RAD2DEG * eps);
   return(0);
   }
/* ========================================================================== */
void usage(const char *mA,                  /* first message string (or NULL) */
           const char *mB) {               /* second message string (or NULL) */
   if (mA || mB) fprintf (stderr, "Error: %s %s\n", mA ? mA : "\0", mB ? mB : "\0");
   fprintf (stderr, "Usage: %s [options]\n", progName);
   fprintf (stderr, "Option:\n");
   fprintf (stderr, " -h(elp)      to print this usage help and exit\n");
   exit(1);
   }
/* ========================================================================== */
#include "../scullions/clFileOpt.c"
#include "../scullions/nemoStrings.c"
/* ========================================================================== */
