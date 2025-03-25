/* listP8b: List .p8b file consisting of an array of Us8 records.
   Nemo Library API example, cf.: https://www.lukatela.com/nemoLibrary

Commonly run as:

   listP8b xyz.p8b n > xyz.pts

   First argument is the input .p8b (UniSpherical 8-byte binary records) file,
   second is an optional integer to limit the number of output text lines.
   Output lines consist of blank-separated φ, λ and a 16-digit unsigned
   hexadecimal number (for instance, "-21.2333 -45.0000 1038e9d52b9dcc56".
   The first hexadecimal digit (in the above example, "1" if the number of
   digiNental plate - thus anything other than digits 1-6 would represent an
   invalid coordinate. Input file is assumed to be ordered on the UniSpherical
   coordinate and contain no duplicate points. The program will report and
   abort if any of those condition is violated.

   Note that any code that reads or writes binary UniSpherical coordinates
   is endianness specific. By convention, binary coordinate filec in mixed
   hardware environments should be assumed to be of little-Endian variety.

   Output should be redirected if further processing is anticipated.
 */

#define PGM_DSCR "List coordinates in .p8b file"
#define PGM_LAST_EDIT_DATE "2025.075"
#include <stdio.h>
#include <nemo.h>
#include "../scullions/scullions.h" /* include after nemo.h has been included */

static const char *progName;    /* for error logging by this source file only */
/* ========================================================================== */
int main (int argc,
          const char *argv[],
          const char *envr[]) {
   int k, n, m;
   FILE *inFp;                       /* input binary file, locations to visit */
   nemoPtUs8 ptUs8, prevPtUs8;
   nemoPtEll locEll;
   nemoPtNcs locNcs;
   int iPlate;
/* -------------------------------------------------------------------------- */
   progName = strrchr(argv[0], '/');                                 /* POSIX */
   if (progName == NULL) progName = strrchr(argv[0], '\\');         /* MS Win */
   if (progName == NULL) progName = argv[0];                      /* neither? */
   else progName += 1;                        /* strip leading path separator */
   fprintf(stderr, "\n[%s]: %s\nsource as of: %s, Nemo Library: %.3f\n",
                   progName, PGM_DSCR, PGM_LAST_EDIT_DATE, NEMO_LIBRARY_DATE);

   if (argc < 2) errorExit(progName, __LINE__,
                 "usage: %s xyzName.p8b [n]\n", progName);

   inFp = fopen(argv[1], "rb");                            /* Open input file */
   if (inFp == NULL) errorExit(progName, __LINE__,
                                "Can't open [%s] for reading\n", argv[1]);

   if (argc > 2) k = atoi(argv[2]);           /* limit number of output lines */
   else k = 0;                                               /* list them all */

   m = fread(&ptUs8, sizeof(nemoPtUs8), 1, inFp);
/* fprintf(stderr, "first read: %d\n", m); */

   prevPtUs8 = n = 0;
   while (m) {
      if ((k) && (n >= k)) break;                             /* want no more */
      nemo_Us8ToNcs(ptUs8, &locNcs);
      if (ptUs8 == prevPtUs8) errorExit(progName, __LINE__,
                 "input line %d: duplicate coordinates.\n", n);
      if (ptUs8 < prevPtUs8) errorExit(progName, __LINE__,
                 "input line %d: file sort order?.\n", n);
      iPlate = (int)((ptUs8 & 0xf000000000000000) >> 60);
      if ((iPlate < 1) || (iPlate > 6)) errorExit(progName, __LINE__,
                 "input line %d: invalid digiNental plate number [%016lx].\n", n, ptUs8);

      nemo_NcsToEll(nemo_ElrWgs84(), &locNcs, &locEll);
      printf("%8.4f %9.4f %016lx\n", NEMO_RAD2DEG * locEll.a[0],
                                     NEMO_RAD2DEG * locEll.a[1], ptUs8);
      n++;
      m = fread(&ptUs8, sizeof(nemoPtUs8), 1, inFp);
      }
   fclose(inFp);
   fprintf(stderr, "%s done, locations:  %d\n", progName, n);
   return(0);
   }
/* ========================================================================== */
#include "../scullions/errorExit.c"
/* ========================================================================== */
