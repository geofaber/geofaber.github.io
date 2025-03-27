/* Read the lines of .csv text file that has at
   least two items on the line, φ and λ - and create a CreateP8b binary file from csv text file: Traverse  Us8 coordinate binary
   version of it. Before writing the array of UniSpherical coordinates, sort
   the array.

   For instance:

   csvToP8b w1904711.csv w1904711.p8b
 */

#define PGM_DSCR "From .csv (φ, λ) create (Us8 format) .p8b file"
#define PGM_LAST_EDIT_DATE "2025.085"

#include <stdio.h>
#include <nemo.h>
#include "../scullions/scullions.h" /* include after nemo.h has been included */

static int compU640(const void *, const void *);
#define LINE_MAX      256

static const char *progName;    /* for error logging by this source file only */
/* ========================================================================== */
int main (int argc,
          const char *argv[],
          const char *envr[]) {
   int n;
   int lCount;
   char textLine[LINE_MAX + 2];
   char *lineRead;
   char *token;
   nemoPtUs8 *locations;
   nemoPtEll locEll;
   nemoPtNcs locNcs;
   FILE *inFp;                                                  /* input .csv */
   FILE *outFp;                                     /* output "canonical" ptb */
/* -------------------------------------------------------------------------- */
   progName = strrchr(argv[0], '/');                                 /* POSIX */
   if (progName == NULL) progName = strrchr(argv[0], '\\');         /* MS Win */
   if (progName == NULL) progName = argv[0];                      /* neither? */
   else progName += 1;                        /* strip leading path separator */
   fprintf(stderr, "\n[%s]: %s\nsource as of: %s, Nemo Library: %.3f\n",
                   progName, PGM_DSCR, PGM_LAST_EDIT_DATE, NEMO_LIBRARY_DATE);

   if (argc < 3) errorExit(progName, __LINE__,
                 "usage: %s xyz.csv xyzCnc.ptb\n", progName);

   inFp = fopen(argv[1], "rt");                        /* Open input file */
   if (inFp == NULL) errorExit(progName, __LINE__,
                                  "Can't open [%s] for reading\n", argv[1]);

   lineRead = fgets(textLine, LINE_MAX, inFp);
   lCount = 0;                                                 /* count lines */
   while (lineRead) {
/*    if (lCount < 10) fprintf(stderr, "%d:%s", lCount, textLine); */
      lCount++;
      lineRead = fgets(textLine, LINE_MAX, inFp);
      }
   rewind(inFp);
   fprintf(stderr, "in .csv lines %d\n", lCount);

   locations = malloc(lCount * sizeof(nemoPtUs8));
   if (locations == NULL) errorExit(progName, __LINE__, "No memory?\n");

   for (n = 0; n < lCount; n++) {
      fgets(textLine, LINE_MAX, inFp);
      token = strtok(textLine, ",");
/*    token = strtok(NULL, ",\n"); */  /* w1904711 csv includes line number!? */
      locEll.a[0] = NEMO_DEG2RAD * atof(token);
      token = strtok(NULL, ",\n");
      locEll.a[1] = NEMO_DEG2RAD * atof(token);
      nemo_EllToNcs(nemo_ElrWgs84(), &locEll, &locNcs);
      locations[n] = nemo_NcsToUs8(&locNcs);
/*    fprintf(stderr, "%s\n", nemo_StrU64Coords(locations[n])); */
      }
   fclose(inFp);

   fprintf(stderr, "Sort start...");
   qsort(locations, lCount, sizeof(nemoPtUs8), compU640);
   fprintf(stderr, " ...end\n");

   outFp = fopen(argv[2], "wb");                        /* Open input file */
   if (outFp == NULL) errorExit(progName, __LINE__,
                                  "Can't open [%s] for reading\n", argv[2]);
   fwrite(locations, sizeof(nemoPtUs8), lCount, outFp);
   fclose(outFp);
   free(locations);
   fprintf(stderr, "%s done, locations:  %d\n", progName, n);

   return(0);
   }
/* ========================================================================== */
/* compare two blocks of memory starting with a 64-bit unsigned integer       */
static int compU640(const void *p1, const void *p2) {
   if (*(uint64_t *)p1 > *(uint64_t *)p2) return(1);
   if (*(uint64_t *)p1 < *(uint64_t *)p2) return(-1);
   return(0);
   };
/* ========================================================================== */
#include "../scullions/errorExit.c"
#include "../scullions/nemoStrings.c"
/* ========================================================================== */
