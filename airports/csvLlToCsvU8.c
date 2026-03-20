/* csvLlToCsvU8.c Given .csv airport file (Code, φ, λ, Name) create an
   equivalent .csv file, in which the location is one UniSpherical coordinate
   insted φ, λ.

   Since both are text files, the program is invoked as a "filter":
   cvsLlToCsvU8 < xyzLl.cvs > xyzU8.csv

   φ, λ, csv file:
   ---------------
   AAA,-17.353,-145.510,"Anaa"
   AAB,-26.696,141.049,"Arrabury"
   AAC,31.079,33.837,"El Arish"
   ...

   equivalent U8 csv file:
   -----------------------
   AAA,da0e8f9abac6b0c5,"Anaa"
   AAB,c6cd60f647fcc985,"Arrabury"
   AAC,41bffb383d26fb8d,"El Arish"
   ...

 */

#include <stdio.h>
#include <nemo.h>

#define LINE_MAX      512
/* ========================================================================== */
int main (int argc,
          const char *argv[],
          const char *envr[]) {
   int n;
   char textLine[LINE_MAX + 2];
   char *lineRead;
   char *token;
   char *code;
   char *name;
   nemoPtEll locEll;                  /* point location as φ, λ in radians... */
   nemoPtNcs locNcs;                      /* ...on "near-conformal" sphere... */
   nemoPtUs8 locUs8;                /* ... and 8-byte UniSpherical coordinate */
/* -------------------------------------------------------------------------- */
   n = 0;
   lineRead = fgets(textLine, LINE_MAX, stdin);
   while (lineRead) {
      code = strtok(textLine, ",");                          /* get IATA code */
      token = strtok(NULL, ",");                           /* Latitude... */
      locEll.a[0] = NEMO_DEG2RAD * atof(token);
      token = strtok(NULL, ",");                           /* ...Longitude... */
      locEll.a[1] = NEMO_DEG2RAD * atof(token);
      nemo_EllToNcs(nemo_ElrWgs84(), &locEll, &locNcs);
      locUs8 = nemo_NcsToUs8(&locNcs); /* ...converted to 8-byte UniSpherical */
      name = strtok(NULL, ",\n");                                 /* get name */
      printf("%s,%016lx", code, locUs8);
      if (name) printf(",%s\n", name);
      else printf("\n");
      n++;
      lineRead = fgets(textLine, LINE_MAX, stdin);
      }
   fprintf(stderr, "done, lines:  %d\n", n);
   }
/* ========================================================================== */
