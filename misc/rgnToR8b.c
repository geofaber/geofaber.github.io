/* Program to convert a text coordinate file from ./rgn/.lns/.pts GGW
   (Galileo Geodetic Workbench) terrestrial point, line and region
   (0, 1 and 2 dimensional terrestrial surface object text files to the
   equivalent UniSpherical binary files. The input files consist of three
   types of lines:

   1) Coordinate text line, consisting of two blank or comma separated items,
   φ and λ. While GGW allows different angular strings, this implementation
   only deals with both angles specified in signed (+/-) decimal degrees, with
   southerly latitude and westerly longitude negative.

   Coordinate lines are written to the output file as UniSpherical 8-byte
   records.

   2) For line and region object files a line segment or a closed ring is
   terminated by a "marker line" consisting of '*' character in the line's
   first character position. While not required, some files may have two
   additional numeric items following the "marker": the line segment or ring
   (sub-object) unsigned integer identifier, and (possibly in round brackets,
   "()") the count of vertices in the segment of ring that precedes the marker.
   Items on the line are, like the coordinates, blank or comma separated.

   Point (.pts) file have no sub-object terminator lines, but all three
   (.pts, lns and .rgn) files have a "file terminating" terminator text line.

   Sub-object terminator lines are written to the output file as "undefined"
   (or "invalid") Us8 coordinate records). If two additional integers are
   given, they will follow the most significant nibble (i.e., digiNental
   plate 0) as two 32-bit unsigned integers. (The first of those must therefore
   not be greater than 2**28 - 1).

   3) Lines starting with ';' character or completely blank lines are comments
   or "text readibility separator lines. Such lines are ignored and create no
   output file binary records.

   An example of the "head" of input file might be:
---------------------------------------------------
 55.7254490   -4.9423700
 55.7253088   -4.9425472
 55.7251789   -4.9423196
 55.7250998   -4.9422000
 55.7252119   -4.9420023
 55.7255849   -4.9419220
 55.7254490   -4.9423700
* 000000 (7)
 55.7245830   -4.9405109
 55.7249287   -4.9402624
 55.7253580   -4.9400959
 55.7252220   -4.9405430
 55.7249591   -4.9409839
 55.7246960   -4.9414240
 55.7245139   -4.9414586
 55.7243997   -4.9411861
 55.7242284   -4.9410542
 55.7240710   -4.9404860
 55.7245830   -4.9405109
* 000001 (11)
...
   And "tail":
...
 44.7452334  -62.8017838
 44.7453134  -62.8018915
 44.7452035  -62.8022040
 44.7451244  -62.8022483
 44.7450776  -62.8021633
 44.7452334  -62.8017838
* 758474 (6)
 44.8383509  -62.6137192
 44.8381255  -62.6138090
 44.8380239  -62.6136583
 44.8380574  -62.6135872
 44.8381372  -62.6136243
 44.8381738  -62.6136811
 44.8382579  -62.6136707
 44.8383509  -62.6137192
* 758475 (8)
 44.8383207  -62.6099584
 44.8382821  -62.6101048
 44.8381861  -62.6101730
 44.8380670  -62.6101591
 44.8380572  -62.6099925
 44.8381279  -62.6098624
 44.8382131  -62.6098295
 44.8382866  -62.6098535
 44.8383207  -62.6099584
* 758476 (9)
*
-----------------------------------------------------
   Note that the file in the above example is that of the two-dimensional
   objects (i.e., "regions"), the segments are "rings" and their starting
   vertex is (by OSM convention?) repeated as the last one in the sequence.
   (This program will not interfere).

   The program is invoked with two (ordered) arguments - input and output file
   path/names; for instance:

   ./rgnToR8b osmLand.rgn osmLand.
 */

#define PGM_DSCR "Convert .rgn text to .r8b binary file"
#define PGM_LAST_EDIT_DATE "2025.091"         /* format as from 'date +%Y.%j' */

#include <nemo.h>
#include "../scullions/scullions.h" /* include after nemo.h has been included */

/* pending inclusion to nemo.h */
#define NEMO_Us8Plate(u8) ((int)((u8 & 0xf000000000000000) >> 60))

void usage(const char *, const char *);
static const char *progName;    /* for error logging by this source file only */
/* ========================================================================== */
int main (int argc,
          const char *argv[],
          const char *envr[]) {

   int n, nLnIn, nRecOut;                    /* count input/output file lines */
   int nComments, nMarks, nSegPts, nTotalPts, markLast;
   int nCountMismatch, nIdSequence, nRingOpen;              /* OSM violations */
   int iSeg, nSeg, prevSeg, iPlate;
   const char *optKey;
   const char *optVal;
   const char *fnIn;         /* input: OSM rgn text coastline coordinate file */
   const char *fnOut; /* output: OSM coastline coordinate file as .lnb binary */
   FILE *fpIn, *fpOut;
   int lineSize;                                        /* must be signed int */
   int lineBufSize = 0;
   char *lineBuf = NULL;      /* will be allocated in first call to getLine() */
   const char delimiters[] = ", \r\n";
   char *pa;
   char *token;
   int maxVert, minVert;
   int nPtIn, nSegIn;
   nemoPtEll ptEll;                    /* input file angular φ, λ coordinates */
   nemoPtNcs ptNcs;                                       /* as above, in NCS */
   nemoPtUs8 ptUs8;                              /* as above, in UniSpherical */
   nemoPtUs8 ringStartPtUs8;                          /* as above, ring start */
   nemoPtUs8 segEndMark;                                  /* segment end mark */
   int mSeg, mSegPrev;                           /* input file segment number */
/* -------------------------------------------------------------------------- */
   progName = strrchr(argv[0], '/');                                 /* POSIX */
   if (progName == NULL) progName = strrchr(argv[0], '\\');         /* MS Win */
   if (progName) progName += 1;             /* skip found last path separator */
   if (progName == NULL) progName = argv[0];     /* neither? Just use argv[0] */
   fprintf(stderr, "\n[%s]: %s\nsource as of: %s, Nemo Library: %.3f\n",
                   progName, PGM_DSCR, PGM_LAST_EDIT_DATE, NEMO_LIBRARY_DATE);

   while ((optKey = clOption(argc, argv, &optVal))) {
      if (*optKey == 'h') usage(NULL, NULL);
      else usage("unrecognized option", optKey);
      }

   if (argc < 3) usage("missing command line filename arguments", NULL);

/* First file argument: input file path/name */
   fnIn = clFileName(argc, argv);                               /* input file */
   fpIn = fopen(fnIn, "rt");
   if (fpIn == NULL) errorExit(progName, __LINE__,
                               "Can't open [%s] for reading\n", argv[1]);
/* fprintf(stderr, "Reading from [%s]\n", argv[1]); */

/* Second file argument: output file path/name */
   fnOut = clFileName(argc, argv);                             /* output file */
   fpOut = fopen(fnOut, "wb");
   if (fpOut == NULL) errorExit(progName, __LINE__,
                                "Can't open [%s] for writing\n", argv[2]);
/* fprintf(stderr, "Writing to [%s]\n", argv[2]); */
   maxVert = 0;
   minVert = 2147483647;
   prevSeg = -1;
   nCountMismatch = nIdSequence = nRingOpen = 0;
   nLnIn = nRecOut = nComments = nMarks = nSegPts = nTotalPts = markLast = 0;
   lineSize = getLine(&lineBuf, &lineBufSize, fpIn);        /* get first line */
   while (lineSize > 0) {
      if (nLnIn%1000000 == 0) fprintf(stderr, "%d M\r", nLnIn/1000000);
      nLnIn++;
//    if (nLnIn > 1000) break;                                   /* debug...! */
      pa = lineBuf;
      while (*pa == ' ') pa++;                    /* skip over leading blanks */
      if ((*pa == '\0') || (*pa == ';') || (*pa == '#')) {    /* comment line */
         nComments++;
         lineSize = getLine(&lineBuf, &lineBufSize, fpIn);
         continue;
         }
      if (*pa == '*') {                                             /* marker */
         nMarks++;
         if (nSegPts) {            /* segment/ring vertices have been written */
            token = strtok(pa, delimiters);                    /* leading '*' */
/*          fprintf(stderr, "%s", token); */
            token = strtok(NULL, delimiters);              /* ring/segment id */
            iSeg = atoi(token);
/*          fprintf(stderr, " %u", iSeg); */
            token = strtok(NULL, delimiters);        /* vertex count, as (n)? */
            if (*token == '(') token++;
            nSeg = atoi(token);
/*          fprintf(stderr, " %u\n", nSeg); */
/*          Optional: count/report OSM convention violations: */
            if (nSeg != nSegPts) {                  /* vertex count mis-match */
               nCountMismatch++;
/*             fprintf(stderr,
               "Vertex count mis-match, input line: %d, counted: %d, given: %d\n",
                nLnIn, nSegPts, nSeg); */
               }
            if (iSeg != (prevSeg + 1)) { /* uninterrupted monotonic numbering */
               nIdSequence++;
/*             fprintf(stderr,
               "Segment numbering, input line: %d, previous: %d, current: %d\n",
                nLnIn, prevSeg, iSeg); */
               }
            if (ptUs8 != ringStartPtUs8) {                      /* open ring? */
               nRingOpen++;
/*             fprintf(stderr, "Open ring? line: %d, ring: %d, vertices: %d gap: %s\n",
               nLnIn, iSeg, nSegPts, nemo_StrUs8Dist(ptUs8, ringStartPtUs8)); */
               }
            prevSeg = iSeg;

            if (nSegPts > maxVert) maxVert = nSegPts;
            if (nSegPts < minVert) minVert = nSegPts;
            nTotalPts += nSegPts;
/*          Construct, write binary marker record: */
            segEndMark = 0L;
            segEndMark = iSeg;
            segEndMark = (segEndMark << 32) + nSegPts;
            iPlate = NEMO_Us8Plate(segEndMark);
            if (iPlate != 0) errorExit(progName, __LINE__,
              "ring/segment-id overflow? line: %d, ring-id: %d\n", nLnIn, iSeg);
            n = fwrite(&segEndMark, sizeof(nemoPtUs8), 1, fpOut);   /* marker */
            if (n != 1) errorExit(progName, __LINE__,
                    "Marker write error, line in, out: %d,%d\n", nLnIn, nRecOut);
            nRecOut++;
            nSegPts = 0;            /* ...and reset segment/ring vertex count */
            }
         else markLast++;         /* it better be the only one at file's end! */

         lineSize = getLine(&lineBuf, &lineBufSize, fpIn);
         continue;
         }

/*    Thus it must be a vertex in the line segment or ring */
      token = strtok(pa, delimiters);                                    /* φ */
      ptEll.a[NEMO_LAT] = NEMO_DEG2RAD * strtod(token, NULL);
      token = strtok(NULL, delimiters);                                  /* λ */
      ptEll.a[NEMO_LNG] = NEMO_DEG2RAD * strtod(token, NULL);
/*    fprintf(stderr, "%+11.7f,%+12.7f\n", NEMO_RAD2DEG * ptEll.a[NEMO_LAT],
                                           NEMO_RAD2DEG * ptEll.a[NEMO_LNG]); */
      nemo_EllToNcs(nemo_ElrWgs84(), &ptEll, &ptNcs);
      ptUs8 = nemo_NcsToUs8(&ptNcs);
      n = fwrite(&ptUs8, sizeof(nemoPtUs8), 1, fpOut);   /* write coordinates */
      if (n != 1) errorExit(progName, __LINE__,
              "Coordinates write error, line in, out: %d,%d\n", nLnIn, nRecOut);
      nRecOut++;
      nSegPts++;
      if (nSegPts == 0) ringStartPtUs8 = ptUs8;

      lineSize = getLine(&lineBuf, &lineBufSize, fpIn);    /* next input line */
      }

   fclose(fpIn);
   fclose(fpOut);
   free(lineBuf);

   fprintf(stderr, "Input file lines:            %8d\n", nLnIn);
   fprintf(stderr, "   comments:                 %8d\n", nComments);
   fprintf(stderr, "   line segment/rings:       %8d\n", nMarks);
   fprintf(stderr, "   minimum vertices/seg:     %8d\n", minVert);
   fprintf(stderr, "   maximum vertices/seg:     %8d\n", maxVert);
   fprintf(stderr, "   segment/ring vertices:    %8d\n", nTotalPts);
   fprintf(stderr, "OSM violations, count:       %8d\n", nCountMismatch);
   fprintf(stderr, "          id-sequence:       %8d\n", nIdSequence);
   fprintf(stderr, "           open rings:       %8d\n", nRingOpen);
   fprintf(stderr, "Output file records:         %8d\n", nRecOut);

   if (markLast > 1) errorExit(progName, __LINE__,
                                       "Terminating markers: %d?\n", markLast);
   return(0);
   }
/* ========================================================================== */
void usage(const char *mA,                  /* first message string (or NULL) */
           const char *mB) {               /* second message string (or NULL) */
   if (mA || mB) fprintf (stderr, "Error: %s %s\n", mA ? mA : "\0", mB ? mB : "\0");
   fprintf (stderr, "Usage: %s [option] inFile outFile\n", progName);
   fprintf (stderr, "  inFile:  .csv coordinate input file\n");
   fprintf (stderr, "  outFile: .lnb coordinate output file\n");
   fprintf (stderr, "Option:\n");
   fprintf (stderr, " -h(elp)      to print this usage help and exit\n");
   exit(1);
   }
/* ========================================================================== */
#include "../scullions/errorExit.c"
#include "../scullions/clFileOpt.c"
#include "../scullions/getLine.c"
#include "../scullions/nemoStrings.c"
/* ========================================================================== */
