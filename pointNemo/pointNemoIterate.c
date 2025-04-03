/* pointNemoIterate.c: find the ellipsoid coordinates of a point equidistant
   from three given, distant ellipsoid points - for instance, three Point Nemo
   proximity vertices. Distances between points are the length of geodesic.

   Given point coordinates are read from a text file. (see below).
   The coordinates are assumed to be the first two comma- or blank-separated
   line items, ellipsoidal φ and λ (respectively) in decimal degrees.

   The program can be compiled simply as:
   gcc path-to-nemo/programs/iterateForNemo.c -lm -o iterateForNemo

   Once compiled, the program will most likely be invoked by redirecting a
   text-file with at least three lines of coordinates to standard input, e.g.:

   ./pointNemoIterate < proximityVertices.pts

   where "proximityVertices.pts" is an ascii text file consisting of
   at least three (in this example five) lines, for example:
# approximate point Nemo φ, λ: -49.020468146,-123.436312526
# three proximity vertices and their distances:
-73.1904914,-127.0394759,  2702770.510
-24.6889471,-124.7868065,  2703129.654
-27.2022152,-109.4535548,  2704912.663

  (Any item following φ and λ is ignored, as are any additional input lines)

  N.B.: three given points (proximity vertices) may be given in any order.
  However, the program assumes that Point Nemo is a point on the same
  "relative hemisphere" as the proximity vertices (i.e., it is the "near",
  instead of the "far" pole of the small circle defined by three vertices.
  (Behavior of Vincenty calculation for geodesics greater than spheroidal
  quadrant needs additional research if this assumption is invalid: as,
  for instance, if attempting to determine the point of the globe
  (sea ~or~ land) that is at maximum geodesic distance from any coastline
  point of the continent of Australia).

  No input validation is performed and no protection against unstable
  or deteriorated geometry (for instance, proximity vertex distances
  approaching the spheroid quadrant) is offered.

  Coordinates and distances are reported to an approximate equivalent of
  at least a millimeter, regardless of the number of decimal fraction
  digits provided.

  The solution is obtained by first determining the position as a centre
  of the small circle defined by proximity vertices on near-conformal sphere,
  then iterating for a more precise location that is equidistant to all three
  proximity vertices on the ellipsoid (c.f., nudgeNemo() preamble). In this
  implementation the convergence is very slow, but the convergence criterion
  is quite high, and the process requires only a few, easy to understand
  spherical trigonometry and vector algebra productions.

  Programmer: Hrvoje Lukatela, <www.lukatela.com/hrvoje> 2022.
 */
#define PGM_DSCR "Iterative trilateration for Point Nemo"
#define PGM_LAST_EDIT_DATE "2024.211"         /* format as from 'date +%Y.%j' */

#include <nemo.h>
#include "../scullions/scullions.h" /* include after nemo.h has been included */

#define MAX_STEPS   1024                  /* the solution failed to converge? */
#define MAX_DIFF  0.0005            /* iteration criterion: half a millimeter */

void findGeoDists(nemoPtEnr *);         /* to find proximities and their mean */
void nudgeNemo(nemoPtEnr *);          /* nudge point to a better equidistance */

static nemoPtEnr proxVtxEl[3];   /* invar. proximity vertices ellipsoid i,j,k */
static nemoPtNcs proxVtxNs[3];                  /* as above, on "Nemo Sphere" */
static double glDist[3];        /* iteration-steep variant geodesic distances */
static double glMean;                          /* as above, mean of all three */
static const char *progName;
#define IN_LINE_LENGTH 256

int main (int argc,
          const char *argv[],
          const char *envr[]) {

   int ipv, iDir, nIter;
   char *pa, *token;
   char textLine[IN_LINE_LENGTH + 2];
   nemoPtEnr pointNemo;                                            /* φ and λ */
   nemoPtNcs ncsAux;                 /* auxiliary spherical point coordinates */
   nemoPtEnr elrAux;                                /* as above, on ellipsoid */
   nemoPtEll inPtEll;                          /* proximity vertex from input */
   double diff, d;
   char strLat[NEMO_DEBUG_STRING_LENGTH + 2],
        strLng[NEMO_DEBUG_STRING_LENGTH + 2];
/* -------------------------------------------------------------------------- */
   progName = strrchr(argv[0], '/');                                 /* POSIX */
   if (progName == NULL) progName = strrchr(argv[0], '\\');         /* MS Win */
   if (progName) progName += 1;             /* skip found last path separator */
   if (progName == NULL) progName = argv[0];     /* neither? Just use argv[0] */
   fprintf(stderr, "\n[%s]: %s\nsource as of: %s, Nemo Library: %.3f\n",
                   progName, PGM_DSCR, PGM_LAST_EDIT_DATE, NEMO_LIBRARY_DATE);

/* get φ and λ of three proximity vertices: */

   ipv = 0;
   pa = fgets(textLine, IN_LINE_LENGTH, stdin);
   while (ipv < 3) {
      if (pa == NULL) errorExit(progName, __LINE__, "failed to read 3 proximity vertices\n");
      if ((*pa == '#') || (*pa == '\n')) {         /*  comment or blank line? */
         pa = fgets(textLine, IN_LINE_LENGTH, stdin);              /* skip it */
         continue;
         }
      token = strtok(textLine, " ,");
      inPtEll.a[NEMO_LAT] = NEMO_DEG2RAD * atof(token);
      token = strtok(NULL, " ,");
      inPtEll.a[NEMO_LNG] = NEMO_DEG2RAD * atof(token);
      nemo_LatLongToDcos3(inPtEll.a, proxVtxEl[ipv].dc);    /* to ell. normal */
      nemo_EnrToNcs(nemo_ElrWgs84(), proxVtxEl + ipv, proxVtxNs + ipv);
      ipv++;
      pa = fgets(textLine, IN_LINE_LENGTH, stdin);
      }

 for (ipv = 0; ipv < 3; ipv++) {
      nemo_NcsToEll(nemo_ElrWgs84(), proxVtxNs + ipv, &inPtEll);
      fprintf(stderr, "%13.9f, %14.9f\n",
       NEMO_RAD2DEG * inPtEll.a[0], NEMO_RAD2DEG * inPtEll.a[1]);
      }
/* Prepare the iteration process. First, initialize Point Nemo as the
   circumcentre of proximity vertices on the Nemo Sphere.
 */
   iDir = nemo_SphereCircumcenter(proxVtxNs + 0, proxVtxNs + 1, proxVtxNs + 2,
                                  &ncsAux);
/* fprintf(stderr, "direction indicator: %d\n", iDir); */
   if (iDir == -1) {              /* must reverse the order of given vertices */
      ncsAux = proxVtxNs[0]; proxVtxNs[0] = proxVtxNs[2]; proxVtxNs[2] = ncsAux;
      elrAux = proxVtxEl[0]; proxVtxEl[0] = proxVtxEl[2]; proxVtxEl[2] = elrAux;
      iDir = nemo_SphereCircumcenter(proxVtxNs + 0,
                                     proxVtxNs + 1, proxVtxNs + 2, &ncsAux);
      }
   if (iDir != 1) errorExit(progName, __LINE__,
                            "ill-defined geometry of proximity vertices\n");

/* transfer the preliminary location back to the ellipsoid: */
   nemo_NcsToEnr(nemo_ElrWgs84(), &ncsAux, &pointNemo);
   findGeoDists(&pointNemo);          /* preliminary geodesics and their mean */

   diff = NEMO_DOUBLE_HUGE;               /* initialize convergence criterion */
   nIter = 0;
   while (diff > MAX_DIFF) {                               /* start iteration */
      if (nIter++ > MAX_STEPS) errorExit(progName, __LINE__,
       "failed to converge in %d iterations\n", MAX_STEPS);
      nudgeNemo(&pointNemo);
      diff = 0.0;
      for (ipv = 0; ipv < 3; ipv++) {
         d = fabs(glDist[ipv] - glMean);
         if (d > diff) diff = d;
         }
/*    fprintf(stderr, "Iter: %2d, max diff: %.3f (%.3f %.3f %.3f)\n",
       nIter, diff, glDist[0] - glMean, glDist[1] - glMean, glDist[2] - glMean); */
      }

/* Report Point Nemo coordinates and mean geodesic length */
   fprintf(stdout, "# %s iterations: %d\n", progName, nIter);
   nemo_Dcos3ToLatLong(pointNemo.dc, inPtEll.a);
   strncpy(strLat, nemo_StrSexagesimal(NEMO_RAD2DEG * inPtEll.a[0]),
           NEMO_DEBUG_STRING_LENGTH);
   strncpy(strLng, nemo_StrSexagesimal(NEMO_RAD2DEG * inPtEll.a[1]),
           NEMO_DEBUG_STRING_LENGTH);

   printf("# Point Nemo φ, λ and distance:\n");
   printf("%13.9f, %14.9f, (%s, %s), %12.3f\n",
           NEMO_RAD2DEG * inPtEll.a[0], NEMO_RAD2DEG * inPtEll.a[1],
           strLat, strLng, glMean);

/* Report Vertices and length of geodesic to each */
   printf("# Proximity Vertices  φ, λ and distance:\n");
   for (ipv = 0; ipv < 3; ipv++) {
      nemo_Dcos3ToLatLong(proxVtxEl[ipv].dc, inPtEll.a);
      strncpy(strLat, nemo_StrSexagesimal(NEMO_RAD2DEG * inPtEll.a[0]),
              NEMO_DEBUG_STRING_LENGTH);
      strncpy(strLng, nemo_StrSexagesimal(NEMO_RAD2DEG * inPtEll.a[1]),
              NEMO_DEBUG_STRING_LENGTH);
      printf("%13.9f, %14.9f, (%s, %s), %12.3f\n",
              NEMO_RAD2DEG * inPtEll.a[0], NEMO_RAD2DEG * inPtEll.a[1],
              strLat, strLng, glDist[ipv]);
      }

   return(0);
   }
/* ========================================================================== */
/* Populate the array of three geodesic distances and their mean value */
void findGeoDists(nemoPtEnr *ptNemo) {
   int ipv;
   for (ipv = 0; ipv < 3; ipv++) {
      glDist[ipv] = nemo_GeodesicSzpila(nemo_ElrWgs84(), ptNemo, proxVtxEl + ipv,
                                        NULL);
      }
   glMean = (glDist[0] + glDist[1] + glDist[2]) / 3.0;
   return;
   }
/* ========================================================================== */
/* "Nudge" the last best position of Point Nemo toward or away from that
   distant proximity vertex, from which the difference between the geodesic
   to it and the mean of all three geodesics is the greatest.

   The process is carried out on Nemo Sphere, where vector algebra productions
   are much simpler than they would be on an ellipsoid of rotation. Once
   the process is completed, the position is returned to the ellipsoid: the
   natural data domain of the Point Nemo and its three proximity vertices.
 */
void nudgeNemo(nemoPtEnr *ptNemo) {                      /* given and updated */

   int ipv;                             /* outer, proximity vertex loop index */
   int idc;                             /* inner, direction cosine loop index */
   int imx;              /* proximity vertex with maximum distance difference */
   double d, diff;
   nemoPtNcs ncsAux;                 /* auxiliary spherical point coordinates */
   double localScale;
   double nudge;              /* amount of nudge, in meters on planet surface */
   double dirToPrxVx[3];                     /* direction to proximity vertex */
/* -------------------------------------------------------------------------- */
/* fprintf(stderr, "%s ItStep\n", nemo_StrEnrCoords(ptNemo)); */
   nemo_EnrToNcs(nemo_ElrWgs84(), ptNemo, &ncsAux);   /* Point Nemo on sphere */
   localScale = nemo_NcsElrScale(nemo_ElrWgs84(), &ncsAux);
   diff = 0;     /* first, find the vertex with greatest difference from mean */
   for (ipv = 0; ipv < 3; ipv++) {    /* do for each distant proximity vertex */
      d = glDist[ipv] - glMean;
      if (fabs(d) > fabs(diff)) {
         diff = d;
         nudge = d;                                /* record nudge difference */
         imx = ipv;         /* record index of vertex with maximum difference */
         }
      }
/* fprintf(stderr, "nudge from %d by %6.3f\n", imx, nudge); */
/* Find vector from Point Nemo towards the vertex it will be nudged to or from */
   for (idc = 0; idc < 3; idc++) dirToPrxVx[idc] = proxVtxNs[imx].dc[idc] - ncsAux.dc[idc];
   nemo_NormalizeV3(dirToPrxVx);

/* Nudge Nemo along that vector, by ~nudge~, scaled down to unit sphere */
   for (idc = 0; idc < 3; idc++) {
      ncsAux.dc[idc] = ncsAux.dc[idc] + (dirToPrxVx[idc] * nudge / localScale);
      }
   nemo_NormalizeV3(ncsAux.dc);

   nemo_NcsToEnr(nemo_ElrWgs84(), &ncsAux, ptNemo);   /* back to ellipsoid... */
   findGeoDists(ptNemo);  /* ...and update ellipsoid geodesics and their mean */

   return;
   }
/* ========================================================================== */
#include "../scullions/errorExit.c"
#include "../scullions/nemoStrings.c"
/* ========================================================================== */
