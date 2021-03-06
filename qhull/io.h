/* io.h -- declarations of Input/Output functions

   Author:   Hannu Huhdanpaa
   Internet: hannu@geom.umn.edu

   The Geometry Center

   07/09/93 - created (hh)
   */


void errprint(facetT *facet, ridgeT *ridge, vertexT *vertex, pointT *point);
void printfacet(FILE *fp, facetT *facet);
void printfacets(FILE *fp, facetT *facetlist);
void printdegeneratehelp(FILE *fp);
void printhashtable(FILE *fp);
void printincidences(FILE *fp, facetT *facetlist);
void printoff2(FILE *fp, facetT *facetlist);
void printoff3(FILE *fp, facetT *facetlist, pointT *points, int numpoints);
void printoff4(FILE *fp, facetT *facetlist, pointT *points, int numpoints);
void printpoint(FILE *fp, pointT *point);
void printridge(FILE *fp, ridgeT *ridge);
void printsummary(FILE *fp, facetT *facetlist);
void printtriangles(FILE *fp, facetT *facetlist);
void printvertex(FILE *fp, vertexT *vertex);
pointT *readpoints (int *numpoints, int *dimension);









































