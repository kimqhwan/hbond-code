/*************************************************
 * chemnetworks.h                                *
 *                                               *
 * Author: Abdullah Ozkanlar                     *
 *         abdullah.ozkanlar@wsu.edu             *
 *                                               *
 * A. Clark Research Lab, Chemistry Department   *
 * Washington State University, Pullman/WA 99164 *
 *************************************************/


#ifndef CHEMNETWORKS_H
#define CHEMNETWORKS_H

/* string manipulation */

int findf(FILE *fd, int n, ...);

/* graph types */

void graph_ss(double *atmS1, int nd1, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD,
              int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, 
              int pbc, double xside, double yside, double zside, FILE *outputfGraphS1S1, FILE *outputfGeodS1S1);

void graph_sAsB(double *atmS1, double *atmS2, int nd1, int nd2, int nsolvent1, int nsolvent2, int nAtomS1, int nAtomS2,
                int s1s2hbdn, int *s12a, int *s12b, double *s12as12bBD, 
                int s1s2hban, int *s1s2v1, int *s1s2v2, int *s1s2v3, int *s1s2v4, int *s1s2v5, double *s1s2v6, 
                int pbc, double xside, double yside, double zside, FILE *outputfGraphS1S2, FILE *outputfGeodS1S2);

void graph_st(double *atmS1, double *atmT1, int nd1, int nd2, int nsolvent1, int nsolute1, int nAtomS1, int nAtomT1, int s1t1cutoffnum, int *s1t1a, int *s1t1b,
              double *s1t1cutoff, int pbc, double xside, double yside, double zside, FILE *outputfGraphS1T1, FILE *outputfGeodS1T1);

void search_pbc_gss(int boxid, double *atmS1, double *atmS1x, int nd1, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, int s1s1hban,
                    int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, FILE *outputfGraphS1S1, FILE *outputfGeodS1S1);

void search_pbc_gsAsB(int boxid, double *atmS1, double *atmS2x, int nd1, int nd2, int nsolvent1, int nsolvent2, int nAtomS1, int nAtomS2, int s1s2hbdn, 
                      int *s12a, int *s12b, double *s12as12bBD, int s1s2hban, int *s1s2v1, int *s1s2v2, int *s1s2v3, int *s1s2v4, int *s1s2v5, double *s1s2v6,
                      FILE *outputfGraphS1S2, FILE *outputfGeodS1S2);

void search_pbc_gst(int boxid, double *atmT1, double *atmS1x, int nd1, int nd2, int nsolvent1, int nsolute1, int nAtomS1, int nAtomT1, int s1t1cutoffnum,
                    int *s1t1a, int *s1t1b, double *s1t1cutoff, FILE *outputfGraphS1T1, FILE *outputfGeodS1T1);

/* water dipoles */

void graph_st_dip(double *atmS1, double *atmT1, int nd1, int nd2, int nsolvent1, int nsolute1, int nAtomS1, int nAtomT1, int s1t1cutoffnum, int *s1t1a, int *s1t1b,
                  double *s1t1cutoff, int pbc, double xside, double yside, double zside, int opos, int h1pos, int h2pos, 
                  FILE *outputfGraphS1T1, FILE *outputfGeodS1T1, FILE *outputfS1T1Dip);

void search_pbc_gst_dip(int boxid, double *atmT1, double *atmS1x, int nd1, int nd2, int nsolvent1, int nsolute1, int nAtomS1, int nAtomT1, int s1t1cutoffnum, 
                        int *s1t1a, int *s1t1b, double *s1t1cutoff, int opos, int h1pos, int h2pos, 
                        FILE *outputfGraphS1T1, FILE *outputfGeodS1T1, FILE *outputfS1T1Dip);

void graph_sAsB_dip(double *atmS1, double *atmS2, int nd1, int nd2, int nsolvent1, int nsolvent2, int nAtomS1, int nAtomS2,
                    int s1s2hbdn, int *s12a, int *s12b, double *s12as12bBD, 
                    int s1s2hban, int *s1s2v1, int *s1s2v2, int *s1s2v3, int *s1s2v4, int *s1s2v5, double *s1s2v6, 
                    int pbc, double xside, double yside, double zside, int opos, int h1pos, int h2pos, int watid, int solid, 
                    FILE *outputfGraphS1S2, FILE *outputfGeodS1S2, FILE *outputfS1S2Dip);

void search_pbc_gsAsB_dip(int boxid, double *atmS1, double *atmS2x, int nd1, int nd2, int nsolvent1, int nsolvent2, int nAtomS1, int nAtomS2, int s1s2hbdn,
                          int *s12a, int *s12b, double *s12as12bBD, int s1s2hban, int *s1s2v1, int *s1s2v2, int *s1s2v3, int *s1s2v4, int *s1s2v5, double *s1s2v6,
                          int opos, int h1pos, int h2pos, int watid, int solid,
                          FILE *outputfGraphS1S2, FILE *outputfGeodS1S2, FILE *outputfS1S2Dip);

double dipole_angle(double *atmT, double *atmS, int i, int j, int *stb, int *sta, int crt, int opos, int h1pos, int h2pos);

/* structures */

void hringsearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, 
                 int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                 FILE *outputfHRing, char *foutputHRingIso);

void hbooksearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD,
                 int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                 FILE *outputfHBook, char *foutputHBookIso);

void hprismsearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, 
                  int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                  FILE *outputfHPrism, char *foutputHPrismIso);

void hcagesearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, 
                 int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                 FILE *outputfHCage, char *foutputHCageIso);

void hbagsearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, 
                int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                FILE *outputfHBag, char *foutputHBagIso);

void hboatsearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, 
                 int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                 FILE *outputfHBoat, char *foutputHBoatIso);

void hchairsearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, 
                  int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                  FILE *outputfHChair, char *foutputHChairIso);

void hprismbooksearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, 
                      int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                      FILE *outputfHPrismBook, char *foutputHPrismBookIso);

int check_isolated_hexamer(double *atmS, int nAtomS, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, int s1s1hban,
                           int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6,
                           int atma, int atmb, int atmc, int atmd, int atme, int atmf);

void hpentamersearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, 
                     int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                     FILE *outputfPentamer, char *foutputPentamerIso);

int check_isolated_pentamer(double *atmS, int nAtomS, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, int s1s1hban, 
                            int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, 
                            int atma, int atmb, int atmc, int atmd, int atme);

void htetramersearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, 
                     int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                     FILE *outputfTetramer, char *foutputTetramerIso);

int check_isolated_tetramer(double *atmS, int nAtomS, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, int s1s1hban, 
                            int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, 
                            int atma, int atmb, int atmc, int atmd);

void htrimersearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, 
                   int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                   FILE *outputfTrimer, char *foutputTrimerIso);

int check_isolated_trimer(double *atmS, int nAtomS, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, int s1s1hban, 
                          int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, 
                          int atma, int atmb, int atmc);

int prismbook_check(double *atmSX, int ndX, int opos, int sXsXhbdn, int *sXa, int *sXb, double *sXasXbBD, 
                    int sXsXhban, int *sXsXv1, int *sXsXv2, int *sXsXv3, int *sXsXv4, int *sXsXv5, double *sXsXv6,
                    int atmaX, int atmbX, int atmcX, int atmdX, int atmeX, int atmfX);

int ring_check(double *atmSX, int ndX, int opos, int sXsXhbdn, int *sXa, int *sXb, double *sXasXbBD, 
               int sXsXhban, int *sXsXv1, int *sXsXv2, int *sXsXv3, int *sXsXv4, int *sXsXv5, double *sXsXv6,
               int atmaX, int atmbX, int atmcX, int atmdX, int atmeX, int atmfX);

/* distance */

double distance(double *atmS, int i, int j,int *sa, int *sb, int crt);

double distanceMix(double *atmS1, double *atmS2, int i, int j,int *sa, int *sb, int crt);

double distanceOxOy(double *atmS, int i, int j, int opos);

/* angle */

double angle(double dist, double dist2, double hyptns);

/* dihedral */

double dihedral(double *atmS, int atm1, int atm2, int atm3, int atm4, int opos);

double dihedral2(double *atmS, int atm1, int atm2, int atm3, int atm4, int opos);

/* permutations */

int perm6(int v6[], int n, int i, int v2[], int *e);

void swap6(int v6[], int i, int j);

int perm5(int v5[], int n, int i, int v2[], int *e);

void swap5(int v5[], int i, int j);

int perm4(int v4[], int n, int i, int v2[], int *e);

void swap4(int v4[], int i, int j);

int perm3(int v3[], int n, int i, int v2[], int *e);

void swap3(int v3[], int i, int j);


// geodesics

void geodesics_ss(int nsolvent, int nAtomS, int EucDistS, int EucRefS, double xside, double yside, double zside, char *foutputGeodSS, char *finput);

void geodesics_sAsB(int nsolventA, int nAtomSA, int nsolventB, int nAtomSB, int EucDistSASB, int EucRefSASBsa, int EucRefSASBsb, double xside, double yside, double zside, 
                    char *foutputGeodSASB, char *finputA, char *finputB);

void geodesics_sAsBsC(int nsolventA, int nAtomSA, int nsolventB, int nAtomSB, int nsolventC, int nAtomSC, int EucDistSASBSC, int EucRefSASBSCsa, int EucRefSASBSCsb, int EucRefSASBSCsc, 
                      double xside, double yside, double zside, char *foutputGeodSASBSC, char *finputA, char *finputB, char *finputC);

void GetPath(int i, int j, int **next, int **gdmatrix, FILE *fout);

struct edge{
    int nodea;
    int nodeb;
    int change[3];
};

int FindChange(int startnode, int endnode, struct edge* edgelist, int numedgelist);

#endif // CHEMNETWORKS_H
