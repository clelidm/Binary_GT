#include <iostream>
#include <list>
//#include <string>
#include "data.h"
using namespace std;

/********************************************************************/
/*******************   GENERATOR  and TOOLS   ***********************/
/********************************************************************/
void init_rand();
void print_bits(uint64_t a, int n);	  // a integers, n = number of bits

/******************************************************************************/
/******************   ModelF2 and GTF2    for n > 64   ************************/
/******************************************************************************/
bool** matrice(int n, int m);                       //create random bool Mat of n x m
void randomize_matrice(bool** M, int n, int m);     //randomize a bool Mat of n x m
void print_matrice(bool** M, int n, int m);         //print

bool invertable(bool** M, int n);                   //check if M is invertable: for square M
double proba_invert(int n, int Nit=1e2);            //  return the pba that a random binary matrice is invertable on Nit trials --- allow testing the speed
double sampling_GTF2(int n, int Nit=1e2);             // sample Nit GT - return the pba that a random bin mat is invertable

int REF(bool** M, int n, int m);                    //  REF M and return the rank
list<unsigned int> RREF(bool** M, int n, int m);    // RREF M and return a list of the non-lead positions (i.e. column index of the circuits)
void print_loops(bool** M, int n, list<unsigned int> no_lead_positions);  // print the circuits of the RREF M

/******************************************************************************/
/*******************   Model64 and GT64   for n<=64   *************************/
/******************************************************************************/
Model64 rand_Model64(int m, int n);   //select randomly m operators, m must be <= 2^n -1
void print_t_Model64(Model64 M);      //print the transpose of Model
void GT_Model64(Model64 *M, uint64_t *GT);

GT64 init_GT64(unsigned int n);
void new_GT64(GT64 *G);
void Nit_randGT_Model64(Model64 *M, int Nit);
void Nit_randGT_Model64_with_print(Model64 *M, int Nit);
void print_t_GT64(GT64 G);


/******************************************************************************/
/*************************    mat64   for  n<=64   ****************************/
/******************************************************************************/
//
// **** !!! n x m  binary matrix = mat64 of n lines, each line is encoded by a single integer on m bits
//
//  return the pba that a random binary matrice is invertable on Nit trials --- allow testing the speed
double sampling_GT64(int n, int Nit=1e2);             //sample Nit GT --> for testing

/*
uint64_t *init_mat64(unsigned int n, unsigned int m); //create random mat n x m =  n integers of m bits
void randomize_mat64(uint64_t * M, int n, int m);     //randomize mat n x m
void print_t_mat64(uint64_t *M, int n, int m);          //print mat n integers of m bits

bool mat64_invertable(uint64_t *M, int n);            //check if M is invertable: for square mat M --> n integers of n bits

int mat64_REF(uint64_t *M, int n, int m);             // REF M and return the rank
list<unsigned int> mat64_RREF(uint64_t *M, int n, int m);
*/

