#include <iostream>
//#include <fstream> // ecrire-lire dans un fichier
#include <ctime>
#include <cmath>
#include <random>  // Mersen Twister
#include <set> 
#include "data.h"

using namespace std;


/********************************************************************/
/**************************    GENERATORS    ************************/
/********************************************************************/
uint64_t rand_op(uint64_t N);
void print_bits(uint64_t a, int n_bits);
bool mat64_invertable(uint64_t *M, int n);

/********************************************************************/
/***************************     MODELS    **************************/
/********************************************************************/
Model64 rand_Model64(int m, int n)//select randomly m operators, m must be <= N
//Each integer corresponds exactly to one operator phi_i
{ 
  set<uint64_t> model; 
  uint64_t N = 0;
  if (n==64) {  N = ~N ; }
  else {  uint64_t un = 1;   N = (un << n) - un;  }

  if (m >= 2*N/3)
    { 
    for (uint64_t i = 1; i <= N; i++) 
      {  model.insert(i);  }
    while (model.size() > m)
      { model.erase(rand_op(N));  }
    }
  else
    {
    while(model.size() < m)
      {  model.insert(rand_op(N));  }
    }

  Model64 M;
  M.n=n;
  M.m=m;
  M.phi = (uint64_t *)malloc(m*sizeof(uint64_t));

  set<uint64_t>::const_iterator it1(model.begin()), it2(model.end()); 
  uint64_t i=0;

  for(;it1 != it2; ++it1) 
    {   M.phi[i] = *it1;  i++;  }
  if ( i != m ) { cout << "Error function rand_model" << endl;  }

  return M;
}

void print_t_Model64(Model64 M)   // print M^t = the transpose of M
{
  //cout << "|M| = " << M.m << endl << endl;
  cout << "M^t = " << endl;

  for(uint64_t i=0; i<M.m; i++)
  {
    print_bits(M.phi[i], M.n);
  }
  cout << endl; 
}

/********************************************************************/
/****************************     GT64    ***************************/
/********************************************************************/
void randomize_bin_GT64(GT64 *G)   // put 2 identical random binary matrix in *sig and *sig_work  --> not GT!
{
  uint64_t a;
  for (int i=0; i<(*G).n; i++)   
    {   a = rand_op((*G).N);   (*G).sig[i] = a;   (*G).sig_work[i] = a;  }
}

void new_GT64(GT64 *G)    //search for a new GT
{
  randomize_bin_GT64(G);

  while (!mat64_invertable((*G).sig_work, (*G).n))
    {  randomize_bin_GT64(G);  }
}

GT64 init_GT64(unsigned int n)  //initialize one GT
{
  GT64 G;     G.n=n;

  G.sig = (uint64_t *)malloc(n * sizeof(uint64_t));
  G.sig_work = (uint64_t *)malloc(n * sizeof(uint64_t));

  uint64_t N=0;
  if (n<=64)
    {
    if (n == 64) {  G.N = ~N ; }  else {  uint64_t un = 1;   G.N = (un << n) - un;  } // = 2^m-1

    new_GT64(&G);
    }
  else {  cout << "Error init_GT64: n is too large (n > 64)." << endl;   }

  return G;
}

void print_t_GT64(GT64 G)   //print transpose of GT
{
  cout << "GT^t = " << endl;
  for (int i=0; i<G.n; i++)
    {   print_bits(G.sig[i], G.n);  }
}

void free_GT64(GT64 G)
{
  free(G.sig);    free(G.sig_work);
}

/********************************************************************/
/**************************   GT of MODELS   ************************/
/********************************************************************/
// *** Needed here: 
//        -  GT64 init_GT64(unsigned int n);
//        -  void new_GT64(GT64 *G);
// ******

void GT_Model64(Model64 *M, uint64_t *GT)
{
  // initial model: m monomials phi[0 to m-1]
  // GT: vec of n independent monomials
  uint64_t phi;
  uint64_t phi_prim = 0;
  uint64_t i;

  /*cout << endl << "GT = ";
  for(int i=0; i<(*M).n; i++)
      { cout << GT[i] << " ";  }
  cout << ":" << endl;
  */

  // for each phi is associated the new monomial phi_prim
  for(uint64_t j=0; j<(*M).m ; j++) 
    { 
      phi_prim = 0;
      phi = (*M).phi[j];
      //cout << phi << " --> ";

      //check each bit of phi
      for(i = 0;  i < (*M).n;  i++ )
        {
        if (phi & 1)          //if the i-th bit of phi == 1; i.e. phi[i] == 1
          {  phi_prim = phi_prim^GT[i];   //cout << i << ", "; 
          }
        phi = (phi >> 1);  // a/2
        }

      //cout << phi_prim << endl;
      (*M).phi[j] = phi_prim;
    }
}

void GT_Model64_with_print(Model64 *M, uint64_t *GT)
{
  // initial model: m monomials phi[0 to m-1]
  // GT: vec of n independent monomials
  uint64_t phi;
  uint64_t phi_prim = 0;
  uint64_t i;

  cout << endl << "GT = ";
  for(int i=0; i<(*M).n; i++)
      { cout << GT[i] << " ";  }
  cout << ":" << endl;

  // for each phi is associated the new monomial phi_prim
  for(uint64_t j=0; j<(*M).m ; j++) 
    { 
      phi_prim = 0;
      phi = (*M).phi[j];
      cout << phi << " --> ";

      //check each bit of phi
      for(i = 0;  i < (*M).n;  i++ )
        {
        if (phi & 1)          //if the i-th bit of phi == 1; i.e. phi[i] == 1
          {  phi_prim = phi_prim^GT[i];   //cout << i << ", "; 
          }
        phi = (phi >> 1);  // a/2
        }

      cout << phi_prim << endl;
      (*M).phi[j] = phi_prim;
    }
}


void Nit_randGT_Model64(Model64 *M, int Nit)
{
  clock_t temps = clock();
  int compt=0;

  GT64 G = init_GT64((*M).n); //initialize and sample one GT
  cout << "Time for performing Nit = " << Nit << " GT (using uint64) starting from the model:  \t" << endl;
  print_t_Model64(*M);

  for (int i=0; i<Nit-1; i++)
  {
    if( i%(Nit/10) == 0 )  {cout << "avance : " << compt << "0%" << endl; compt++;}
    //print_t_Model64(M);
    //print_t_GT64(G);  cout << endl;
    GT_Model64(M, G.sig);   //GT(M)
    new_GT64(&G);          //find a new GT
  }
  GT_Model64(M, G.sig);   //last GT(M)

    //temps de tourne du programme
  temps = clock()- temps;
  cout << "Time = " << (double) temps/CLOCKS_PER_SEC << endl << endl;
}

void Nit_randGT_Model64_with_print(Model64 *M, int Nit)
{
  clock_t temps = clock();
  int compt=0;

  GT64 G = init_GT64((*M).n); //initialize and sample one GT

  for (int i=0; i<Nit-1; i++)
  {
    if( i%(Nit/10) == 0 )  {cout << "avance : " << compt << "0%" << endl; compt++;}
    //print_t_Model64(M);
    //print_t_GT64(G);  cout << endl;
    GT_Model64_with_print(M, G.sig);   //GT(M)
    new_GT64(&G);          //find a new GT
  }
  GT_Model64_with_print(M, G.sig);   //last GT(M)

    //temps de tourne du programme
  temps = clock()- temps;
  cout << endl << "Time = " << (double) temps/CLOCKS_PER_SEC << endl;

  free_GT64(G);
}
