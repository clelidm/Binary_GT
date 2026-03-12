//g++ main.cpp RREF_F2.cpp RREF_64bits.cpp Model_GT64.cpp
//time ./a.out
//
//#include <stdio.h>
//#include <cstdlib>
#include <iostream>
//#include <fstream>
//#include <cmath>
#include <list>
//#include <vector>
//#include <string>
//#include <cinttypes> //for type int64_t in C++11 -- not necessary. Necessary for calling INT64_MAX
#include <random>  // Mersen Twister
#include <set> 
//#include "data.h"
#include "library.h"


using namespace std;


/******************************************************************************/
/************************** MAIN **********************************************/
/******************************************************************************/
//Rem : ---> 2^30 ~ 10^9

int main()
{
    // INITIALISATION
    init_rand();
    
    //CONSTANTS
    const unsigned int n = 10;  // number of spins
    cout << "n=" << n << endl;
    unsigned int m = 4;  // number of operators   // m must be:   0 < m <= N


    if ( n <= 64)       //number of spins   n<=64 ---> using integers encoded on n bits
    {
      uint64_t N = 0;
      if (n==64) {  N = ~N ; }
      else {  uint64_t un = 1;   N = (un << n) - un;  }
      cout << "N = ";  print_bits(N, n); cout << endl; 

      // *** model
      if (m <= N)    // m must be:   0 < m <= N
        {
        Model64 M = rand_Model64(m, n);

        //ex. for a single GT with prints:
        cout << endl << "*********** Example of a GT on a model: " << endl << endl;
        print_t_Model64(M); 
        GT64 G = init_GT64(n);  print_t_GT64(G);
        GT_Model64(&M, G.sig);  print_t_Model64(M); 
        
        cout << endl << "*********** 10 runs with prints: " << endl << endl;
        Nit_randGT_Model64_with_print(&M, 10);

        cout << endl << "*********** 10^6 runs -- timing: " << endl << endl;
        Nit_randGT_Model64(&M, 1e6);      // run 1e6 GT successively on M
        }
    }
    else
    {
      // *****   method 1 : M = bool matrix  *****  // binary matrix code 
      // ***********************   REF step-by-step   ***********************
      /*bool**M = matrice(n, n);  print_matrice(M, n, n);   cout << endl; 

      cout << "Invertable?  " << invertable(M, n) << endl;   //print_matrice(M, n, m);   cout << endl;

      REF(M, n, n);     // print_matrice(M, n, m);   cout << endl;

      free(M);
      */
    }

//  *******************   SPEED TEST and PROBA INVERTABLE   ************/
//  ******* test: run Nit GT, print the proba to sample a GT matrix and the speed:
    cout << endl << "*********** probability to sample a GT over 10^5 runs -- timing: " << endl << endl;
    sampling_GT64(n, 1e5);     // GT64
    sampling_GTF2(n, 1e5);     // GTF2

  return 0;
}
