#include <iostream>
#include <ctime>
#include <list>
#include <cmath>

using namespace std;

bool rand_bin()  { return (drand48()<0.5); }

/********************************************************************/
/***************************    MATRICES   **************************/
/********************************************************************/
bool** matrice(int n, int m)
{
  bool** M = (bool**) malloc(n*sizeof(bool*));

  for (int i=0; i<n; i++)
    { M[i] = (bool*) malloc(m*sizeof(bool));
      for (int j=0; j<m; j++)
          { M[i][j] = rand_bin(); }
    }

  return M;
}

void randomize_matrice(bool** M, int n, int m)
{
  for (int i=0; i<n; i++)
    {
      for (int j=0; j<m; j++)
          { M[i][j] = rand_bin(); }
    }
}

void print_matrice(bool** M, int n, int m)
{
  for (int i=0; i<n; i++)
    { 
      for (int j=0; j<m; j++)  {  cout << M[i][j] << " "; }
        cout << endl;
    }
}
/********************************************************************/
/***************************    OPERATIONS   ************************/
/********************************************************************/
void swap_row(bool** M, int i1, int i2, int n, int m)   //swap L_i1 and L_i2  in the matrix M
{
  if (i1>=n || i2>=n) { cout << "error swap" << endl; }
  else {
    bool* temp = M[i1]; //(bool*) malloc(m*sizeof(bool));
    M[i1]=M[i2];  M[i2]=temp;
    }
}

void add_row(bool** M, int i1, int i2, int n, int m)   //L_i2 <---- L_i1 XOR L_i2
{
  if (i1>=n || i2>=n) { cout << "error add" << endl; }
  else {
    for (int k = 0; k < m; k++)   {   M[i2][k] = ( M[i2][k] != M[i1][k] );   }
    }
}

/********************************************************************/
/********    INVERTIBILITY  --  for square matrices only  ***********/
/********************************************************************/
// This algo check if there is a zero "lead":
bool invertable(bool** M, int n)    // (i_lead, j_lead) = positions of the lead
{
  int j_lead = 0;
  bool test_stop=false;

  int i=0;

  for (int i_lead = 0; i_lead < n  &&  j_lead < n; i_lead++)   // 
    {
      //cout << "j_lead = " << j_lead << endl;
    i = i_lead;

    //search for the 1rst non-zero element in the column
    while (! M[i][j_lead])  // while M[i][j_lead] = 0
      {
        //cout << "test: " << M[i][j_lead] << endl;
      i++;
      if (i == n)   //all elements are 0  --> no lead in this column, 1.go to the next column and 2.restart
        {
        //cout << "no real j_lead = " << j_lead << "--> not invertable" << endl;
        test_stop=true; break;
        }
      }
    if(test_stop)   {  break;  }
    //cout << "real j_lead =" << j_lead << endl;
    //has found a non-zero element --->  1. increase rank; and 2. swap row i and row i_lead

    if (i != i_lead)  { swap_row(M, i_lead, i, n, n); }
      //print_matrice(M, n, m);  cout << endl;


    //put to zero every element under the lead (starting from i = i_lead + 1)
    i++;    
    while(i<n)
      {
      //cout << "i=" << i << endl;
      if (M[i][j_lead]) 
        { add_row(M, i_lead, i, n, n);   //L_i <---- L_i XOR L_i_lead
          //print_matrice(M, n, n);   cout << endl; 
        }  
      i++;      
      }

    //go to the next column
    j_lead++;
    }

    return (!test_stop);   // if test_stop==true  --> then M is not invertable  (rank < n)
}

double proba_invert(int n, int Nit=1e2)
{
  clock_t temps = clock();
  bool **M = matrice(n, n);

  int p=0;
  int compt=0;

  for (int i=0; i<Nit; i++)
  {
    if(i%(Nit/10)==0) {cout << "avance : " << compt << "0%" << endl; compt++;}
    if(invertable(M, n)) { p++;}
    randomize_matrice(M, n, n);
  }
  cout << "p = " << p << endl; 

  //temps de tourne du programme
    temps = clock()- temps;
    cout << "time generic n:  \t" << (double) temps/CLOCKS_PER_SEC << endl;

  return ((double) p)/Nit;
}

double sampling_GTF2(int n, int Nit=1e2)  //sample Nit GT
{
  clock_t temps = clock();
  bool **M = matrice(n, n);

  int i_success=0, i_tot=0;
  int compt=0;

  while (i_success<Nit)
  //for (int i=0; i<Nit; i++)
  {
    if(invertable(M, n)) 
      {
      if(i_success%(Nit/10)==0) {cout << "avance : " << compt << "0%" << endl; compt++;}
      i_success++; 
      }
    randomize_matrice(M, n, n);
    i_tot++;
  }
  cout << "p = " << ((double) i_success)/i_tot << endl; 

  //temps de tourne du programme
    temps = clock()- temps;
    cout << "time sampling GTF2 (for generic n):  \t" << (double) temps/CLOCKS_PER_SEC << endl;

  return ((double) i_success)/i_tot;
}

/********************************************************************/
/***************************   matrice REF   ************************/
/********************************************************************/
int REF(bool** M, int n, int m)    // (i_lead, j_lead) = positions of the lead
{
  int j_lead = 0;
  int rank=0;    //rank = final number of leads
  bool test_stop=false;

  int i=0;

  for (int i_lead = 0; i_lead < n  &&  j_lead < m; i_lead++)   // 
    {
      //cout << "j_lead = " << j_lead << endl;
    i = i_lead;

    //search for the 1rst non-zero element in the column
    while (! M[i][j_lead])  // while M[i][j_lead] = 0
      {
        //cout << "test: " << M[i][j_lead] << endl;
      i++;
      if (i == n)   //all elements are 0  --> no lead in this column, 1.go to the next column and 2.restart
        {
        cout << "no real j_lead = " << j_lead << endl;
        j_lead++;   //1. go to the next column
          if (j_lead >= m)  { //cout << "break j_lead = " << j_lead << endl; 
                            test_stop=true; break; }  // reduction finished
        i = i_lead; //2. re-start the search for non-zero element from new (i_lead, j_lead)
        }
      }
    if(test_stop)   {  break;  }
    //cout << "real j_lead =" << j_lead << endl;
    //has found a non-zero element --->  1. increase rank; and 2. swap row i and row i_lead
    rank++; 
    if (i != i_lead)  { swap_row(M, i_lead, i, n, m); }
      //print_matrice(M, n, m);  cout << endl;


    //put to zero every element under the lead (starting from i = i_lead + 1)
    i++;
    while(i<n)
      {
      //cout << "i=" << i << endl;
      if (M[i][j_lead]) 
        { add_row(M, i_lead, i, n, m);   //L_i <---- L_i XOR L_i_lead
          //print_matrice(M, n, m);   cout << endl; 
        }  
      i++;      
      }

    //go to the next column
    j_lead++;
    }

    //cout << "Final rank = " << rank << endl; 
    return rank;
}

/********************************************************************/
/**************************   matrice RREF   ************************/
/********************************************************************/
list<unsigned int> RREF(bool** M, int n, int m)    // (i_lead, j_lead) = positions of the lead
{
  int j_lead = 0;
  int rank=0;    //rank = final number of leads
  bool test_stop=false;

  list<unsigned int> no_lead_positions; 
  //int *lead_positions = (int *)malloc( min(n,m) * sizeof(int) );

  int i=0;

  for (int i_lead = 0; i_lead < n  &&  j_lead < m; i_lead++)   // 
    {
      //cout << "j_lead = " << j_lead << endl;
    i = i_lead;

    //search for the 1rst non-zero element in the column
    while (! M[i][j_lead])  // while M[i][j_lead] = 0
      {
        //cout << "test: " << M[i][j_lead] << endl;
      i++;
      if (i == n)   //all elements are 0  --> no lead in this column, 1.go to the next column and 2.restart
        {
        no_lead_positions.push_back(j_lead);
        j_lead++;   //1. go to the next column
          if (j_lead >= m)  { //cout << "break j_lead = " << j_lead << endl; 
                            test_stop=true; break; }  // reduction finished
        i = i_lead; //2. re-start the search for non-zero element from new (i_lead, j_lead)
        }
      }
    if(test_stop)   {  break;  }
    //cout << "real j_lead =" << j_lead << endl;
    //has found a non-zero element --->  1. increase rank; and 2. swap row i and row i_lead
    rank++; 
    if (i != i_lead)  { swap_row(M, i_lead, i, n, m); }
      //print_matrice(M, n, m);  cout << endl;

    //put to zero every element under the lead (starting from i = i_lead + 1)
    for (i=0; i < n; i++)
      {
        if (M[i][j_lead] && i!=i_lead)
          { add_row(M, i_lead, i, n, m);   //L_i <---- L_i XOR L_i_lead
            //print_matrice(M, n, m);   cout << endl; 
          }

      }

    //go to the next column
    j_lead++;
    }

  for (i=j_lead; i<m; i++)
    {
      no_lead_positions.push_back(i);  //cout << i << " ";
    }

  cout << "Final rank = " << rank << endl; 
  return no_lead_positions;
}

void print_loops(bool** M, int n, list<unsigned int> no_lead_positions)
{
  int i=0;
  for (list<unsigned int>::iterator it = no_lead_positions.begin(); it != no_lead_positions.end(); it++)
  {
    cout << (*it) << ": \t";
    for (i=0; i<n; i++)
      { if(M[i][(*it)])  {  cout << "S" << i << " "; }   }
    cout << endl;
  }
  cout << endl;
}

/********************************************************************/
/***************************    TESTINGS    *************************/
/********************************************************************/
/*
int main()
{
  initialise_rand();

  //  *******************   REF step-by-step   *********************** 
    int n=5, m=5;

    bool**M = matrice(n, m);  print_matrice(M, n, m);   cout << endl;
 
    swap_row(M, 0, 1, m);   print_matrice(M, n, m);   cout << endl;

    add_row(M, 0, 1, m);    print_matrice(M, n, m);   cout << endl;

    cout << "Invertable?  " << invertable(M, n) << endl;   //print_matrice(M, n, m);   cout << endl;

    REF(M, n, m);     // print_matrice(M, n, m);   cout << endl;

    free(M);    
 
  //  ********************   LOOPS   ***********************
  int n=5, m=6;
  bool**M = matrice(n, m);  print_matrice(M, n, m);   cout << endl;

  list<unsigned int> loops = RREF(M, n, m);   print_matrice(M, n, m);

  print_loops(M, n, loops);

  free(M);

  //  ***********************   Proba invertible   ***********************
  double p = proba_invert(30, 1e6) ;
  cout << "proba invertable: " << p << endl;   //<< " +/- " << p/sqrt(1e4) 

  return 0;
}
*/
