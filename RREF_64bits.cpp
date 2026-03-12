#include <iostream>
//#include <fstream> // ecrire-lire dans un fichier
#include <ctime>
#include <list>
#include <cmath>
#include <bitset>
#include <random>  // Mersen Twister

using namespace std;

/********************************************************************/
/**************************    GENERATORS    ************************/
/********************************************************************/
std::mt19937_64 gen_64;   // Mersenne Twister

void init_rand()
{
  int seed_val = (unsigned)time(NULL);
  srand48(seed_val);
  gen_64.seed(seed_val);  
}

uint64_t rand_op(uint64_t N)   //return a random integer i \in [1, N],    N included , 0 excluded (not an operator)
// !!!! --> change the probability (a bit higher than 2.88)
{ 
  uniform_int_distribution<uint64_t> uni_int64(1, N);
  return uni_int64(gen_64);
}

/********************************************************************/
/************************    MATRICES mat64   ***********************/
/********************************************************************/
//each line i of the matrice is given by the binary representation of the integer M[i]:
uint64_t *init_mat64(unsigned int n, unsigned int m)
{
  uint64_t *M = (uint64_t *)malloc(n * sizeof(uint64_t));

  uint64_t N=0;
  if (m<=64)
    {
    if (m == 64) {  N = ~N ; }  else {  uint64_t un = 1;   N = (un << m) - un;  } // = 2^m-1

    for (int i=0; i<n; i++)
      {    M[i] = rand_op(N);   }
    }
  else {  cout << "Error init_mat64: m is too large (m > 64)." << endl;   }

  return M;
}

void randomize_mat64(uint64_t * M, int n, int m)
{
  uint64_t N=0;
  if (m == 64) {  N = ~N ; }  else {  uint64_t un = 1;   N = (un << m) - un;  } // = 2^m-1

  for (int i=0; i<n; i++)
  {    M[i] = rand_op(N);   }
}

/**********************    LOCAL functions   **************************/
void randomize_vec64(uint64_t * M, int n, uint64_t N)
{
  for (int i=0; i<n; i++)   
    {    M[i] = rand_op(N);   }
}

/********************************************************************/
/***************************    PRINTS    ***************************/
/********************************************************************/
void print_bits(uint64_t a, int m) // n is the number of bits encoding a --> n= number of bits to print
{
  cout << "\t ";
  uint64_t b = a;
  //cout << a << ": \t" << bitset<2>(a) << endl;

  for( int i = 0;  i < m;  i++ )
    {   cout << (a & 1)  << " ";   a = (a >> 1);  // a/2
    }
  cout << "\t :" << b ;
  cout << endl;
}

void print_t_mat64(uint64_t *M, int n, int m)   //print transpose of GT
{
  for (int i=0; i<n; i++)    {   print_bits(M[i], m);  }
}

void print_mat64(uint64_t *M, int n, int m)   // print m64 of n integers of m bits
{
  for (int i=0; i<n; i++)    {   print_bits(M[i], m);  }
}

/********************************************************************/
/***************************    OPERATIONS   ************************/
/********************************************************************/
void swap_row_mat(uint64_t *M, int i1, int i2, int n)   //swap L_i1 and L_i2  in the matrix M
{
  if (i1>=n || i2>=n) { cout << "error swap" << endl; }
  else {
    uint64_t temp = M[i1];
    M[i1]=M[i2];  M[i2]=temp;
    }
}

void add_row_mat(uint64_t *M, int i1, int i2, int n)   //L_i2 <---- L_i1 XOR L_i2
{
  if (i1>=n || i2>=n) { cout << "error add" << endl; }
  else {
    M[i2] = ( M[i2]^M[i1] );
    }
}

bool bit_is_zero(uint64_t a, int k_bit)  // 0 <= i_bit < 32
//return true if the k-th bits is 0
{
    uint64_t un = 1;
    return ( a & ( un << k_bit ) ) == 0;
}

/********************************************************************/
/********    INVERTIBILITY  --  for square matrices only  ***********/
/********************************************************************/
// This algo check if there is a zero "lead":

bool mat64_invertable(uint64_t *M, int n)    // (i_lead, j_lead) = positions of the lead
{
  int j_lead = 0;
  bool test_stop=false;

  int i=0;

  for (int i_lead = 0; i_lead < n  &&  j_lead < n; i_lead++)   // 
    {
    //cout << "j_lead = " << j_lead << endl;
    i = i_lead;

    //search for the 1rst non-zero element in the column
    while ( bit_is_zero(M[i], j_lead) )  // while M[i][j_lead] = 0
      {
      //cout << "i0 = " << i << endl;
      i++;
      if (i == n)   //all elements are 0  --> no lead in this column, 1.go to the next column and 2.restart
        {
        //cout << "no real j_lead = " << j_lead << "--> not invertable" << endl;
        test_stop=true; break;
        }
      }
    //cout << "i1 = " << i << endl;
    if(test_stop)   {  break;  }
    //cout << "real j_lead =" << j_lead << endl;
    //has found a non-zero element --->  1. increase rank; and 2. swap row i and row i_lead

    if (i != i_lead)  { swap_row_mat(M, i_lead, i, n); }
      //print_mat(M, n, n);  cout << endl;

    i++;
    //put to zero every element under the lead (starting from i = i_lead + 1)
    while(i<n)
      {
      //cout << "i=" << i << endl;
      if ( ! bit_is_zero(M[i], j_lead) )   // if M[i][j_lead] == 1
        { add_row_mat(M, i_lead, i, n);   //L_i <---- L_i XOR L_i_lead
          //print_mat(M, n, n);   cout << endl; 
        }  
      i++;      
      }

    //go to the next column
    j_lead++;
    }

  return (!test_stop);   // if test_stop==true  --> then M is not invertable  (rank < n)
}

//proba_invert_mat64
double sampling_GT64(int n, int Nit=1e2)  //sample Nit GT
{
  clock_t temps = clock();

  uint64_t N = 0;     if (n==64) {  N = ~N ; }  else {  uint64_t un = 1;   N = (un << n) - un;  }
  uint64_t *M = init_mat64(n, n);

  int i_success=0, i_tot=0;
  int compt=0;

  while (i_success<Nit)
  {
    if(mat64_invertable(M, n))
      {
      i_success++; 
      if(i_success%(Nit/10)==0) {cout << "avance : " << compt << "0%" << endl; compt++;}
      //print_mat64(M, n, n); cout << endl;
      }
    randomize_vec64(M, n, N);
    i_tot++;
  }
  cout << "p = " << ((double) i_success)/i_tot << endl; 

  //temps de tourne du programme
    temps = clock()- temps;
    cout << "time sampling GT for uint64:  \t" << (double) temps/CLOCKS_PER_SEC << endl;

  return ((double) i_success)/i_tot;
}

/********************************************************************/
/***************************   matrice REF   ************************/
/********************************************************************/
int mat64_REF(uint64_t *M, int n, int m)    // (i_lead, j_lead) = positions of the lead
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
    while ( bit_is_zero(M[i], j_lead) )  // while M[i][j_lead] == 0
      {
      //cout << "i0 = " << i << endl;
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
    //cout << "i1 = " << i << endl;
    if(test_stop)   {  break;  }
    cout << "real j_lead =" << j_lead << endl;
    //has found a non-zero element --->  1. increase rank; and 2. swap row i and row i_lead
    rank++; 
    if (i != i_lead)  { swap_row_mat(M, i_lead, i, n); }
      print_mat64(M, n, m);  cout << endl;

    i++;
    //put to zero every element under the lead (starting from i = i_lead + 1)
    while(i<n)
      {
      //cout << "i=" << i << endl;
      if ( ! bit_is_zero(M[i], j_lead) )   // if M[i][j_lead] == 1
        { add_row_mat(M, i_lead, i, n);   //L_i <---- L_i XOR L_i_lead
          print_mat64(M, n, m);   cout << endl; 
        }  
      i++;      
      }

    //go to the next column
    j_lead++;
    }

  cout << "Final rank = " << rank << endl; 
  return rank;
}

/********************************************************************/
/**************************   matrice RREF   ************************/
/********************************************************************/
list<unsigned int> mat64_RREF(uint64_t *M, int n, int m)    // (i_lead, j_lead) = positions of the lead
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
    while ( bit_is_zero(M[i], j_lead) )  // while M[i][j_lead] = 0
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
    if (i != i_lead)  { swap_row_mat(M, i_lead, i, n); }
      print_mat64(M, n, m);  cout << endl;

    //put to zero every element under the lead (starting from i = i_lead + 1)
    for (i=0; i < n; i++)
      {
      if ( (! bit_is_zero(M[i], j_lead)) && i!=i_lead) 
          { 
          add_row_mat(M, i_lead, i, n);    //L_i <---- L_i XOR L_i_lead
          print_mat64(M, n, m);   cout << endl; 
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


/********************************************************************/
/**************************   RREF  EXAMPLE   ***********************/
/********************************************************************/
/*
int main()
{
  const unsigned int n=10;
  // ***********************   REF step-by-step   ***********************
      uint64_t *M = init_mat64(n, n);    print_mat64(M, n, n);  cout << endl;

      cout << "Invertable?" << mat64_invertable(M, n) << endl;  //  print_mat(M, n);  cout << endl;

      mat64_REF(M, n, n);   print_mat64(M, n, n);  cout << endl;
      mat64_RREF(M, n, n);  print_mat64(M, n, n);  cout << endl;

      proba_invert_mat64(n, 1e6);    // probability it is invertable and timing
      // !!!! --> probability a bit higher than 2.88  --> because see rand_op()  : operator 0 is directly excluded

  return 0;
}
*/
