#include <iostream>
using namespace std;


//CONSTANTS


//STRUCT
struct Model64
{
    unsigned int n;          // 0 <= n <= 64 
    uint64_t m;     //ranging from 1 to 2^n-1
    uint64_t* phi;  // vector of monomials  
    //set<uint64_t> set_phi;  // phi are ordered from the smallest to the largest --> faster for search for an operator  
};

struct GT64
{
    unsigned int n;          // 0 <= n <= 64
    uint64_t N;		// N= 2^n-1
    uint64_t* sig;  // vector of monomials --> contains GT
    uint64_t* sig_work;  // vector of monomials copy for work, not necessarily GT    
};

struct ModelF2
{
    unsigned int n;          
    uint64_t m;     // ranging from 1 to 2^n-1
    bool** el;  		// bool matrix  
};

struct GTF2
{
    unsigned int n;          
    bool** el;     			// bool squared matrix --> contains GT
    bool** el_work;  		// bool squared matrix --> copy for work, not necessarily GT    
};
