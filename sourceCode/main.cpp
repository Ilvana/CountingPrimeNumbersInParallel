# include <cstdlib>
# include <iostream>
# include <iomanip>
#include <math.h>
# include <omp.h>

using namespace std;

int main ( int argc, char *argv[] );
void prime_number_sweep ( int n_lo, int n_hi, int n_factor, bool parallel );
int prime_number_par ( int n );
int prime_number_seq ( int n );


int main ( int argc, char *argv[] )
{
  int n_factor;
  int n_hi;
  int n_lo;

  cout << "\n";
  cout << "PRIME_OPENMP\n";
  cout << "  C++/OpenMP version\n";

  cout << "\n";
  cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n";
  cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n";

  n_lo = 1;
  n_hi = 131072;
  n_factor = 2;
  n_lo = 1;
  n_hi = 1000000;
  n_factor = 10;

  prime_number_sweep ( n_lo, n_hi, n_factor, false );
  cout << "\n";
  cout << "PRIME_OPENMP\n";
  cout << "  Normal end of execution.\n";

  return 0;
}

void prime_number_sweep ( int n_lo, int n_hi, int n_factor, bool parallel )
{
  int i;
  int n;
  int primes;
  double wtime, wtime_;

  cout << "\n";
  if ( !parallel )
    cout << "TEST - sequential\n";
  else
    cout << "TEST - parallel\n";
  cout << "  Call PRIME_NUMBER to count the primes from 1 to N.\n";
  cout << "\n";
  cout << "         N        Pi          Time\n";
  cout << "\n";

  n = n_lo;

  while ( n <= n_hi )
  {
    wtime = omp_get_wtime ( );

    if ( !parallel )
        primes = prime_number_seq ( n );
    else
        primes = prime_number_par ( n );

    wtime = omp_get_wtime ( ) - wtime;
    if (wtime < 0)
        wtime = 0;
    float wtime_ = roundf(wtime * 10000) / 10000;


    cout << "  " << setw(8) << n
         << "  " << setw(8) << primes
         << "  " << setw(14) << wtime_ << "\n";

    n = n * n_factor;
  }

  return;
}

//kod za paralelno izvrsenje algoritma
int prime_number_par ( int n )
{
  int i;
  int j;
  int prime;
  int total = 0;
  int niz[4] = {0};

# pragma omp parallel \
  shared ( n, niz ) \
  private ( i, j, prime )

// #pragma omp for reduction (+ : total)
# pragma omp for
  for ( i = 2; i <= n; i++ )
  {
    prime = 1;
    for ( j = 2; j < i; j++ )
    {
      if ( i % j == 0 )
      {
        prime = 0;
        break;
      }
    }
    niz[omp_get_thread_num()] += prime;
  }

  //return total;
  return niz[0] + niz[1] + niz[2] + niz[3];
}

//kod za sekvencijalno izvrsenje algoritma
int prime_number_seq ( int n )
{
  int i;
  int j;
  int prime;
  int total = 0;

  for ( i = 2; i <= n; i++ )
  {
    prime = 1;

    for ( j = 2; j < i; j++ )
    {
      if ( i % j == 0 )
      {
        prime = 0;
        break;
      }
    }
    total = total + prime;
  }

  return total;
}

