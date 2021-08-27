#include <iostream>
#include <random>
#include <cmath>
#include <tuple>
#include <queue>
#include <utility>
#include <vector>
#include <cassert>
#include <fstream>
#include <string>

using namespace std;


double V;
double alp, bet, del, gam, eps, rho;
double theta, cumhaz;

int N, S, I;
long int P;

double t, delta_t, t_max, T;
double e, i, r;

vector< double> E_times, I_times, R_times;

random_device rD;
random_device::result_type seed;



//bool host_event_compare( tuple< double, char > he1, tuple< double, char > he2 )
//{
// return (get<0>(he1) > get<0>(he2));
//}

auto host_event_compare = []( const tuple< double, char >& he1, const tuple< double, char >& he2 )
{ return get<0>(he1) > get<0>(he2); };

typedef priority_queue< tuple< double, char >, vector< tuple< double, char >>, decltype(host_event_compare)> host_event_pq;
host_event_pq hepq(host_event_compare);

 
long int one_step(
		  long int P
		  ,double eps
		  ,double rho
		  ,int I
		  ,double delta_t
		  ,mt19937_64& rng
		  )
{
  static binomial_distribution< long int > bin;
  static poisson_distribution< long int > poi;

  long int P1, P2;

  if( P == 0 ) P1 = 0;
  else
    {
        bin.param(binomial_distribution<long int>::param_type{ P, exp( -rho * delta_t )});
	P1 = bin(rng);
    }
  if( I == 0 ) P2 = 0;
  else
    {
        poi.param(poisson_distribution<long int>::param_type{ eps * I / rho * ( 1 - exp( -rho * delta_t ))});
	P2 = poi(rng);
    }

  //cout << P1 + P2 << "\n";
  return P1 + P2;
}

char g(
       tuple< double, long int, double >& Q // t, P_t, cumhaz
       ,double alp
       ,double bet
       ,double eps
       ,double rho
       ,int S
       ,int I
       ,double delta_t
       ,double T
       ,mt19937_64& rng
       )
{
  double t = get<0>(Q);
  long int P = get<1>(Q);
  double cumhaz = get<2>(Q);
  double delta_cumhaz;
  
  if( I > 0 )
    {
      //stop and return
      // (a) when t hits T ...= min(hepq.top().time(), t_max)
      // (b) when cumhaz hits theta
      while(true)
	{
	  if(t + delta_t >= T)
	    {
	      get<0>(Q) = T;
	      get<1>(Q) = P;
	      get<2>(Q) = cumhaz;
	      return 'a';
	    }
	  
	  delta_cumhaz = (alp*P + bet*I) * S * delta_t;
	  if( cumhaz + delta_cumhaz >= theta )
	    {
	      get<0>(Q) = t + 0.5*delta_t;
	      get<1>(Q) = P;
	      get<2>(Q) = cumhaz + 0.5*delta_cumhaz;
	      return 'c';
	    }

	  cumhaz += delta_cumhaz;
	  P = one_step( P, eps, rho, I, delta_t, rng );
	  t += delta_t;
	  //cout << t << "\t" << P << "\t" << cumhaz << "\t" << theta << "\n";
        }
    }
  else if( I == 0 )
    {
      //stop and return
      // (a) when t hits T... = min(hepq.top().time(), t_max)
      // (c) when cumhaz hits theta
      // (z) when P_t = 0
      while(true)
	{
	  if(t + delta_t >= T)
	    {
	      get<0>(Q) = T;
	      get<1>(Q) = P;
	      get<2>(Q) = cumhaz;
	      return 'a';
	    }
	  
	  delta_cumhaz = (alp*P + bet*I) * S * delta_t;
	  if( cumhaz + delta_cumhaz >= theta )
	    {
	      get<0>(Q) = t + 0.5*delta_t;
	      get<1>(Q) = P;
	      get<2>(Q) = cumhaz + 0.5*delta_cumhaz;
	      return 'c';
	    }

	  if( P == 0 )
	    {
	      get<0>(Q) = t;
	      get<1>(Q) = P;
	      get<2>(Q) = cumhaz;
	      return 'z';
	    }
	  
	  cumhaz += delta_cumhaz;
	  P = one_step( P, eps, rho, I, delta_t, rng );
	  t += delta_t;
	  //cout << t << "\t" << P << "\t" << cumhaz << "\t" << theta << "\n";
	}
    }
  else
    {      return 'x'; //test for fall through but shouldn't get to this
    }
}

int adjust_state( int& I, char ev )
{
  if( ev == 'I' ) {I += 1; return 0;}
  else if( ev == 'R' ) {I -= 1; return 0;}
  else return 1;
}
  

char gret;

int main( int argc, char* argv[] )
{

  if( argc != 12 )
    {
      cout << "Usage:wssv \npath \nfile_no\n alp\t bet\tdel\tgam\teps\trho\ndelta_t\tT\tN\n";
      return 1;
    }

  string path{argv[1]};
  string f_no{argv[2]};

  alp = stod(argv[3]);
  bet = stod(argv[4]);
  del = stod(argv[5]);
  gam = stod(argv[6]);
  eps = stod(argv[7]);
  rho = stod(argv[8]);

  delta_t = stod(argv[9]);
  T = stod(argv[10]);
  N = stoi(argv[11]);
    
  
  ofstream f;
  
  seed = rD();
  //seed = 2300143884;
  cout << "RNG seed: " << seed << "\n";
  
  mt19937_64 ranD{seed};
  
  //alp = 3.33e-11;
  //bet = 3.33e-10;
  //del = 0.005;
  //gam = 0.05;

  //N = 80;
  S = N - 1;
  I = 0;

  //eps = 2e5;
  //rho = 0.005;

  //t_max = 5000.0;
  //delta_t = 0.01;

  exponential_distribution< double > E_sampler{del};
  exponential_distribution< double > I_sampler{gam};
  exponential_distribution< double > theta_sampler{1.0};

  tuple< double, long int, double > Q{0.0, 0L, 0.0};
  bool save_and_exit{false};
  
  e = get<0>(Q);
  i = e + E_sampler(ranD);
  r = i + I_sampler(ranD);
  
  //hepq.push( tuple< double, char >{ e, 'E' });
  hepq.push( tuple< double, char >{ i, 'I' });
  hepq.push( tuple< double, char >{ r, 'R' });

  E_times.push_back(e);
  I_times.push_back(i);
  R_times.push_back(r);

  /*
  cout << e << "\t" << i << "\t" << r << "\n";
  cout << "===========================\n"
       << "S: " << S << "\n"
       << "I: " << I << "\n"
       << "t: " << get<0>(Q) << "\n"
       << "P: " << get<1>(Q) << "\n"
       << "H: " << get<2>(Q) << "\n";*/

  theta = theta_sampler(ranD);

  while( S > 0 )
    {
      if(save_and_exit) break;
  
      T = hepq.empty() ? t_max : get<0>(hepq.top());
      gret = g( Q, alp, bet, eps, rho, S, I, delta_t, T, ranD );

      
      if( hepq.empty() and gret == 'z' )
	{
	  //cout << "z : pq empty\n";
	  save_and_exit = true;
	}
      
      else if( !hepq.empty() and gret == 'z' )
	{
	  //cout << "z : pq not empty\n";
	  get<0>(Q) = get<0>( hepq.top() );
	  adjust_state( I, get<1>( hepq.top() ));
	  hepq.pop();		
	}
      else if( hepq.empty() and gret == 'a' )
	{
	  //cout << "a : pq empty\n";
	  save_and_exit = true;
	}
      else if( !hepq.empty() and gret == 'a' )
	{
	  //cout << "a : pq not empty\n";
	  get<0>(Q) = get<0>( hepq.top() );
	  adjust_state( I, get<1>( hepq.top() ));
	  hepq.pop();		
	}
      else if( gret == 'c' )
	{
	  //cout << "c cumhaz\n";
	  S -= 1;
	  e = get<0>(Q);
	  i = e + E_sampler(ranD);
	  r = i + I_sampler(ranD);
	  
	  //cout << "(e,i,r) = " << e << "," << i << "," << r << "\n";
	  hepq.push( tuple< double, char >{ i, 'I' });
	  hepq.push( tuple< double, char >{ r, 'R' });
	  
	  E_times.push_back(e);
	  I_times.push_back(i);
	  R_times.push_back(r);

	  theta = theta_sampler(ranD);
	  get<2>(Q) = 0.0;
        }
      else
	{
	  //cout << "Error - fall through\n";
	}

      /*
      cout << "===========================\n"
	   << "S: " << S << "\n"
	   << "I: " << I << "\n"
	   << "t: " << get<0>(Q) << "\n"
	   << "P: " << get<1>(Q) << "\n"
	   << "H: " << get<2>(Q) << "\n";*/
      
    }
  
  assert( E_times.size() == I_times.size() and I_times.size() == R_times.size() );

  f.open(path + "/out_E" + "_" + f_no);
  for( double x : E_times ) f << x << "\n";
  f.close();

  f.open(path + "/out_I" + "_" + f_no);
  for( double x : I_times ) f << x << "\n";
  f.close();

  f.open(path + "/out_R" + "_" + f_no);;
  for( double x : R_times ) f << x << "\n";
  f.close();

  
      
}


/*
	    S -= 1;
	    e = get<0>(Q);
	    i = e + E_sampler(ranD);
	    r = i + I_sampler(ranD);

	    cout << "(e,i,r) = " << e << "," << i << "," << r << "\n";
	    hepq.push( tuple< double, char >{ i, 'I' });
	    hepq.push( tuple< double, char >{ r, 'R' });
	    
	    E_times.push_back(e);
	    I_times.push_back(i);
	    R_times.push_back(r);


	    cout << "S: " << S << "\n"
	    << "I: " << I << "\n"
	    << "t: " << get<0>(Q) << "\n"
	    << "P: " << get<1>(Q) << "\n"
	    << "H: " << get<2>(Q) << "\n";
	    
*/
