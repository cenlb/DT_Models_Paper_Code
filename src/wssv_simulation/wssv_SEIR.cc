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

double bet, del, gam;
double theta;

int N, S, I;
long int P;

double t, T;
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

 
int adjust_state( int& I, char ev )
{
  if( ev == 'I' ) {I += 1; return 0;}
  else if( ev == 'R' ) {I -= 1; return 0;}
  else return 1;
}
  
int main( int argc, char* argv[] )
{

  if( argc != 7 )
    {
      cout << "Usage:wssv path file_no bet del gam N\n";
      return 1;
    }

  string path{argv[1]};
  string f_no{argv[2]};

  bet = stod(argv[3]);
  del = stod(argv[4]);
  gam = stod(argv[5]);

  N = stoi(argv[6]);
    
  
  ofstream f;
  
  seed = rD();
  //seed = 2300143884;
  cout << "RNG seed: " << seed << "\n";
  
  mt19937_64 ranD{seed};
  
  S = N - 1;
  I = 0;

  exponential_distribution< double > E_sampler{del};
  exponential_distribution< double > I_sampler{gam};
  exponential_distribution< double > theta_sampler;

  bool save_and_exit{false};

  t = 0.0;
  e = t;
  i = e + E_sampler(ranD);
  r = i + I_sampler(ranD);
  
  hepq.push( tuple< double, char >{ i, 'I' });
  hepq.push( tuple< double, char >{ r, 'R' });

  E_times.push_back(e);
  I_times.push_back(i);
  R_times.push_back(r);

  
  while( S > 0 )
    {
      if(save_and_exit) break;

      if( I == 0 and hepq.empty() )
	{
	  //cout << "I == 0 and hepq.empty()\n";
	  save_and_exit = true;
	}
      else if( I == 0 and hepq.empty() )
	{
	  //cout << "I == 0 and !hepq.empty()\n";
	  t = get<0>(hepq.top());
	  adjust_state( I, get<1>(hepq.top()) );
	}
      else //I > 0
	{
	  //cout << "I > 0\n";
	  theta_sampler.param(exponential_distribution<double>::param_type{ bet * S * I });
	  T = t + theta_sampler(ranD);
	  if( !hepq.empty() and T >= get<0>(hepq.top()) )
	    {
	      t = get<0>(hepq.top());
	      adjust_state( I, get<1>(hepq.top()) );
	      hepq.pop();
	    }
	  else
	    {
	      S -= 1;
	      t = T;
	      e = t;

	      i = e + E_sampler(ranD);
	      r = i + I_sampler(ranD);
  
	      hepq.push( tuple< double, char >{ i, 'I' });
	      hepq.push( tuple< double, char >{ r, 'R' });

	      E_times.push_back(e);
	      I_times.push_back(i);
	      R_times.push_back(r);

	      /*cout << e << "\t" << i << "\t" << r << "\n";*/

	    }
	  
	}

      
      /*cout << "===========================\n"
	   << "S: " << S << "\n"
	   << "I: " << I << "\n"
	   << "t: " << t << "\n";*/
      
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
