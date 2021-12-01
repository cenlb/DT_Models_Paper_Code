#ifndef he_pq_H
#define he_pq_H

#include <algorithm>
#include <random>
#include <iostream>
#include <queue>
#include <tuple>

typedef std::tuple< double, char, int > host_event; //time, type, host-id
//typedef std::tuple< double, long int, double > pret; //time, P_t, cumhaz
typedef std::tuple< bool, double, long int > pret; //exposure?, t, P_t

/*
std::priority_queue does not have an erase member fn for arb. element

solution:

alexm:
https://stackoverflow.com/questions/19467485/how-to-remove-element-not-at-top-from-priority-queue

21 July 2020
*/

//auto mycomp = [](const host_event& he1, const host_event& he2){return (std::get<0>(he1) < std::get<0>(he2));};

class hecomp
{
 public:
  bool operator()(const host_event& he1, const host_event& he2)
  { return (std::get<0>(he1) > std::get<0>(he2)) or
      ((std::get<0>(he1) == std::get<0>(he2)) and (std::get<1>(he1) > std::get<1>(he2)) );
      }
};



class he_priority_queue : public std::priority_queue< host_event, std::vector<host_event>, hecomp >
{
 public:
  //he_priority_queue(hecomp& myco)
  // {std::priority_queue< host_event, std::vector<host_event>, hecomp >::comp = myco;};
    
  
  
  host_event remove(int host_id, char ev_type)
  {
    host_event ret;
    std::vector< host_event >::iterator it
      = std::find_if(c.begin(), c.end()
		     ,[host_id, ev_type](host_event he)
		     { return ( std::get<1>(he) == ev_type and std::get<2>(he) == host_id ); });

    if (it != this->c.end())
      {
	ret = *it;
	this->c.erase(it);
	std::make_heap(this->c.begin(), this->c.end(), this->comp);
	return ret;
      }

    else
      {
	return host_event{};
      }
  }

    void print( std::ostream& os )
  {
    std::vector< host_event > tmp;
    for( std::vector< host_event >::const_iterator it{c.cbegin()}, itend{c.cend()};
	 it != itend;
	 ++it
	 )
      tmp.push_back(*it);

    sort(tmp.begin(), tmp.end());

    for( host_event he : tmp ) os << std::get<0>(he) << "\t" << std::get<1>(he) << "\t"
				  << std::get<2>(he) << "\n";
    os << "\n";
  }
};

long int one_step(
		  long int P
		  ,double eps
		  ,double rho
		  ,int I
		  ,double delta_t
		  ,std::mt19937_64& rng
		  );


pret draw_exposure(
		   pret& pr
		   ,he_priority_queue& hepq
		   ,double alp
		   ,double bet
		   ,double eps
		   ,double rho
		   ,int S
		   ,int I
		   ,double delta_t
		   ,std::mt19937_64& rng
		   );

/*! Removes host_id from vector of EXPOSED indices and adds it to vector of INFECTED indices
  Sorts indI
*/
void adjust_state_I(
		    int host_id
		    ,std::vector< int >& indE
		    ,std::vector< int >& indI
		    );

/*! Removes host_id from vector of INFECTED indices
*/
void adjust_state_R(
		    int host_id
		    ,std::vector< int >& indI
		    );
/*! Removed nrem_E uniformly & randomly from indE and nrem_I similarly from indI
  and all related host events from the queue
 */
void adjust_state_X(
		    std::size_t nrem_E
		    ,std::size_t nrem_I
		    ,std::vector< int >& indE
		    ,std::vector< int >& indI
		    ,he_priority_queue& hepq
		    ,std::mt19937_64& rng
		    );
#endif //he_pq_H

/*
char g(
       std::tuple< double, long int, double >& Q // t, P_t, cumhaz
       ,double alp
       ,double bet
       ,double eps
       ,double rho
       ,int S
       ,int I
       ,double delta_t
       ,double theta
       ,double T
       ,std::mt19937_64& rng
       );
*/
