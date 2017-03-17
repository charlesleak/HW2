// for defining helper functions and holding problem data

#ifndef _SIMULATION_HEADER_
#define _SIMULATION_HEADER_

#include <string>
#include <memory>
#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>

#include "pugixml.hpp"
#include "Distribution.h"
#include "Reaction.h"
#include "Nuclide.h"
#include "Material.h"
#include "Surface.h"
#include "Cell.h"
#include "Source.h"
#include "Particle.h"
#include "Point.h"


// function that returns an item from a vector of objects of type T by name provided
// the object class has a string and a method called name() allowing for it to be returned
template< typename T >
std::shared_ptr< T > findByName( std::vector< std::shared_ptr< T > > vec, std::string name ) {
  for ( auto v : vec ) {
    if ( v->name() == name ) { return v; }
  }
  return nullptr;
}

class simulation {
  private:
    unsigned long long starthist;                                                   // start number of histories to simulate
    unsigned long long endhist;                                                     // total number of histories to simulate
    std::vector< std::shared_ptr< distribution<double> > >  double_distributions;   // all double distributions
    std::vector< std::shared_ptr< distribution<int> > >  int_distributions;         // all int distributions
    std::vector< std::shared_ptr< distribution<point> > >  point_distributions;     // all point distributions
    std::vector< std::shared_ptr< nuclide > > nuclides;                             // all nuclides
    std::vector< std::shared_ptr< material > > materials;                           // all materials
    std::vector< std::shared_ptr< surface > > surfaces;                             // all surfaces
    std::vector< std::shared_ptr< cell > > cells;                                   // all cells

  public:
    std::vector< std::shared_ptr<estimator > > estimators; // BAD PRACTICE TO HAVE PUBLIC DATA I'M SO SORRY
    std::shared_ptr< source > src;                         // the source
    std::string problemName;                               // I MEAN IT I'M VERY SORRY

    simulation( std::string input_file_name );             // constructor takes xml filename and initiates problem
    ~simulation() {};                                      // destructor

    unsigned long long histories() {return endhist; };     // accessing numhist
    void roulette( particle* p, double Ir );               // uses the importance ratio Ir to roulette a particle
    void split( particle* p, double Ir, std::stack< particle >* bank );             // uses the importance ratio to split a particle
    void findResidency( particle* p );                     // find cell the particle is in, changes p_cell
    void changeResidency( particle* p, std::stack< particle >* bank );              // calls findResidency, changes p_wgt, kills particle if necessary
};

#endif
