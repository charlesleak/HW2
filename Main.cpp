// work in progress main for HW2
#include <memory>
#include <iostream>
#include <stack>
#include <vector>
#include <limits>
#include <string>
#include <cmath>
#include <ctime>

#include "Random.h"
#include "Distribution.h"
#include "Point.h"
#include "Particle.h"
#include "Reaction.h"
#include "Nuclide.h"
#include "Material.h"
#include "Surface.h"
#include "Cell.h"
#include "Simulation.h"

int main() {

  // user enters the XML file name
  std::string input_file_name;
  std::cout << " Enter XML input file name: " << std::endl;
  std::cin >> input_file_name;

  // load and initialize problem
  simulation sim( input_file_name );

  // for timing
  std::clock_t start = std::clock();

  // simulation loop through all histories
  double sci1 = sim.histories() / std::pow( 10, std::floor( std::log10( sim.histories() ) ) ); // to print scientific notation
  double sci2 = std::floor( std::log10( sim.histories() ) );                                   // to print scientific notation
  std::cout << " Running " << sim.problemName << " for " << sci1 << "E" << sci2 << " histories." << std::endl;
  for ( unsigned long long history = 0 ; history < sim.histories() ; history++ ) {

    // create a new particle from source distributions, make bank, and deposit it in bank
    std::stack< particle > bank = sim.src->sample();

    // loop for a single history
    while ( ! bank.empty() ) {

      // take a particle from the bank
      particle p = bank.top();
      sim.findResidency( &p ); //determine and assign p_cell
      bank.pop();

      while ( p.alive() ) { // particle loop

        // determine its next action, either media interaction or boundary crossing
        double dist_collision = -std::log( Urand() ) / p.cellPointer()->macro_xs();
        std::pair< std::shared_ptr< surface >, double > S = p.cellPointer()->surfaceIntersect( p.getRay() );
        double dist_surface = S.second;
        double distance = std::fmin( dist_collision, dist_surface );

        // move particle, calling cell estimators
        p.cellPointer()->moveParticle( &p, distance );

        // check if particle left cell
        if ( distance == dist_surface ) {
          // cross surface, calling estimator
          S.first->crossSurface( &p );
          // find which cell particle's in, change p_cell, roulette or split, or kill if void
          sim.changeResidency( &p, &bank );
        }

        // if it didn't leave cell, it had a collision in the cell
        else {
          // sample nuclide and reaction
          p.cellPointer()->sampleCollision( &p, &bank );
        }
        
      } // end particle loop

    } // end history loop

    // print timer
    if ( ( fmod( std::log10( history + 1 ), 1 ) == 0 ) || ( history + 1 == sim.histories() ) ) {
      double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
      if ( duration != 0.0 ) {
	// to print scientific notation
        sci1 = ( history + 1 ) / std::pow( 10, std::floor( std::log10( history + 1 ) ) );
        sci2 = std::floor( std::log10( history + 1 ) );
        std::cout << "  " << sci1 << "E" << sci2 << " histories took " << duration << " seconds to run.";
	if ( history + 1 == sim.histories() ) {std::cout << " Simulation finished around "; }
	else {std::cout << " Simulation should finish around "; }
        // predict time left
        double timeLeft = duration / (history + 1) * ( sim.histories() - ( history + 1 ) );
        // print local time + time left
        time_t rawtime;
        struct tm * timeinfo;
        time(&rawtime);
        rawtime += timeLeft;
        timeinfo = localtime (&rawtime);
        std::cout << asctime(timeinfo);
      }
    }

    // tally closeout: a history has been completed
    for ( auto e : sim.estimators ) { e->endHistory(); }

  } // end simulation loop

  std::cout << " Done." << std::endl;
  for ( auto e : sim.estimators ) { e->report(); }

  return 0;
}
