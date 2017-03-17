#include <cmath>
#include <iostream>

#include "Estimator.h"
#include "Material.h"
#include "Particle.h"
#include "Cell.h"

void surface_current_estimator::score( particle* p ) { tally_hist += p->wgt(); }

void track_length_estimator::score( particle* p ) {
  // instead of scoring a weighted binary value, score a weighted track length divided by volume of cell*
  // dividing by volume of cell is nonsensical for problem 5, so I'll assume you don't actually want flux
  point reverse = point( -1.0 * p->dir().x, -1.0 * p->dir().y, -1.0 * p->dir().z ); // direction opposite particle path
  std::pair< std::shared_ptr< surface >, double > S = p->cellPointer()->surfaceIntersect( ray(p->pos(), reverse ) );
  tally_hist += p->wgt() * S.second; // for flux, would divide by cell volume here
}

void counting_estimator::score( particle* p ) { count_hist++; }

void counting_estimator::endHistory() {
  if ( tally.size() < count_hist + 1 ) { tally.resize( count_hist + 1, 0.0 ); }
  tally[ count_hist ] += 1.0;
  nhist++;
  count_hist = 0;
}

void counting_estimator::report() {
  std::cout << " " << name() << std::endl;
  double s1 = 0.0, s2 = 0.0;
  for ( int i = 0 ; i < tally.size() ; i++ ) {
    double p = tally[i] / nhist;
    std::cout << " " << i << " " << p << "   " << std::sqrt( p * ( 1.0 - p ) / nhist ) / p <<  std::endl;
    s1 += p * i;
    s2 += p * i * i;
  }
  std::cout << "   mean = " << s1 << std::endl;
  std::cout << "   var  = " << s2 - s1*s1 << std::endl;
}

void track_estimator::score( particle* p ) { ntracks ++; }
void track_estimator::endHistory() {  }
void track_estimator::report() { 
  std::cout << " " << name() << "   " << ntracks << std::endl; 
}
