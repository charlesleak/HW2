#ifndef _CELL_HEADER_
#define _CELL_HEADER_

#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <limits>

#include "Point.h"
#include "Surface.h"
#include "Material.h"
#include "Estimator.h"

class cell {
  private:
    std::string cell_name;                                                // name of cell
    std::vector< std::pair< std::shared_ptr< surface >, int > > surfaces; // surfaces defining cell and respective orientation
    std::shared_ptr< material > cell_material;                            // pointer to material in cell
    std::vector< std::shared_ptr< estimator > > cell_estimators;          // estimators tracking in cell
    double importance;                                                    // importance of cell to decide particle weights
  public:

    cell( std::string label ) : cell_name(label) { importance = 1.0; };   // constructor takes name and assumes importance 1.0
    ~cell() {};                                                           // destructor

    std::string name() { return cell_name; };                             // return cell name
    void setMaterial( std::shared_ptr< material > M ) { cell_material = M; };                   // set pointer to material in cell
    std::shared_ptr< material > getMaterial() { return cell_material; }   // return pointer to material in cell
    void setImportance( double imp ) { importance = imp; };               // set importance of cell
    double getImportance() { return importance; }                         // return importance of cell
    void addSurface( std::shared_ptr< surface > S, int sense );           // add a surface defining the cell
    void attachEstimator( std::shared_ptr< estimator > E ) { cell_estimators.push_back( E ); }; // add an estimator
    bool testPoint( point p );                                            // true if point p is inside the cell
    std::pair< std::shared_ptr< surface >, double > surfaceIntersect( ray r );                  // return first surface ray r will intersect and distance to intersection
    double macro_xs() {                                                   // return macro xs of the material in the cell
      if ( cell_material ) { return getMaterial()->macro_xs(); }
      else { return 0.0; }
    };
    void moveParticle( particle* p, double s );                           // move particle to cell edge and scores estimators
    void sampleCollision( particle* p, std::stack<particle>* bank );      // sample collision according to material method
//    double volume();                                                      // return volume of cell
    void scoreEstimators( particle* p );                                  // score cell estimators
};

#endif
