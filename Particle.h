#ifndef _PARTICLE_HEADER_
#define _PARTICLE_HEADER_

#include <memory>

#include "Point.h"

class cell; // forward declaration

class particle {
  private:
    point  p_pos, p_dir;              // position and direction of particle
    double p_wgt;                     // particle weight
    bool   exist;                     // true means particle is alive
    std::shared_ptr< cell > p_cell;   // pointer to cell the particle is in
  public:
    particle( point p, point d );     // constructor with position and direction
    ~particle() {};                   // destructor

    point pos() { return p_pos; };    // return particle position
    point dir() { return p_dir; };    // return particle direction 
    double wgt() { return p_wgt; };   // return particle weight
    bool alive() { return exist; };   // return particle state flag
    ray getRay() { return ray( p_pos, p_dir ); }               // return particle position and direction as ray
    std::shared_ptr< cell > cellPointer() { return p_cell; }   // return pointer to cell the particle is in
    void move( double s );            // move particle s units in its current direction
    void scatter( double mu0 );       // change particle direction by cos_t0=mu0 and uniformly sampled azimuth
    void kill();                      // change exist to false
    void setDirection( point p );     // change p_dir and normalize p_dir again
    void adjustWeight( double f );    // multiply weight by f
    void recordCell( std::shared_ptr< cell > cel );            // change p_cell to cel
};

#endif
