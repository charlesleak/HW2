#ifndef _SOURCE_HEADER_
#define _SOURCE_HEADER_

#include <stack>
#include <memory>

#include "Point.h"
#include "Distribution.h"
#include "Particle.h"

class source {
  private:
    std::shared_ptr< distribution< point > > dist_pos; // distribution of position of source particles
    std::shared_ptr< distribution< point > > dist_dir; // distribution of direction of source particles
  public:
     source( std::shared_ptr< distribution<point> > pos, std::shared_ptr< distribution<point> > dir )
       : dist_pos(pos), dist_dir(dir) {};            // constructor takes dist_pos and dist_dir
    ~source() {};                                    // destructor

    std::stack< particle > sample();                 // returns a bank with one source particle in it
};

#endif
