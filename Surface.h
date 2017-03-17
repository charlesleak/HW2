#ifndef _SURFACE_HEADER_
#define _SURFACE_HEADER_

#include <string>
#include <vector>
#include <limits>

#include "Point.h"
#include "Particle.h"
#include "Estimator.h"

class surface {
  private:
    std::string surface_name;  // name of surface
    bool reflect_bc;           // true if reflecting boundary
    std::vector< std::shared_ptr< estimator > > surface_estimators;             // estimators
  public:
    surface( std::string label ) : surface_name(label) { reflect_bc = false; }; // constructor takes name
    ~surface() {};

    virtual std::string name()    final { return surface_name; };               // return name
    virtual void makeReflecting() final { reflect_bc = true; };                 // make reflector

    virtual void attachEstimator( std::shared_ptr< estimator > E ) final {      // add estimator
      surface_estimators.push_back( E );
    }
    virtual void scoreEstimators( particle* p ) final {                         // uses estimator's score method
      for ( auto e : surface_estimators ) { e->score( p ); }
    }

    virtual void crossSurface( particle* p ) final {                            // scores estimators, reflects, nudges particle
      // score estimators
      for ( auto e : surface_estimators ) { e->score( p ); }

      // reflect if needed
      if ( reflect_bc ) { 
        point d = reflect( p->getRay() );
        p->setDirection( d );
      }

      // advance particle off the surface
      p->move( std::numeric_limits<float>::epsilon() );
    }

    virtual double eval( point p )   = 0;       // pure virtual
    virtual double distance( ray r ) = 0;       // pure virtual
    virtual point  reflect( ray r )  = 0;       // pure virtual
};

class plane : public surface {
  private:
    double a, b, c, d;
  public:
    plane( std::string label, double p1, double p2, double p3, double p4 ) :   // constructor takes name and plane equation
      surface(label), a(p1), b(p2), c(p3), d(p4) {};
    ~plane() {};               // destructor

    double eval( point p );    // return positive, zero or negative
    double distance( ray r );  // return min positive distance to intersection
    point  reflect( ray r );   // return new reflected direction
};

class sphere : public surface {
  private:
    double x0, y0, z0, rad;
  public:
    sphere( std::string label, double p1, double p2, double p3, double p4 ) :   // constructor takes name and sphere equation
      surface(label), x0(p1), y0(p2), z0(p3), rad(p4) {};
    ~sphere() {};             // destructor

    double eval( point p );   // return positive, zero or negative
    double distance( ray r ); // return min positive distance to intersection
    point  reflect( ray r );  // return new reflected direction
};

class cylinderx : public surface { // cylinder parrallel to x axis
  private:
    double y0, z0, rad;
  public:
    cylinderx( std::string label, double p1, double p2, double p3 ) :   // constructor takes name and cylinder equation
      surface(label), y0(p1), z0(p2), rad(p3) {};
    ~cylinderx() {};             // destructor

    double eval( point p );   // return positive, zero or negative
    double distance( ray r ); // return min positive distance to intersection
    point  reflect( ray r );  // return new reflected direction
};

class cylinderz : public surface { // cylinder parrallel to z axis
  private:
    double x0, y0, rad;
  public:
    cylinderz( std::string label, double p1, double p2, double p3 ) :   // constructor takes name and cylinder equation
      surface(label), x0(p1), y0(p2), rad(p3) {};
    ~cylinderz() {};             // destructor

    double eval( point p );   // return positive, zero or negative
    double distance( ray r ); // return min positive distance to intersection
    point  reflect( ray r );  // return new reflected direction
};

#endif
