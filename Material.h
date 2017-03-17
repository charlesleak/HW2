#ifndef _MATERIAL_HEADER_
#define _MATERIAL_HEADER_

#include <vector>
#include <string>
#include <utility>
#include <memory>

#include "Nuclide.h"

class material {
  private:
    std::string material_name;         // name of material
    double      material_atom_density; // atom density b-1 cm-1
    std::vector< std::pair< std::shared_ptr< nuclide >, double > > nuclides;                           // pairs of nuclide and atom fractions
    double micro_xs();                 // returns micro xs of material for use by macro_xs
  public:
    material( std::string label, double aden ) : material_name(label), material_atom_density(aden) {}; // contructor takes name and atom density
    ~material() {};                    // destructor

    std::string name() { return material_name; }                      // return material name
    double atom_density() { return material_atom_density; }           // return atom density of material
    std::vector< std::pair< std::shared_ptr< nuclide >, double > > getNuclides() { return nuclides; }; // returns the paired list of nuclides
    void   addNuclide( std::shared_ptr< nuclide >, double );          // add a nuclide with its at%
    double macro_xs();                // return the material's macro xs
    std::shared_ptr< nuclide > sample_nuclide();                      // sample nuclide based on cross sections and atom fractions
    std::string sample_collision( particle* p, std::stack<particle>* bank ); // samples nuclide, samples reaction from nuclide, calls reaction's sample method, returns reaction name
};


#endif
