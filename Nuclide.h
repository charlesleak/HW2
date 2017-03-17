#ifndef _NUCLIDE_HEADER_
#define _NUCLIDE_HEADER_

#include <vector>
#include <string>
#include <memory>

#include "Reaction.h"

class nuclide {
  private:
    std::string nuclide_name;                        // name of nuclide
    std::vector< std::shared_ptr< reaction > > rxn;  // list of reactions
  public:
    nuclide( std::string label ) : nuclide_name(label) {};                    // constructor takes name
    ~nuclide() {};                                   // destructor

    std::string name() { return nuclide_name; }      // return name of nuclide
    std::vector< std::shared_ptr< reaction > > getReactions() {return rxn;} ; // return list of reactions
    void addReaction( std::shared_ptr< reaction > ); // add a reaction to the list of reactions
    double total_xs();                               // return the total micro xs
    std::shared_ptr< reaction > sample_reaction();   // returns a random reaction based on micro xs
};


#endif
