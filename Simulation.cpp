#include "Simulation.h"

// constructor reads in the xml file
simulation::simulation( std::string input_file_name ) {
  // attempt to load
  pugi::xml_document input_file;
  pugi::xml_parse_result load_result = input_file.load_file( input_file_name.c_str() );

  // check to see if result failed and throw an exception if it did
  if ( ! load_result ) {
    std::cout << load_result.description() << std::endl;
    throw;
  }

  // simulation name and number of histories
  pugi::xml_node sim_node = input_file.child("simulation");
  problemName = sim_node.attribute("name").value();
  pugi::xml_node history_node = sim_node.child("histories");
  starthist = history_node.attribute("start").as_ullong();
  endhist = history_node.attribute("end").as_ullong();

  // distributions
  pugi::xml_node input_distributions = input_file.child("distributions");

  // find total number of distributions
  int num_distributions = 0;
  for ( auto d : input_distributions ) { num_distributions++; }

  // since distributions may depend on other distributions, need to iterate
  int set_distributions = 0;
  while ( set_distributions < num_distributions ) {
    int previous_set_distributions = set_distributions;

    for ( auto d : input_distributions ) {
      std::string type = d.name();
      std::string name = d.attribute("name").value();
      std::string data = d.attribute("datatype").value();

      if ( data == "double" ) {
        // skip rest of loop if distribution already done
        if ( findByName( double_distributions, name ) ) { continue; }

        std::shared_ptr< distribution<double> > Dist;
        if ( type == "delta" ) {
          double a = d.attribute("a").as_double();
          Dist = std::make_shared< arbitraryDelta_distribution< double > > ( name, a );
        }
        else if ( type == "uniform" ) {
          double a = d.attribute("a").as_double();
          double b = d.attribute("b").as_double();
          Dist = std::make_shared< uniform_distribution > ( name, a, b );
        }
        else if ( type == "linear" ) {
          double a  = d.attribute("a").as_double();
          double b  = d.attribute("b").as_double();
          double fa = d.attribute("fa").as_double();
          double fb = d.attribute("fb").as_double();
          Dist = std::make_shared< linear_distribution > ( name, a, b, fa, fb );
        }
        else if ( type == "henyeyGreenstein" ) {
          double a = d.attribute("a").as_double();
          Dist = std::make_shared< HenyeyGreenstein_distribution > ( name, a );
        }
        else {
          std::cout << "unsupported distribution with data type " << data << std::endl;
          throw;
        }
        double_distributions.push_back( Dist );
      }
      // integer-valued distributions
      else if ( data == "int" ) {
        // skip rest of loop if distribution already done
        if ( findByName( int_distributions, name ) ) { continue; }

        std::shared_ptr< distribution<int> > Dist;
        if ( type == "delta" ) {
          double a = d.attribute("a").as_int();
          Dist = std::make_shared< arbitraryDelta_distribution< int > > ( name, a );
        }
        else if ( type == "meanMultiplicity" ) {
          double nubar = d.attribute("nubar").as_double();
          Dist = std::make_shared< meanMultiplicity_distribution > ( name, nubar );
        }
        else if ( type == "terrellFission" ) {
          double nubar = d.attribute("nubar").as_double();
          double sigma = d.attribute("sigma").as_double();
          double b     = d.attribute("b").as_double();
          Dist = std::make_shared< TerrellFission_distribution > ( name, nubar, sigma, b );
        }
        else {
          std::cout << "unsupported distribution with data type " << data << std::endl;
          throw;
        }
        int_distributions.push_back( Dist );
      }
      else if ( data == "point" ) {
        // skip rest of loop if distribution already done
        if ( findByName( point_distributions, name ) ) { continue; }

        std::shared_ptr< distribution< point > > Dist;
        if ( type == "delta" ) {
          double x = d.attribute("x").as_double(); 
          double y = d.attribute("y").as_double(); 
          double z = d.attribute("z").as_double();         
          Dist = std::make_shared< arbitraryDelta_distribution< point > > ( name, point( x, y, z ) );
        }
        else if ( type == "isotropic" ) {
          Dist = std::make_shared< isotropicDirection_distribution > ( name );
        }
        else if ( type == "anisotropic" ) {
          double u = d.attribute("u").as_double(); 
          double v = d.attribute("v").as_double(); 
          double w = d.attribute("w").as_double();         
          std::shared_ptr< distribution<double> > angDist = 
            findByName( double_distributions, d.attribute("distribution").value() );
      
          // in the angular distribution has not come yet, skip to the end of the loop
          if ( ! angDist ) { continue; }

          Dist = std::make_shared< anisotropicDirection_distribution > ( name, point( u, v, w ), angDist );
        }
        else if ( type == "independentXYZ" ) {
          std::shared_ptr< distribution<double> > distX = findByName( double_distributions, d.attribute("x").value() ); 
          std::shared_ptr< distribution<double> > distY = findByName( double_distributions, d.attribute("y").value() ); 
          std::shared_ptr< distribution<double> > distZ = findByName( double_distributions, d.attribute("z").value() ); 

          // if any of these distributions have not yet been resolved, skip to the end of the loop
          if ( !distX || !distY || !distZ ) { continue; }

          Dist = std::make_shared< independentXYZ_distribution > ( name, distX, distY, distZ );
        }
	else if ( type == "discrete" ) {
	  std::vector< std::pair< point, double > > pairs;
	  for ( auto dp : d.children() ) {
            point ptemp( dp.attribute("x").as_double(), dp.attribute("y").as_double(), dp.attribute("z").as_double() );
	    std::pair< point, double > pairtemp;
            pairtemp = std::make_pair( ptemp, dp.attribute("p").as_double() );
	    pairs.push_back( pairtemp );
	  }        
          Dist = std::make_shared< arbitraryDiscrete_distribution< point > > ( name, pairs );
	}
        else {
          std::cout << "unsupported " << data << " distribution of type " << type << std::endl;
          throw;
        }
        point_distributions.push_back( Dist );
      }
      else {
        std::cout << "unsupported distribution with data type " << data << std::endl;
        throw;
      }
      // if we reach here, assume distribution has been set
       set_distributions++;
    }
    // check to see if number of distributions has increased, if not, caught in an infinite loop
    if ( previous_set_distributions == set_distributions ) { 
      std::cout << "distributions could not be resolved. " << std::endl;
      throw;
    }
  }

  // iterate over nuclides
  pugi::xml_node input_nuclides = input_file.child("nuclides");
  for ( auto n : input_nuclides ) {
    std::string name = n.attribute("name").value();

    std::shared_ptr< nuclide > Nuc = std::make_shared< nuclide > ( n.attribute("name").value() );
    nuclides.push_back( Nuc );

    // iterate over its reactions
    for ( auto r : n.children() ) {
      std::shared_ptr< reaction > Rxn;
      std::string rxn_type = r.name();

      double xs = r.attribute("xs").as_double();
      if ( rxn_type == "capture" ) {
        Nuc->addReaction( std::make_shared< capture_reaction > ( xs ) );
      }
      else if ( rxn_type == "scatter" ) {
        std::string dist_name = r.attribute("distribution").value();
        std::shared_ptr< distribution<double> > scatterDist = findByName( double_distributions, dist_name );
        if ( scatterDist ) {
          Nuc->addReaction( std::make_shared< scatter_reaction > ( xs, scatterDist ) );
        }
        else {
          std::cout << " unknown scattering distribution " << dist_name << " in nuclide " << name << std::endl;
          throw;
        }
      }
      else if ( rxn_type == "fission" ) {
        std::string mult_dist_name = r.attribute("multiplicity").value();
        std::shared_ptr< distribution<int> > multDist = findByName( int_distributions, mult_dist_name );
        if ( multDist ) {
          Nuc->addReaction( std::make_shared< fission_reaction > ( xs, multDist ) );
        }
        else {
          std::cout << " unknown multiplicity distribution " << mult_dist_name << " in nuclide " << name << std::endl;
          throw;
        }
      }
      else {
        std::cout << "unknown reaction type " << rxn_type << std::endl;
        throw;
      }
    }
  } 

  // iterate over materials
  pugi::xml_node input_materials = input_file.child("materials");
  for ( auto m : input_materials ) {
    std::string name = m.attribute("name").value();
    double      aden = m.attribute("density").as_double();
    
    std::shared_ptr< material > Mat = std::make_shared< material > ( name, aden );    
    materials.push_back( Mat );

    // iterate over nuclides
    for ( auto n : m.children() ) {
      if ( (std::string) n.name() == "nuclide" ) {
        std::string nuclide_name = n.attribute("name").value();
        double      frac         = n.attribute("frac").as_double();
        
        Mat->addNuclide( findByName( nuclides, nuclide_name ), frac );
      }
    }
  }

  // iterate over surfaces
  pugi::xml_node input_surfaces = input_file.child("surfaces");
  for ( auto s : input_surfaces ) {
    std::string type = s.name();

    std::shared_ptr< surface > S;
    if ( type == "plane" ) {
      std::string name = s.attribute("name").value();
      double      a    = s.attribute("a").as_double();
      double      b    = s.attribute("b").as_double();
      double      c    = s.attribute("c").as_double();
      double      d    = s.attribute("d").as_double();
      S = std::make_shared< plane > ( name, a, b, c, d );
    }
    else if ( type == "sphere" ) {
      std::string name = s.attribute("name").value();
      double      x0    = s.attribute("x0").as_double();
      double      y0    = s.attribute("y0").as_double();
      double      z0    = s.attribute("z0").as_double();
      double      rad    = s.attribute("rad").as_double();
      S = std::make_shared< sphere > ( name, x0, y0, z0, rad );
    }
    else if ( type == "cylinderx" ) {
      std::string name = s.attribute("name").value();
      double      y0    = s.attribute("y0").as_double();
      double      z0    = s.attribute("z0").as_double();
      double      rad    = s.attribute("rad").as_double();
      S = std::make_shared< cylinderx > ( name, y0, z0, rad );
    }
    else if ( type == "cylinderz" ) {
      std::string name = s.attribute("name").value();
      double      x0    = s.attribute("x0").as_double();
      double      y0    = s.attribute("y0").as_double();
      double      rad    = s.attribute("rad").as_double();
      S = std::make_shared< cylinderz > ( name, x0, y0, rad );
    }
    else {
      std::cout << " unkown surface type " << type << std::endl;
      throw;
    }

    if ( (std::string) s.attribute("bc").value() == "reflect" ) {
      S->makeReflecting();
    }
    surfaces.push_back( S );
  }

  // iterate over cells
  pugi::xml_node input_cells = input_file.child("cells");
  for ( auto c : input_cells ) {
    std::string name = c.attribute("name").value();

    std::shared_ptr< cell > Cel = std::make_shared< cell > ( name );
    cells.push_back( Cel );

    // cell material
    if ( c.attribute("material") ) {
      std::shared_ptr< material > matPtr = findByName( materials, c.attribute("material").value() );
      if ( matPtr ) {
        Cel->setMaterial( matPtr );
      }
      else {
        std::cout << " unknown material " << c.attribute("material").value() << " in cell " << name << std::endl;
        throw;
      } 
   }

    // cell importance
    if ( c.attribute("importance") ) {
      Cel->setImportance( c.attribute("importance").as_double() );
    }
   
    // iterate over surfaces
    for ( auto s : c.children() ) {
      if ( (std::string) s.name() == "surface" ) {
        std::string name  = s.attribute("name").value();
        int         sense = s.attribute("sense").as_int();

        std::shared_ptr< surface > SurfPtr = findByName( surfaces, name );
        if ( SurfPtr ) {
          Cel->addSurface( findByName( surfaces, name ), sense );
        }
        else {
          std::cout << " unknown surface with name " << name << std::endl;
          throw;
        }
      }
      else {
        std::cout << " unknown data type " << s.name() << " in cell " << name << std::endl;
        throw;
      }
    } 
  }

  // iterate over estimatators
  pugi::xml_node input_estimators = input_file.child("estimators");
  for ( auto e : input_estimators ) {
    std::string type = e.name();
    std::string name = e.attribute("name").value();
    std::shared_ptr< estimator > Est;
    if ( type == "current" ) {
      Est = std::make_shared< surface_current_estimator > ( name );
      // get the surfaces
      for ( auto s : e.children() ) {
        if ( (std::string) s.name() == "surface" ) {
          std::string name = s.attribute("name").value();
          std::shared_ptr< surface > SurfPtr = findByName( surfaces, name );
          if ( SurfPtr ) {
            SurfPtr->attachEstimator( Est );
          }
          else {
            std::cout << " unknown surface label " << name << " in estimator " << e.attribute("name").value() << std::endl;
          }
        }
      } 
    }
    else if ( type == "countingSurface" ) {
      Est = std::make_shared< counting_estimator > ( name );
      // get the surfaces
      for ( auto s : e.children() ) {
        if ( (std::string) s.name() == "surface" ) {
          std::string name = s.attribute("name").value();
          std::shared_ptr< surface > SurfPtr = findByName( surfaces, name );
          if ( SurfPtr ) {
            SurfPtr->attachEstimator( Est );
          }
          else {
            std::cout << " unknown surface label " << name << " in estimator " << e.attribute("name").value() << std::endl;
          }
        }
      } 
    }
    else if ( type == "trackLength" ) {
      Est = std::make_shared< track_length_estimator > ( name );
      // get the cells
      for ( auto s : e.children() ) {
        if ( (std::string) s.name() == "cell" ) {
          std::string name = s.attribute("name").value();
          std::shared_ptr< cell > CellPtr = findByName( cells, name );
          if ( CellPtr ) {
            CellPtr->attachEstimator( Est );
          }
          else {
            std::cout << " unknown cell label " << name << " in estimator " << e.attribute("name").value() << std::endl;
          }
        }
      }
    } else if ( type == "track" ) {
      Est = std::make_shared< track_estimator > ( name );
      // get the cells
      for ( auto c : cells ) {
        c->attachEstimator( Est );
      }
    }
    else {
      std::cout << "unknown estimator type " << name << std::endl;
      throw;
    }
    estimators.push_back( Est );
  }

  // create source
  pugi::xml_node input_source = input_file.child("source");
  pugi::xml_node input_source_position  = input_source.child("position");
  pugi::xml_node input_source_direction = input_source.child("direction");

  std::string pos_dist_name = input_source_position.attribute("distribution").value();
  std::string dir_dist_name = input_source_direction.attribute("distribution").value();

  std::shared_ptr< distribution< point > > posDist = findByName( point_distributions, pos_dist_name );
  std::shared_ptr< distribution< point > > dirDist = findByName( point_distributions, dir_dist_name );

  if ( posDist && dirDist ) {
    src = std::make_shared< source > ( posDist, dirDist );  
  }
  else {
    if ( ! posDist ) { std::cout << " unknown position distribution "  << pos_dist_name << " in source " << std::endl; }
    if ( ! dirDist ) { std::cout << " unknown direction distribution " << dir_dist_name << " in source " << std::endl; }
    throw;
  }

}

// rouletting a particle
void simulation::roulette( particle* p, double Ir ) {
  if ( Urand() < Ir ) { p->kill(); }
  else { p->adjustWeight( 1.0 / Ir ); }
}

// splitting a particle
void simulation::split( particle* p, double Ir, std::stack< particle >* bank ) {
  double N = std::floor( Ir + Urand() ); // split particle into N particles
  for (int i = 0; i < N-1; i++ ) {           // make N-1 new particles
    particle pTemp( p->pos(), p->dir() );
    pTemp.adjustWeight( p->wgt() / N );  // with reduced weight
    bank->push( pTemp );
  }
  p->adjustWeight( 1.0 / N );            // reduce weight of current (Nth) particle
}

// find the new residency of particle and sets p_cell
void simulation::findResidency( particle* p ) {
  for ( auto c : cells ) {
    if ( c->testPoint( p->pos() ) ) {
      p->recordCell(c);
    }
  }
}

// change residency of particle function
void simulation::changeResidency( particle* p, std::stack< particle >* bank ) {
  double I1 = p->cellPointer()->getImportance(); // importance of resident cell before move
  findResidency( p );                            // changes the p_cell
  double Ir = p->cellPointer()->getImportance() / I1; // ratio of importances of resident cells before and after move
  if ( Ir == 0 ) { p->kill(); }                       // particle entered a void and needed to be killed
  else if ( Ir < 1.0 ) { roulette( p, Ir ); }
  else if ( Ir > 1.0 ) { split( p, Ir, bank ); }
}
