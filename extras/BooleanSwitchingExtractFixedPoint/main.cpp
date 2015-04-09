
// Serialization CPP code
#include <boost/serialization/export.hpp>
#include "database/structures/Grid.h"
#include "database/structures/PointerGrid.h"
#include "database/structures/SuccinctGrid.h"
#include "database/structures/UniformGrid.h"
#include "database/structures/EdgeGrid.h"
#include "database/structures/ParameterSpace.h"
#include "database/structures/EuclideanParameterSpace.h"
#include "database/structures/AbstractParameterSpace.h"
BOOST_CLASS_EXPORT_IMPLEMENT(PointerGrid);
BOOST_CLASS_EXPORT_IMPLEMENT(SuccinctGrid);
BOOST_CLASS_EXPORT_IMPLEMENT(UniformGrid);
BOOST_CLASS_EXPORT_IMPLEMENT(EdgeGrid);
BOOST_CLASS_EXPORT_IMPLEMENT(EuclideanParameter);
BOOST_CLASS_EXPORT_IMPLEMENT(EuclideanParameterSpace);
BOOST_CLASS_EXPORT_IMPLEMENT(AbstractParameterSpace);

#define BS_DEBUG_MODELMAP
#include <iostream>
#include "Model.h"
#include "database/structures/MorseGraph.h"
#include "database/program/jobs/Compute_Morse_Graph.h"
#include "boost/foreach.hpp"
#include "boost/shared_ptr.hpp"
#include "Parameter/BooleanSwitchingParameterSpace.h"
#include "Parameter/FactorGraph.h"
#include "database/structures/ParameterSpace.h"
#include "database/structures/Database.h"

#include <boost/algorithm/string.hpp>


template <typename T>
std::string NumberToString ( T Number );

template <typename T>
std::string NumberToString ( T Number )
{
  std::stringstream ss;
  ss << Number;
  return ss.str();
}

//arguments :

//path/database.mdb path mynetwork.txt

int main ( int argc, char * argv [] ) {
  Database database;
  database . load ( argv [ 1 ] );

  const AbstractParameterSpace & space =
  dynamic_cast<const AbstractParameterSpace&> (database . parameter_space ());

  boost::unordered_set < uint64_t > parameters;
  boost::unordered_set < uint64_t > mgr_indices;


  const std::vector < MGCC_Record > & mgcc_records =
  database . MGCC_Records ();

  // // the index of interest comes from the last two arguments
  // uint64_t mgcc;
  // mgcc = atoi(argv[argc-2]);
  // uint64_t incc;
  // incc = atoi(argv[argc-1]);

  // To store the correspondence betwen the space dimensions and variable names
  std::vector<std::string> variables;
  // To store the parameter index
  std::vector<uint64_t> parameterIndex;
  // To store the parameter inequalities
  std::vector<std::string> parameterInequalities;
  // To store the information of all the walls
  // here switch the <key,value> from Model.h class, to allow search easier
  std::unordered_map <size_t,Wall> walls;
  std::unordered_map <size_t,Wall>::iterator itw;

  // To store the wall information (independent of parameter)
  std::unordered_map<uint64_t,std::string> wallInformation;
  std::unordered_map<uint64_t,std::string>::iterator itWallInfo;

  // To store the outedges between walls for the given incc for each parameter
  std::vector< std::unordered_map< uint64_t, std::vector<uint64_t> > > outedges;

  // Initialize the model
  Model model;
  model . initialize ( argc-1, argv+1 );
  boost::shared_ptr<Grid> phase_space = model . phaseSpace ();
  if ( not phase_space ) {
    throw std::logic_error ( "Clutching_Graph_Job. model.phaseSpace() failed"
                             " to return a valid pointer.\n");
  }
  //
  BooleanSwitchingParameterSpace & boolean_space = *
  boost::dynamic_pointer_cast<BooleanSwitchingParameterSpace> ( model . parameterSpace () );

  walls = model.getWalls();

  // Extract the correspondence network node and variable names
  BooleanSwitching::Network network;
  std::string str1 = argv[2];
  std::string str2 = argv[3];
  std::string str3 = str1 + "/" + str2;
  network.load(str3.c_str());
  // index starts at 1
  for ( unsigned int i=1; i<=network.size(); ++i ) { variables . push_back ( network.name(i) ); }
	std::ofstream ofile;
	ofile . open ( "Variables.txt" );
	for ( unsigned int i=0; i<variables.size(); ++i ) { ofile << i << " " << variables[i] << "\n"; }
  ofile . close();


	ofile . open ( "FixedPoints.txt" );

	for ( unsigned int mgcc=0; mgcc<mgcc_records.size(); ++mgcc ) {

		// store the fixed point wall without repetition
		std::unordered_map < Wall, uint64_t > fpWall;
		std::unordered_map < Wall, uint64_t >::const_iterator itfpWall;

		// Loop through all the mgccp_records
	  for ( unsigned int imgccp=0; imgccp<mgcc_records[mgcc].mgccp_indices.size(); ++imgccp ) {

	    const MGCCP_Record & mgccp_record =
	    database . MGCCP_Records () [ mgcc_records[mgcc].mgccp_indices[imgccp] ];

	    uint64_t morsegraph_index = mgccp_record . morsegraph_index;
	    mgr_indices . insert ( morsegraph_index );
	    const MorseGraphRecord & morsegraph_record = database . morsegraphData () [ morsegraph_index ];

	    uint64_t dag_index = morsegraph_record . dag_index;
	    const DAG_Data & dag_data = database . dagData () [ dag_index ];

	    std::vector < uint64_t > pindex = mgccp_record . parameter_indices;

	    // Loop through all the parameters for a given mgccp
	    for ( unsigned int ip=0; ip<pindex.size(); ++ip ) {
	      //
	      // store the parameter index
	      parameterIndex . push_back ( pindex[ip] );
	      //
	      boost::shared_ptr<BooleanSwitchingParameter> p =
	      boost::dynamic_pointer_cast<BooleanSwitchingParameter> ( boolean_space . parameter ( pindex[ip] ) );
	      //
	      // store the parameter inequalities
	      parameterInequalities . push_back ( boolean_space . prettyPrint(p) );

	      // Get the wall maps (to be used for the outedges )
	      std::vector < std::pair<int64_t,int64_t> > wallMaps = model . getWallMaps ( p );

	//
	//
	// Save the parameter inequalities
	//      parameterfile << boolean_space . prettyPrint ( p );
	//      parameterfile << "\n\n";

	      //
	      boost::shared_ptr<const Map> map = model . map ( p );
	      if ( not map ) {
	        std::cout << "No map with mgccp index " << imgccp << "and with parameter " <<
	        *p << "!.\n";
	      }
	      // Perform Morse Graph computation
	      MorseGraph mg;
	      Compute_Morse_Graph
	       ( & mg,
	         phase_space,
	         map,
	         0,
	         0,
	         0,
	         0 );

				// MAY NEED TO PUT BACK MODEL INIT HERE, MAY HAVE OVERLAP IN THE ANNOTATION
	      model . annotate ( & mg );

	      typedef std::vector<Grid::GridElement> CellContainer;
	      typedef  MorseGraph::VertexIterator VI;
	      VI it, stop;
	      for (boost::tie ( it, stop ) = mg . Vertices (); it != stop;  ++ it ) {
	        // get the incc
	        CS_Data my_cs_data;
	        my_cs_data . vertices . push_back ( *it );
	        uint64_t my_cs_index = database . csIndex ( my_cs_data );
	        // Now we make an INCCP record by hand:
	        INCCP_Record my_inccp_record;
	        my_inccp_record . cs_index = my_cs_index;
	        my_inccp_record . mgccp_index = mgcc_records[mgcc].mgccp_indices[imgccp];
	        uint64_t my_inccp_index = database . inccpIndex ( my_inccp_record );
	        //

					boost::shared_ptr<const Grid> my_subgrid ( mg . grid ( *it ) );
					if ( not my_subgrid ) {
						std::cout << "Abort! This vertex does not have an associated grid!\n";
						abort ();
					}
					CellContainer my_subset = phase_space -> subset ( * my_subgrid );

					BOOST_FOREACH ( Grid::GridElement ge, my_subset ) {
						if ( not boost::dynamic_pointer_cast < AtlasGeo > ( phase_space -> geometry ( ge ) ) ) {
							std::cout << "Unexpected null response from geometry\n";
						}
						AtlasGeo geo = * boost::dynamic_pointer_cast < AtlasGeo > ( phase_space -> geometry ( ge ) );
						RectGeo box =  geo . rect ();
						int id = geo . id ();

						itw = walls . find ( id );
						if ( itw != walls.end() ) {
	          	if ( itw -> second . isFixedPoint() ) {
								// to avoid repetition
								itfpWall = fpWall . find ( itw -> second );
								if ( itfpWall == fpWall.end() ) {
									fpWall [ itw -> second ] = pindex[ip];
									std::string ss;
									ss = "";
									int mydim = itw->second.rect().dimension();
									std::vector<double> lb = itw->second.rect().lower_bounds;
									std::vector<double> ub = itw->second.rect().upper_bounds;
									// for ( unsigned int i=0; i<mydim-1; ++i )  {
									// 	ss += "[" + NumberToString(lb[i]) + ", " + NumberToString(ub[i]) + "]x";
									// }
									// ss += "[" + NumberToString(lb[mydim-1]) + ", " + NumberToString(ub[mydim-1]) + "]";
									//
									for ( unsigned int i=0; i<mydim; ++i )  {
										ss += NumberToString(lb[i]) + " ";
									}

									ofile << mgcc << " " << ss << "\n";
								}
							}
						} else {
							std::cout << "Error could not find wall.\n";
							abort();
						}
					} // end loop mysubset
	      } // end loop over mg.Vertices
	// END DEBUG
	//

	  } // end of loop through parameters

		} // end loop through mgccp

	}


ofile . close ();

  return 0;
}
