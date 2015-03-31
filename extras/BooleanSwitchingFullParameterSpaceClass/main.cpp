
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

#include "GI.h"

//arguments :

//path/database.mdb path2 mynetwork.txt

std::string equivalencyString ( const std::string & mystr ) {

  std::unordered_map < std::string, std::string > annotationsEquivalency ( {
              {"FP ","FP "},
              {"FP ON ", "FP "},
              {"FP OFF ", "FP "},
              {"FC ", "FC "}
  } );

  std::string result  = annotationsEquivalency[mystr];

  return result;

}

// ==============================================================
// When extracting the Symbol for the Morse Graph nodes, we apply
// a rule for equivalence between the existing annotations.
// ==============================================================
//
// split a string into different fields delimited by a colon ":"
// return the first field found
std::string extractSymbol ( const std::string & string ) {
  std::vector<std::string> fields;
  boost::split ( fields, string, boost::is_any_of (":") );
  if ( fields . size () >= 1 ) return "\""+equivalencyString(fields[0])+"\"";
  else return "\"\""; // DEBUG
}
std::string extractSymbolWithoutQuotes ( const std::string & string ) {
  std::vector<std::string> fields;
  boost::split ( fields, string, boost::is_any_of (":") );
  if ( fields . size () >= 1 ) return equivalencyString(fields[0]);
  else return std::string (); // DEBUG
}

std::string extractConditionalString ( const std::string & string ) {
  std::vector<std::string> fields;
  boost::split ( fields, string, boost::is_any_of (":") );
  // convention for the condition string from AnnotationConditions.h
  // symbol : description : extra information
  // the conditional string is simply symbol : description
  if ( fields . size () >= 2 ) return fields[0] + ":" + fields[1];
  else return std::string (); // DEBUG
}

std::string constructLabel ( const std::string & string ) {
  std::vector<std::string> fields;
  boost::split ( fields, string, boost::is_any_of (":") );
  if ( fields.size() == 2 ) {
    // return "&#92;n " + extractSymbolWithoutQuotes(string);
    return extractSymbolWithoutQuotes(string);
  }
  if ( fields.size() == 3 ) {
    std::string str;
    // This case should be checked
    str = "&#92;n " + extractSymbolWithoutQuotes(string) + "&#92;n { ";
    std::vector<std::string> extrafields;
    boost::split ( extrafields, fields[2], boost::is_any_of (" ") );
    // Warning : need -2 because there is an extra space at the end of the string
    for ( unsigned int i=0; i<extrafields.size()-2; ++i ) {
      str += extrafields[i] + ", ";
    }
    str += extrafields[extrafields.size()-2] + " }";
    return str;
  } else {
    std::cout << "string = " << string << "\n";
    std::cout << "condition string format not implemented\n";
    return std::string (); // TODO: Arnaud, could you take a look at this? 9-3-2014
    //abort();
  }
}

//
// annotation convention   symbol : text
// example : "FP : Morse set is a fixed point"
// we scan the string to find the symbol to label the node in graphviz
//
std::string makeLabel ( const std::vector < std::string > & myannotations ) {
  std::string output;
  output = "";
  if ( myannotations.size() == 0 ) {
    std::cout << "No annotations.\n";
    abort();
  }
  for ( unsigned int i=0; i<myannotations.size(); ++i ) {
    // we concatenate all the labels from all the annotations
    output += constructLabel(myannotations[i]);
  }
  return output;
}




DAG makeDAG ( const Database & database, uint64_t mgcc ) {
	DAG result;
	MGCC_Record const& mgcc_record = database . MGCC_Records () [ mgcc ];
  uint64_t mgccp_index = mgcc_record . mgccp_indices [ 0 ];
  const MGCCP_Record & mgccp_record = database . MGCCP_Records () [ mgccp_index ];
  uint64_t morsegraph_index = mgccp_record . morsegraph_index;
  const MorseGraphRecord & morsegraph_record =
  database . morsegraphData () [ morsegraph_index ];
  //
  // Annotation of the morse graph
  //
  int64_t const&annotation_index = morsegraph_record . annotation_index;
  Annotation_Record const& ar = database . annotationData () [annotation_index];
	// int64_t string_index = * ar . string_indices . begin ();
  std::vector<std::string> annotation_string;
  for ( int64_t string_index : ar . string_indices ) {
    annotation_string . push_back ( database . stringData () [ string_index ] );
  }
  result . annotation = annotation_string;
  //
  // Annotation of the morse graph vertices
  //
  std::vector<uint64_t> const& annotation_index_by_vertex =  morsegraph_record . annotation_index_by_vertex;
  std::vector<std::vector<std::string> > mystring;
  for ( uint64_t index : annotation_index_by_vertex ) {
    Annotation_Record const &arv = database . annotationData () [index];
		// int64_t string_index = * arv . string_indices . begin ();
    std::vector<std::string> newstring;
    for ( int64_t string_index : arv . string_indices ) {
      newstring . push_back ( database . stringData () [ string_index ] );
    }

    for ( std::string mys : newstring ) {
      std::vector<std::string> ss;
      ss . push_back  ( mys );
      mys = makeLabel(ss); // Write the proper annotations for the vertex
    }

    mystring . push_back ( newstring );
  }
  result . annotation_vertex = mystring;
  //
  //
  //
  uint64_t dag_index = morsegraph_record . dag_index;
  const DAG_Data & dag = database . dagData () [ dag_index ];
  // Vertices
  result . num_vertices_ = dag . num_vertices;
  result . labels_ . resize ( dag . num_vertices );
  for ( int i = 0; i < dag . num_vertices; ++ i ) {
    CS_Data cs_data;
    cs_data . vertices . push_back ( i );
    uint64_t cs_index = database . csIndex ( cs_data );
    INCCP_Record inccp_record;
    inccp_record . cs_index = cs_index;
    inccp_record . mgccp_index = mgccp_index;
    uint64_t inccp_index = database.inccpIndex ( inccp_record );
    uint64_t incc_index = database.inccp_to_incc () [ inccp_index ];
#ifndef NOCONLEYINDEX
    // uint64_t conley = database . incc_conley () [ incc_index ];
    // const CI_Data & ci = database . ciData () [ conley ];
    // result . labels_ [ i ] = conleyStringForZoo ( database, ci );
    result . labels_ [ i ] = std::string ( "" );
#else
    //    result . labels_ [ i ] = std::string ( "Unknown Conley Index" );
    result . labels_ [ i ] = std::string ( "" );
#endif
  }
 	// Edges
  typedef std::pair<int, int> Edge;
  for ( const Edge & e : dag . partial_order ) {
    DAG::Edge new_edge ( e . first, e . second );
    result . edges_ . insert ( new_edge );
  }
  return result;
}







int main ( int argc, char * argv [] ) {
  Database database;
  database . load ( argv [ 1 ] );

  const AbstractParameterSpace & space =
  dynamic_cast<const AbstractParameterSpace&> (database . parameter_space ());

  boost::unordered_set < uint64_t > parameters;
  boost::unordered_set < uint64_t > mgr_indices;

  const std::vector < MGCC_Record > & mgcc_records =
  database . MGCC_Records ();

  std::cout << "Number of records : " << mgcc_records.size() << "\n";

  Model model;
  model . initialize ( argc-1, argv+1 );
  //
  BooleanSwitchingParameterSpace & boolean_space = *
  boost::dynamic_pointer_cast<BooleanSwitchingParameterSpace> (
  model . parameterSpace () );

  typedef std::pair<uint64_t, uint64_t> Edge;
  std::vector<Edge> edges;

  // will store the nodes belonging to each mgcc
  std::vector< std::pair< uint64_t, std::vector<uint64_t> > > mgccNodes;

  // Extract a subpace of the parameter space
  std::unordered_set<uint64_t> subgraph;

  // Will store the different DAG class
  std::set<DAG> DAGclasses;
  std::set<DAG>::iterator itDAG;

  for ( uint64_t mgcc=0; mgcc<mgcc_records.size(); ++mgcc ) {

    std::vector<uint64_t> mynodes;
    int position;
    DAG myDAG = makeDAG ( database, mgcc );

    // check if the DAG is already in the list of classes
    itDAG = DAGclasses . find ( myDAG );
    //
    // if not, add it to the list of DAG classes
    // the position will be used as an index for the colormap
    if ( itDAG == DAGclasses . end() ) {
      DAGclasses . insert ( myDAG );
    }

    position = std::distance(DAGclasses.begin(), itDAG);

    // Loop through the MGCCPs
    for ( unsigned int mgccpi=0; mgccpi<mgcc_records[mgcc].mgccp_indices.size(); ++mgccpi ) {
      //
      const MGCCP_Record & mgccp_record =
      database . MGCCP_Records () [ mgcc_records[mgcc].mgccp_indices[mgccpi] ];
      //
      uint64_t morsegraph_index = mgccp_record . morsegraph_index;
      mgr_indices . insert ( morsegraph_index );
      const MorseGraphRecord & morsegraph_record =
      database . morsegraphData () [ morsegraph_index ];
      //
      uint64_t dag_index = morsegraph_record . dag_index;
      const DAG_Data & dag_data = database . dagData () [ dag_index ];
      //

      std::vector < uint64_t > pindex = mgccp_record . parameter_indices;
      //
      //
     for ( uint64_t p=0; p<pindex.size()-1; ++p ) {
        subgraph.insert(pindex[p]);
        mynodes.push_back(pindex[p]);
      }
      subgraph.insert(pindex[pindex.size()-1]);
      mynodes.push_back(pindex[pindex.size()-1]);
    }

    // we use the position for a colormap index
    mgccNodes . push_back ( std::pair<uint64_t,std::vector<uint64_t> >(position,mynodes) );

  }


  std::cout << "Number of DAG classes found : " << DAGclasses.size() << "\n";


  for ( uint64_t p : subgraph ) {
    std::vector<uint64_t> adj = boolean_space . adjacencies ( p );
    for ( uint64_t q : adj ) {
      if ( subgraph . count ( q ) ) {
        if ( p < q ) { edges . push_back ( std::make_pair ( p, q ) ); }
      }
    }
  }

  std::ofstream myfile;
  myfile . open ( "fullParameterGraph.gv" );
  myfile << "graph {\n";
  for ( std::pair< uint64_t,std::vector<uint64_t> > n : mgccNodes ) {
    for ( uint64_t value : n.second ) { // value is the index/label of the node
      myfile << value << "[shape=circle, style=filled,colorscheme=paired12, fillcolor="<< n.first+1 <<" ]\n";
    }
  }

 for ( Edge e : edges ) {
    myfile << e.first << " -- " << e.second <<"\n";
 }

// Add the legend : consider 12 colors maximum from paired12 colormap in graphviz
std::string colormap[] = {"#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6"};

myfile << "{ rank = sink;\n";
myfile << "Legend [shape=none, margin=0, label=<\n";
myfile << "<TABLE BORDER=\"0\" CELLBORDER=\"1\" CELLSPACING=\"0\" CELLPADDING=\"4\">\n";
myfile << "<TR>\n";
myfile << "<td> MGCC: </td>\n";
for ( uint64_t i=0; i<DAGclasses.size(); ++i ) {
  myfile << "<td bgcolor=\"" << colormap[i] << "\">" << i << "</td>\n";
}
myfile << "</TR>\n";
myfile << "</TABLE>\n";
myfile << ">];\n";
myfile << "}\n";

 myfile << "}";
 myfile.close();

  return 0;
}