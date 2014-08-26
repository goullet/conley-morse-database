#ifndef BOOLEANSWITCHINGPARAMETERSPACE_H
#define BOOLEANSWITCHINGPARAMETERSPACE_H

#include <iostream>
#include <sstream>
#include <string>
#include <exception>
#include <vector>
#include <stack>
#include <fstream>
#include <sstream>
#include <unordered_map>

#include <boost/shared_ptr.hpp>

#include "database/structures/ParameterSpace.h"
#include "database/structures/AbstractParameterSpace.h"

#include "Network.h"
#include "Parameter/FactorGraph.h"
#include "Parameter/BooleanSwitchingParameter.h"
#include "Parameter/Polytope.h"

/// class BooleanSwitchingParameterSpace
class BooleanSwitchingParameterSpace : public AbstractParameterSpace {
public:
  // typedef
  typedef uint64_t ParameterIndex;
  typedef boost::counting_iterator < ParameterIndex > iterator;
  typedef iterator const_iterator;
  // Constructor/Deconstructor
  BooleanSwitchingParameterSpace ( void ) {}
  virtual ~BooleanSwitchingParameterSpace ( void ) {}

  /// initialize
  ///    Create the ParameterSpace given the configuration specified
  virtual void initialize ( int argc, char * argv [] );

  /// adjacencies
  ///    Return a vector of adjacent vertices.
  virtual std::vector<ParameterIndex> adjacencies ( ParameterIndex v ) const;
  
  /// size
  ///    Return the number of vertices
  virtual uint64_t size ( void ) const;

  /// parameter
  ///    Return the parameter object associated with a vertex
  virtual boost::shared_ptr<Parameter> parameter ( ParameterIndex v ) const;
  
  /// search
  ///    Given a parameter, find the vertex associated with it
  ///    (This can be used to find a parameter which might contain the other)
  virtual uint64_t search ( boost::shared_ptr<Parameter> parameter ) const;
  
  /// closestFace
  ///   given a domain, return the closest face 
  ///   closest face is output in the following form:
  ///   There are d entries in an std::vector<int>
  ///     0 means lower bound, 1 means between, 2 means upper bound
  std::vector<int> closestFace ( boost::shared_ptr<Parameter> parameter, 
                                 std::vector<size_t> const& domain ) const;
  /// domainLimits
  ///    Return a vector containing the number of thresholds plus one in each dimension
  ///    This gives us the number of bins in each dimension, which is needed
  ///    for multidimensional iteration through the "domains" (i.e. regions between
  ///    thresholds)
  std::vector<size_t> domainLimits ( void ) const;

  /// dimension
  ///    Return dimension
  int dimension ( void ) const;

  /// factorGraph
  ///    Return factor graph
  const FactorGraph &
  factorGraph ( int i ) const;

  /// prettyPrint
  ///    Print out parameter in human readable format
  std::string
  prettyPrint ( boost::shared_ptr<Parameter> parameter ) const;

private:

  /// dimension_
  ///   dimension of phase space. This is equal to the number of network nodes.
  int dimension_; 

  /// factors_
  ///    The stored FactorGraphs
  std::vector<FactorGraph> factors_;

  /// network_
  BooleanSwitching::Network network_;

  // Serialization
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
  }
};

inline void 
BooleanSwitchingParameterSpace::initialize ( int argc, char * argv [] ) {
  std::cout << "BooleanSwitchingParameterSpace::initialize\n";
  // Load the network file
  std::string filestring ( argv[1] );
  std::string appendstring ( argv[2] );
  std::string loadstring = filestring + appendstring;

  std::cout << "BooleanSwitchingParameterSpace::initialize." << 
    " Loading network file " << loadstring << "\n";
    network_ . load ( loadstring.c_str() );
    // Get dimension
  dimension_ = network_ . size ();
  std::cout << "BooleanSwitchingParameterSpace::initialize." << 
    "dimension_ = " << dimension_ << "\n"; // DEBUG

  // Loop through nodes and create FactorGraphs
  factors_ . resize ( dimension_ );
  for ( BooleanSwitching::Node const& node : network_ ) {
    int d = node . index - 1; // Index of node (minus one to start indexing at 0)
    int n = 0; // Number of in edges for node d
    std::vector<int> logic; // Logic for node d inputs
    for ( int i = 0; i < node . logic . size (); ++ i ) {
        int k = node . logic [ i ] . size ();
        logic . push_back ( k );
            n += k;
        }
    int m = node . out_order . size ();
    factors_ [ d ] . construct ( MonotonicMap ( n, m, logic ) );
    std::cout << "BooleanSwitchingParameterSpace::initialize." << 
      "factors_[" << d << "].size() = " << factors_[d].size() << "\n"; // DEBUG
  }
}

inline std::vector<BooleanSwitchingParameterSpace::ParameterIndex> 
BooleanSwitchingParameterSpace::adjacencies ( ParameterIndex v ) const {
  std::vector<ParameterIndex> result;
  boost::shared_ptr<BooleanSwitchingParameter> p = 
      boost::dynamic_pointer_cast<BooleanSwitchingParameter> ( parameter ( v ) );
  if ( not p ) {
    std::stringstream ss;
    ss << "BooleanSwitchingParameterSpace::adjacencies. ";
    ss << "Invalid ParameterIndex v = " << v << "\n";
    throw std::domain_error ( ss . str () );
  }
  // Loop through coordinates and change monotonic functions by one
  uint64_t multiplier = 1;
  for ( int d = 0; d < dimension_; ++ d ) {
    int digit = p -> monotonic_function_ [ d ];
    const std::vector<int> & neighbors = factors_ [ d ] . adjacencies ( digit );
    for ( int neighbor : neighbors ) {
        result . push_back ( v + multiplier * ( neighbor - digit ) );
    }
    multiplier *= factors_ [ d ] . size ();
  }
  return result;
}
    
inline uint64_t 
BooleanSwitchingParameterSpace::size ( void ) const {
  uint64_t result = 1;
  for ( int d = 0; d < dimension_; ++ d ) {
    result *= factors_ [ d ] . size ();
  }
  return result;
}

inline boost::shared_ptr<Parameter> 
BooleanSwitchingParameterSpace::parameter ( ParameterIndex v ) const {
  boost::shared_ptr<BooleanSwitchingParameter> 
    p ( new BooleanSwitchingParameter(dimension_) );
  for ( int d = 0; d < dimension_; ++ d ) {
    size_t factor_size = factors_ [ d ] . size ();
    p -> monotonic_function_ [ d ] = v % factor_size;
    v /= factor_size;
  }
  return boost::dynamic_pointer_cast<Parameter> ( p );
}
    
inline uint64_t 
BooleanSwitchingParameterSpace::search ( boost::shared_ptr<Parameter> parameter ) const {
  const BooleanSwitchingParameter & p = 
    * boost::dynamic_pointer_cast<BooleanSwitchingParameter> ( parameter );
  uint64_t result = 0;
  uint64_t multiplier = 1;
  for ( int d = 0; d < dimension_; ++ d ) {
    result += multiplier * p . monotonic_function_ [ d ];
    multiplier *= factors_ [ d ] . size ();
  }
  return result;
}


inline std::vector<int> 
BooleanSwitchingParameterSpace::closestFace 
                ( boost::shared_ptr<Parameter> p, 
                  const std::vector<size_t> const& domain ) const {
  boost::shared_ptr<BooleanSwitchingParameter> parameter =
  boost::dynamic_pointer_cast<BooleanSwitchingParameter> ( p ); 
  std::vector<int> result ( dimension_ );
  if ( dimension_ != domain . size () ) { 
    std::cout << "error. BooleanSwitchingParameter::closestFace. "
                 "Inappropriate input domain size.\n";
    throw std::logic_error ( "BooleanSwitchingParameter::closestFace. "
                             "Inappropriate input domain size.\n");
  }
  // Determination of out-states.
  //    Loop through each node of network.
  //    For each node, examine the out-edges in the order they are listed
  //    Depending on the "critical_value", which is the domain bin
  //    of the associated phase space variable, we determine whether this
  //    output should be "on" or "off". We store this data in the 
  //    "state" data structure, which maps (source,target) to "off,on"
  //  Notes:  The indexing of nodes starts at 1, whereas the indexing
  //          of dimension starts at 0. This requires a translation
  //          from network indexing to dimension indexing, which is
  //          just to subtract 1.
  //          We store the "state" via network indexing.
  std::unordered_map < std::pair<int, int>, bool > state;
  for ( BooleanSwitching::Node const& node : network_ ) {
    int critical_value = domain [ node . index - 1 ];
    int count = 0;
    for ( int out_node : node . out_order ) {
      state [ std::make_pair ( node . index, out_node ) ] = 
        ( count ++ < critical_value );
    }
  }
  // Loop through each node of network
  //  Loop through each in-edge (i.e. each factor of logic, 
  //                                  each summand of factor)
  //     Access state corresponding to in-edge.
  //     If this is a down-regulator, flip the state
  //     Append the state as the last bit of a growing code-word
  //  End
  //  Determine monotonic function for node
  //  Call the monotonic function with the code look-up
  //  Write the result component
  // End
  for ( BooleanSwitching::Node const& node : network_ ) {
    uint64_t code = 0;
    for ( std::vector<int> const& factor : node . logic ) {
      for ( int in_node : factor ) {
        code <<= 1;
        bool bit = state [ std::make_pair ( std::abs(in_node), 
                           node . index ) ];
        if ( in_node < 0 ) bit = not bit; // Take into account down-regulation
        if ( bit ) ++ code;
      }
    }
    // Note. Input code for node (node . index) has been established
    int d = node . index - 1; // Index of node (minus one to start indexing at 0)
    int monotonic_function_index = parameter -> monotonic_function_ [ d ];
    const MonotonicMap & monotonic_function = 
      factors_ [ d ] . vertices [ monotonic_function_index ];
    int bin = monotonic_function . data_ [ code ];
    if ( bin < domain [ d ] ) result [ d ] = 0;
    else if ( bin == domain [ d ] ) result [ d ] = 1;
    else if ( bin > domain [ d ] ) result [ d ] = 2;    
  }
  return result;
}

inline std::vector<size_t> 
BooleanSwitchingParameterSpace::domainLimits ( void ) const {
  std::vector<size_t> result ( dimension_ );
  for ( BooleanSwitching::Node const& node : network_ ) {
    result [ node . index - 1 ] = node . out_order . size () + 1;
  }
  return result;
}

inline int
BooleanSwitchingParameterSpace::dimension ( void ) const {
  return dimension_;
}

inline const FactorGraph &
BooleanSwitchingParameterSpace::factorGraph ( int i ) const {
  return factors_ [ i ];
}

std::string BooleanSwitchingParameterSpace::
prettyPrint ( boost::shared_ptr<Parameter> parameter ) const {
  std::stringstream result;
  const BooleanSwitchingParameter & p = 
      * boost::dynamic_pointer_cast<BooleanSwitchingParameter> ( parameter );
  for ( int d = 0; d < dimension_; ++ d ) {
    std::string symbol = network_ . name ( d );
    std::vector<std::string> input_symbols, output_symbols;
    BooleanSwitching::Node node = network_ . node ( d + 1 );
    for ( std::vector<int> const& factor : node . logic ) {
      for ( int input : factor ) {
        input_symbols . push_back ( network_ . name ( std::abs(input) );
      }
    }
    for ( int output : node . out_order ) {
      output_symbols . push_back ( network_ . name ( output );
    }
    MonotonicMap mono = factors_ [ d ] . vertices [ p . monotonic_function_ [ d ] ];
    result << mono . prettyPrint ( symbol, input_symbols, output_symbols );
  }
  return result . str ();  
}


#endif
