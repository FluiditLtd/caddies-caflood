#ifndef _ARGSDATA_HPP_
#define _ARGSDATA_HPP_


//! \file ArgsData.hpp
//!  Contains the structure with the aruments data.
//! \author Michele Guidolin, University of Exeter, 
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2013-03

#include"BaseTypes.hpp"
#include"ca2D.hpp"
#include"Utilities.hpp"
#include<string>
#include<vector>

#define READ_TOKEN(_var,_tok,_sect)					\
  if(!CA::fromString(_var,_tok))					\
  {									\
    std::cerr<<"Error reading '"<<CA::trimToken(_sect)<<"' element"<<std::endl; \
    return 1;								\
  }


//! Remove file extension
inline std::string removeExtension(const std::string& filename) 
{
  size_t lastdot = filename.find_last_of(".");
  if (lastdot == std::string::npos) return filename;
  return filename.substr(0, lastdot); 
}


//! Public data structure that contains all the data set by the
//! arguments.
struct ArgsData
{

  // Create the arguments list and define the prefix which identifies
  // the arguments that are optional. This prefix differs between
  // windows and unix.
  CA::Arguments args;

  //! The working directory where the data and the configuration
  //! files are located.
  std::string working_dir;

  //! The setup file which contain the initial configuration of the ca
  //! algorithm.
  std::string setup_file;

  //! The output directory where the output data files are saved.
  std::string output_dir;

  //! Contain the string used as separator of directory.
  std::string sdir;

  //! The direcotry where temporary GRID/BUFFER data is saved/loaded.
  std::string data_dir;

  //! If true, print any information into console.
  bool info;		

  //! If true, perform preprocessing.
  bool pre_proc;		

  //! If true, perform postprocessing.
  bool post_proc;

  //! If true, perform WCA2D model.
  bool WCA2D;

  // Constructor
  ArgsData():
#if defined _WIN32 || defined __CYGWIN__   
    args("/"),
#else
    args("-"),
#endif
    working_dir("."),
#if defined _WIN32 || defined __CYGWIN__   
    sdir("\\"),
#else
    sdir("/"),
#endif
    data_dir(),
    info(false),
    pre_proc(false),
    post_proc(false),
    WCA2D(false)
  {}
  
  ~ArgsData(){}
  
};


//! Identifies different physical variables, i.e. the buffer of data.
struct PV
{
  enum Type
  {
    UNKNOWN = 0,
    WD,				//!< The water depth. 
    WL,				//!< The water level (depth + elv).
    VEL,			//!< The velocity saved as Vertical, Horizonatl.
  };
};

//! Transform an input string into a physical variable enum.
std::istream& operator>>(std::istream& in,  PV::Type& pv);
//! Transform a physical variable enum into an output string.
std::ostream& operator<<(std::ostream& out, PV::Type& pv);


//! The structure contains two temporary CellBuffer of the
//! velocity used to compute horizontal and vertical one.
struct VELBuff
{
  bool init;

  cpp11::shared_ptr<CA::CellBuffReal> V;//!< Vertical velocities cell buffer.
  cpp11::shared_ptr<CA::CellBuffReal> U;//!< Horizontal velocities cell buffer.

  CA::Real tupdated;		//!< The last time the buffer has been updated. 

  VELBuff():
    init(false),V(),U(), tupdated(0)
  {}

  ~VELBuff()
  {}
};


//! Initialise the temporary velocity cell buffer.
//! \return A non zero value if there was an error.
int initVELBuff(CA::Grid& GRID, double nodata, VELBuff& vb);


#endif
