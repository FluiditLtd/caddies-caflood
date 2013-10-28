#ifndef _TIMEPLOT_HPP_
#define _TIMEPLOT_HPP_


//! \file TimePlot.hpp
//! Contains the structure(s) of the various time plot(s) of the CA2D model.
//! \author Michele Guidolin, University of Exeter, 
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2013-03


#include"ca2D.hpp"
#include"BaseTypes.hpp"
#include"ArgsData.hpp"
#include"PointList.hpp"
#include<string>
#include<vector>
#include<fstream>


//! The configuration of output of a time plot of a physical variable
//! at a given coordinates/points using a given frequency.
struct TimePlot
{
  std::string    name;		        //!< Name of the plot. 
  PV::Type       pv;		        //!< The physical variable.
  std::vector<std::string> pnames;	//!< The list of points name.
  std::vector<CA::Real> xcoos;	        //!< The list of X coordinates.
  std::vector<CA::Real> ycoos;	        //!< The list of Y coordinates.
  CA::Real period;		        //!< The period in second.
};


//! The structure used during the model computation to store the time
//! plot data. 
struct TPData
{
  std::string   filename;	        //!< The name of the file. 
  cpp11::shared_ptr<std::ofstream> file;//!< The file where to output the time plot data.
  CA::PointList pl;	                //!< The coordinate of the points to plot.
  std::vector<CA::Real> pvals;	        //!< Buffer with the values of the points.
  std::vector<CA::Real> pelvs;	        //!< Eventual buffer with the values of the elevations.
  CA::Real      time_next;		//!< The time of the next output. 
};


//! Initialise the time plot structure usign a CSV file. 
//! Each row represents a new "variable" where the 
//! first column is the name of the element 
//! and the following columns have the multiple/single values.
//! \attention The order of elements is not important.
//! \param[in]  filename This is the file where the data is read.
//! \param[out] tp       The structure containing the read data.
//! \return A non zero value if there was an error.
int initTimePlotFromCSV(const std::string& filename, TimePlot& tp);


//! Initialise the time plot data that is used during the
//! computation from the time plot configuration.
//! \return A non zero value if there was an error.
int initTPData(const std::string& filename, const CA::Grid&  GRID, CA::CellBuffReal&  ELV,
	       const TimePlot& tp, TPData& tpdata);

#endif
