#ifndef _RASTERGRID_HPP_
#define _RASTERGRID_HPP_


//! \file RasterGrid.hpp
//! Contains the structure(s) of the various output raster grid(s) of the CA2D model.
//! \author Michele Guidolin, University of Exeter, 
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2013-03

#include"ca2D.hpp"
#include"BaseTypes.hpp"
#include"ArgsData.hpp"
#include<string>
#include<vector>


//! The configuration of the output of a raster grid of a physical
//! variable at a given coordinates/points using a given period.
struct RasterGrid
{
  std::string    name;		        //!< Name of the grid. 
  PV::Type       pv;		        //!< The physical variable.
  bool           peak;  		//!< If true, output the peak values of the physical variable. 
  CA::Real       period;		//!< The period in second.
};


//! The structure is used during the model computation to store the raster
//! grid data.
struct RGData
{
  std::string filename;			      //!< The base name of the file where the raster grid(s) are saved.
  CA::Real      time_next;		      //!< The time of the next output. 
};

//! The structure is used during the model computation to store the
//! peak raster grid data.
struct RGPeak
{
  cpp11::shared_ptr<CA::CellBuffReal> WD;  //!< Cell buffer with water depth peak values.
  cpp11::shared_ptr<CA::CellBuffReal> V;   //!< Cell buffer with velocity peak values.
};

//! Initialise the raster grid structure usign a CSV file. 
//! Each row represents a new "variable" where the 
//! first column is the name of the element 
//! and the following columns have the multiple/single values.
//! \attention The order of elements is not important.
//! \param[in]  filename This is the file where the data is read.
//! \param[out] setup    The structure containing the read data.
//! \return A non zero value if there was an error.
int initRasterGridFromCSV(const std::string& filename, RasterGrid& tp);


//! Initialise the raster grid data that is used during the
//! computation from the raster grid configuration.
//! \return A non zero value if there was an error.
int initRGData(const std::string& filename, CA::Grid& GRID, double nodata, const RasterGrid& rg, 
	       RGData& rgdata, RGPeak& rgpeak);


#endif
