/*
    
Copyright (c) 2013 Centre for Water Systems,
                   University of Exeter

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*/

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
  bool           final;  		//!< If true, output the final extend values of the physical variable. 
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



//! Class that manage all Raster grid outputs
class RGManager
{
private:
  
  //! The structure used during the model computation to store the raster grid data.
  struct Data
  {
    std::string filename;            //!< The base name of the file where the raster grid(s) are saved.
    CA::Real      time_next;	     //!< The time of the next output. 
  };

  //! The structure is used during the model computation to store the
  //! peak raster grid data.
  struct Peak
  {
    cpp11::shared_ptr<CA::CellBuffReal> WD;  //!< Cell buffer with water depth peak values.
    cpp11::shared_ptr<CA::CellBuffReal> V;   //!< Cell buffer with velocity peak values.
  };
  

public:

  //! Construct a Raster Grid manager
  //! \param base  This is the base for all the output filenames of the various raster grid.
  //! \param names This is a list of the names for the raster grid outut files.
  RGManager(CA::Grid&  GRID, const std::vector<RasterGrid>& rgs, 
	    const std::string& base, std::vector<std::string> names);

  //! Destroy a Raster Grid Manager.
  ~RGManager();
  
  
  //! Update the peak values.
  //! \param  domain     The are to update the peak.
  //! \params WD         The cell buffer with the water depth.
  //! \params V          The cell buffer with the velocity magnitude.
  //! \params MASK       The cell buffer with the mask
  //! \return True if the peak were updated.
  bool updatePeak(const CA::BoxList&  domain, CA::CellBuffReal& WD, CA::CellBuffReal& V, CA::CellBuffState& MASK);


  //! Output only the peak raster grids
  //! \params t          The simulation time.
  //! \params WD         The cell buffer with the water depth.
  //! \params V          The cell buffer with the velocity magnitude.
  //! \params saveid     The id to use to save the buffers.
  //! \params output     If true, output information to console.
  //! \return True if the rasters were outputed.
  bool outputPeak(CA::Real t,CA::CellBuffReal& WD, CA::CellBuffReal& V, const std::string& saveid, bool output);

  
  //! Output all the raster grids 
  //! \params t          The simulation time.
  //! \params WD         The cell buffer with the water depth.
  //! \params V          The cell buffer with the velocity magnitude.
  //! \params A          The cell buffer with the velocity angle.
  //! \params saveid     The id to use to save the buffers.
  //! \params output     If true, output information to console.
  //! \param  final      If true, this is the final iteration.
  //! \return True if the rasters were outputed.
  bool output(CA::Real t, CA::CellBuffReal& WD, CA::CellBuffReal& V, CA::CellBuffReal& A,
	      const std::string& saveid, bool output, bool final = false);
  

protected:


  //! Initialise the raster grid data that is used during the
  //! computation from the raster grid configuration.
  int initData(const std::string& filename, const RasterGrid& rg, Data& rgdata, Peak& rgpeak);

private:

  //! Reference to the grid.
  CA::Grid& _grid;

  //! Reference to the List of raster grids
  const std::vector<RasterGrid>& _rgs;

  //! List of raster grid data.
  std::vector<Data> _datas;

  // Peak buffers
  Peak _peak;

};


#endif
