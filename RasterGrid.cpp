//! \file RasterGrid.cpp
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2013-03


#include"ArgsData.hpp"
#include"RasterGrid.hpp"
#include"Utilities.hpp"
#include<iostream>
#include<fstream>
#include<limits>


// Initialise the RasterGrid structure usign a CSV file. 
int initRasterGridFromCSV(const std::string& filename, RasterGrid& rg)
{
  std::ifstream ifile(filename.c_str());
  
  if(!ifile)
  {
    std::cerr<<"Error opening CSV file: "<<filename<<std::endl;
    return 1;
  }
  
  // Parse the file line by line until the end of file 
  // and retrieve the tokens of each line.
  while(!ifile.eof())
  {
    std::vector<std::string> tokens( CA::getLineTokens(ifile, ',') );
    
    // If the tokens vector is empty we reached the eof or an
    // empty line... continue.
    if(tokens.empty())
      continue;       

    if(CA::compareCaseInsensitive("Raster Grid Name",tokens[0],true))
    {
      std::string str;
      READ_TOKEN(str,tokens[1],tokens[0]);
      
      rg.name = CA::trimToken(str);
    }

    if(CA::compareCaseInsensitive("Physical Variable",tokens[0],true))
      READ_TOKEN(rg.pv,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Peak",tokens[0],true))
    {
      std::string str( CA::trimToken(tokens[1]) );
      READ_TOKEN(rg.peak,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Period",tokens[0],true))
      READ_TOKEN(rg.period,tokens[1],tokens[0]);

  }

  return 0;
}


int initRGData(const std::string& filename, CA::Grid& GRID, double nodata, const RasterGrid& rg, 
	       RGData& rgdata, RGPeak& rgpeak)
{

  rgdata.filename = filename;
  
  switch(rg.pv)
  {
  case PV::VEL:	 
    if(rgpeak.V.get() == 0)
    {
      rgpeak.V.reset( new CA::CellBuffReal(GRID) ); 
      rgpeak.V->clear(0.0);
    }
    // ATTENTION! The break is removed since in order to
    // post-process VA we need WD. 
    //break;
  case PV::WD:
  case PV::WL:
    if(rgpeak.WD.get() == 0)
    {
      rgpeak.WD.reset( new CA::CellBuffReal(GRID) ); 
      rgpeak.WD->clear(0.0);
    }
    break;
  }

  if(rg.period > 0.0)
    rgdata.time_next = rg.period;
  else
    rgdata.time_next = std::numeric_limits<CA::Real>::max(); 

  return 0;
}
