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

//! \file RasterGrid.cpp
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2013-03


#include"ArgsData.hpp"
#include"RasterGrid.hpp"
#include"Utilities.hpp"
#include<iostream>
#include<fstream>
#include<limits>


// -------------------------//
// Include the CA 2D functions //
// -------------------------//
#include CA_2D_INCLUDE(updatePEAKC)
#include CA_2D_INCLUDE(updatePEAKE)


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


RGManager::RGManager(CA::Grid&  GRID, const std::vector<RasterGrid>& rgs, 
		     const std::string& base, std::vector<std::string> names):
  _grid(GRID),
  _rgs(rgs),
  _datas(rgs.size()),
  _peak()
{
  for(size_t i = 0; i<_rgs.size(); ++i)
  {
    std::string filename(base+"_"+names[i]);
    initData(filename, _rgs[i], _datas[i], _peak);
  }
}
  

RGManager::~RGManager()
{

}


bool RGManager::updatePeak(const CA::BoxList&  domain, 
			   CA::CellBuffReal& WD, CA::CellBuffReal& V, CA::CellBuffState& MASK)
{
  // This variable make sure that the peak are updated only once.
  bool VAPEAKupdated = false;
  bool WDPEAKupdated = false;
  
  for(size_t i = 0; i<_datas.size(); ++i)
  {
    // Check if the peak values need to be updated.
    if(_rgs[i].peak == true)
    {
      switch(_rgs[i].pv)
      {
      case PV::VEL:
	// Update the absolute maximum velocity.
	// ATTENTION need to be tested.
	if(!VAPEAKupdated)
	{
	  CA::Execute::function(domain, updatePEAKC, _grid, (*_peak.V), V, MASK);
	  VAPEAKupdated = true;
	}
	// ATTENTION! The break is removed since in order to
	// post-process VA we need WD. 
	// break;	  
      case PV::WL:
      case PV::WD:	  
	// Update the absolute maximum water depth only once.
	if(!WDPEAKupdated)
	{
	  CA::Execute::function(domain, updatePEAKC, _grid, (*_peak.WD), WD, MASK);
	  WDPEAKupdated = true;
	}
	break;
      }
    }
  }

  return (WDPEAKupdated || VAPEAKupdated);
}


bool RGManager::outputPeak(CA::Real t, CA::CellBuffReal& WD, CA::CellBuffReal& V, 
			   const std::string& saveid,  bool output)
{
  // This variables is used to indicates if the output to console
  // happen in the case of time plot.
  bool outputed    = false;
  
  // Since WL variable does not exist. We need to save WD and then use
  // ELV to compute WL. The following variable is used to save WD only
  // once if both WD and WL are requested. The same is the case for
  // VEL, which needs WD.
  bool VAPEAKsaved   = false;
  bool WDPEAKsaved   = false;
  
  for(size_t i = 0; i<_datas.size(); ++i)
  {	
    // Check if the peak values need to be saved.
    if(_rgs[i].peak == true)
    {
      if(!outputed && output)
      {
	std::cout<<"Write Raster Grid (MIN "<<t/60<<"): ";
	outputed = true;
      }

      // Non velocity raster grid.
      switch(_rgs[i].pv)
      {
      case PV::VEL:	  
	if(!VAPEAKsaved)
	{
	  if(output)
	    std::cout<<" VAPEAK";
	  
	  _peak.V->saveData(saveid+"_V","PEAK");
	  VAPEAKsaved = true;
	}
	// ATTENTION! The break is removed since in order to
	// post-process VA we need WD. 
	// break;
      case PV::WL:
      case PV::WD:
	if(!WDPEAKsaved)
	{
	  if(output)
	    std::cout<<" WDPEAK";
	  
	  _peak.WD->saveData(saveid+"_WD","PEAK");
	  WDPEAKsaved = true;
	}
	break;
      }
      
    }
    // Update the next time to save a raster grid.
    _datas[i].time_next += _rgs[i].period;
  }

  if(outputed && output)
    std::cout<<std::endl;

  return (WDPEAKsaved || VAPEAKsaved);
}


bool RGManager::output(CA::Real t, CA::CellBuffReal& WD, 
		       CA::CellBuffReal& V, CA::CellBuffReal& A, 
		       const std::string& saveid,
		       bool output)
{
  // This variables is used to indicates if the output to console
  // happen in the case of time plot.
  bool outputed    = false;
  
  // Since WL variable does not exist. We need to save WD and then use
  // ELV to compute WL. The following variable is used to save WD only
  // once if both WD and WL are requested. The same is the case for
  // VEL, which needs WD.
  bool VAsaved       = false;
  bool VAPEAKsaved   = false;
  bool WDsaved       = false;
  bool WDPEAKsaved   = false;
  
  for(size_t i = 0; i<_datas.size(); ++i)
  {
    // Check if it is time to plot!
    if(t >= _datas[i].time_next)
    {	
      if(!outputed && output)
      {
	std::cout<<"Write Raster Grid (MIN "<<t/60<<"): ";
	outputed = true;
      }

      // Retrieve the string of the time.
      std::string strtime;
      CA::toString(strtime,std::floor(t+0.5));

      // Save the buffer using direct I/O where the main ID is the
      // buffer name and the subID is the timestep.
      switch(_rgs[i].pv)
      {
      case PV::VEL:	  
	if(!VAsaved)
	{
	  if(output)
	    std::cout<<" VA";	  
	  V.saveData(saveid+"_V",strtime);
	  A.saveData(saveid+"_A",strtime);
	  VAsaved = true;
	}
	// ATTENTION! The break is removed since in order to
	// post-process VA we need WD. 
	//break;
      case PV::WL:
      case PV::WD:
	if(!WDsaved)
	{
	  if(output)
	    std::cout<<" WD";
	  
	  WD.saveData(saveid+"_WD",strtime);
	  WDsaved = true;
	}
	break;
      }
	
      // Check if the peak values need to be saved.
      if(_rgs[i].peak == true)
      {
	// Non velocity raster grid.
	switch(_rgs[i].pv)
	{
	case PV::VEL:	  
	  if(!VAPEAKsaved)
	  {
	    if(output)
	      std::cout<<" VAPEAK";
	  
	    _peak.V->saveData(saveid+"_V","PEAK");
	    VAPEAKsaved = true;
	  }
	  // ATTENTION! The break is removed since in order to
	  // post-process VA we need WD. 
	  // break;
	case PV::WL:
	case PV::WD:
	  if(!WDPEAKsaved)
	  {
	    if(output)
	      std::cout<<" WDPEAK";
	  
	    _peak.WD->saveData(saveid+"_WD","PEAK");
	    WDPEAKsaved = true;
	  }
	  break;
	}
	  
      }
      // Update the next time to save a raster grid.
      _datas[i].time_next += _rgs[i].period;
    }
  }

  if(outputed && output)
    std::cout<<std::endl;

  return (WDPEAKsaved || VAPEAKsaved || WDsaved || VAsaved);
}


int RGManager::initData(const std::string& filename, const RasterGrid& rg, Data& rgdata, Peak& rgpeak)
{
  rgdata.filename = filename;
  
  switch(rg.pv)
  {
  case PV::VEL:	 
    if(rgpeak.V.get() == 0)
    {
      rgpeak.V.reset( new CA::CellBuffReal(_grid) ); 
      rgpeak.V->clear(0.0);
    }
    // ATTENTION! The break is removed since in order to
    // post-process VA we need WD. 
    //break;
  case PV::WD:
  case PV::WL:
    if(rgpeak.WD.get() == 0)
    {
      rgpeak.WD.reset( new CA::CellBuffReal(_grid) ); 
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
