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

//! \file WaterLevel.cpp
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2013-03


#include"ArgsData.hpp"
#include"WaterLevel.hpp"
#include"Utilities.hpp"
#include<iostream>
#include<fstream>

// -------------------------//
// Include the CA 2D functions //
// -------------------------//
#include CA_2D_INCLUDE(computeArea)
#include CA_2D_INCLUDE(addRaise)


int initWLEventFromCSV(const std::string& filename, WLEvent& wle)
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
    // If true the token was identified;
    bool found_tok = false;

    std::vector<std::string> tokens( CA::getLineTokens(ifile, ',') );
    
    // If the tokens vector is empty we reached the eof or an
    // empty line... continue.
    if(tokens.empty())
      continue;       

    if(CA::compareCaseInsensitive("Event Name",tokens[0],true))
    {
      std::string str;
      READ_TOKEN(found_tok,str,tokens[1],tokens[0]);
      
      wle.name = CA::trimToken(str," \t\r");
    }

    if(CA::compareCaseInsensitive("Water Level",tokens[0],true))
    {
      found_tok=true;
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(found_tok,value,tokens[i],tokens[0]);

	wle.wls.push_back(value);
      }
    }

    if(CA::compareCaseInsensitive("Time",tokens[0],true))
    {
      found_tok=true;
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(found_tok,value,tokens[i],tokens[0]);

	wle.times.push_back(value);
      }
    }

    if(CA::compareCaseInsensitive("Area",tokens[0],true))
    {
      found_tok=true;
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(found_tok,value,tokens[i],tokens[0]);

	wle.area.push_back(value);
      }
    }

    if(CA::compareCaseInsensitive("Zone",tokens[0],true))
    {
      found_tok=true;
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(found_tok,value,tokens[i],tokens[0]);

	wle.zone.push_back(value);
      }
    }

    if(CA::compareCaseInsensitive("Analytical Solution U",tokens[0],true))
    {
      found_tok=true;
      READ_TOKEN(found_tok,wle.u,tokens[1],tokens[0]);
    }

    if(CA::compareCaseInsensitive("Analytical Solution N",tokens[0],true))
    {
      found_tok=true;
      READ_TOKEN(found_tok,wle.n,tokens[1],tokens[0]);
    }

    // If the token was not identified stop!
    if(!found_tok)
    {
      std::cerr<<"Element '"<<CA::trimToken(tokens[0])<<"' not identified"<<std::endl; \
      return 1;
    }
  }
  
  return 0;
}


WaterLevelManager::WaterLevelManager(CA::Grid&  GRID, const std::vector<WLEvent>& wles):
  _grid(GRID),
  _wles(wles),
  _datas(wles.size())
{
  for(size_t i = 0; i<_wles.size(); ++i)
  {
    initData(_wles[i], _datas[i]);
  }
}


WaterLevelManager::~WaterLevelManager()
{
  
}


//! Add the computational domain of the WaterLevel events into the given domain.
void WaterLevelManager::addDomain(CA::BoxList& compdomain)
{
  for(size_t i = 0; i<_datas.size(); ++i)
  {
    compdomain.add(_datas[i].box_area);
  }
}


void WaterLevelManager::analyseArea(CA::CellBuffReal& TMP, CA::CellBuffState& MASK, CA::BoxList&  domain)
{
  for(size_t i = 0; i<_datas.size(); ++i)
  {
    TMP.fill(domain, 0.0);
    CA::Execute::function(_datas[i].box_area, computeArea, _grid, TMP, MASK);    
    TMP.sequentialOp(_datas[i].box_area, _datas[i].grid_area, CA::Seq::Add); 
  }
}


void WaterLevelManager::getElevation(CA::CellBuffReal& ELV)
{
  // Loop through the WaterLevel event(s).
  for(size_t i = 0; i<_wles.size(); ++i)
  {
    // Retrieve the minimum elevation of the given area.
    ELV.sequentialOp(_datas[i].box_area, _datas[i].last_level, CA::Seq::MinAbs);	
  }
}


void WaterLevelManager::prepare(CA::Real t, CA::Real period_time_dt, CA::Real next_dt)
{
  // Loop through the WaterLevel event(s).
  for(size_t i = 0; i<_wles.size(); ++i)
  {
    // Set the volume to zero.
    _datas[i].volume = 0.0;    

    size_t index = _datas[i].index;

    // If the index is larger than the available rain/time, do
    // nothing.
    if(index >= _wles[i].wls.size() )
      continue;

    // Compute the water level at specific are using
    // interpolation. Check if the index is the last available
    // one. In this case use only one value.
    CA::Real level      = 0;
    if(index == _wles[i].wls.size() -1)
    {
      level = _wles[i].wls[index];
    }
    else
    {	
      CA::Real y0 = _wles[i].wls[index];
      CA::Real y1 = _wles[i].wls[index+1];
      CA::Real x0 = _wles[i].times[index];
      CA::Real x1 = _wles[i].times[index+1];
      level = y0 + (y1-y0) * ( ((t+period_time_dt) - x0)/(x1 - x0) );
    }

    // ATTENTION At the moment to work the last level must be the original elevation.
    //_datas[i].last_level = level; 

  }
}


CA::Real WaterLevelManager::volume()
{
  CA::Real wl_volume = 0.0;
  
  // Loop through the water level event(s).
  for(size_t i = 0; i<_wles.size(); ++i)
    wl_volume  += _datas[i].volume;

  return wl_volume;
}


void WaterLevelManager::add(CA::CellBuffReal& WD, CA::CellBuffReal& ELV,
			    CA::CellBuffState& MASK, CA::Real t, CA::Real next_dt)
{
  // Loop through the WaterLevel event(s).
  for(size_t i = 0; i<_wles.size(); ++i)
  {
    // If the value of U is set (differ from zero). Then the
    // analytical solution is compute with C=0
    if(_wles[i].u!=0)
    {
      CA::Real pn = std::pow(_wles[i].n,2);
      CA::Real pu = std::pow(_wles[i].u,3);
      CA::Real level = std::pow((7.0/3.0)*(0-pn*pu*(-_wles[i].u*(t+next_dt))) ,(3.0/7.0));
      // Given the way the CA2D model work, we need to set the water
      // depth instead of the water level. Thus the water depth value
      // at specific location is the value of the water level event
      // minus the elevation.
      CA::Execute::function(_datas[i].box_area, addRaise, _grid, WD, ELV, MASK, level);     
      continue;
    }

    size_t index = _datas[i].index;

    // If the index is larger than the available rain/time, do
    // nothing.
    if(index >= _wles[i].wls.size() )
      continue;
    
    // Compute the water level at specific are using
    // interpolation. Check if the index is the last available
    // one. In this case use only one value.
    CA::Real level      = 0;
    if(index == _wles[i].wls.size() -1)
    {
      level = _wles[i].wls[index];
    }
    else
    {	
      CA::Real y0 = _wles[i].wls[index];
      CA::Real y1 = _wles[i].wls[index+1];
      CA::Real x0 = _wles[i].times[index];
      CA::Real x1 = _wles[i].times[index+1];
      level = y0 + (y1-y0) * ( (t - x0)/(x1 - x0) );
    }
        
    // Given the way the CA2D model work, we need to set the water
    // depth instead of the water level. Thus the water depth value
    // at specific location is the value of the water level event
    // minus the elevation.
    CA::Execute::function(_datas[i].box_area, addRaise, _grid, WD, ELV, MASK, level);     
    
    // Check if the simulation time now is equal or higher than the
    // time of the NEXT index.
    if(t >= _wles[i].times[index+1])
      index++;
    
    // Update index.
    _datas[i].index = index;
  }
}


CA::Real WaterLevelManager::potentialVA(CA::Real t, CA::Real period_time_dt)
{
  CA::Real potential_va = 0.0;
  
  // Use the differnce in water level between the previous water level
  // set and the new one that is going to happen. If it is the first
  // call then the previous level is the elevation value.
  
  // Loop through the water level event(s).
  for(size_t i = 0; i<_wles.size(); ++i)
  {
    size_t index = _datas[i].index;
    
    // If the index is larger than the available time, do
    // nothing.
    if(index >= _wles[i].wls.size() )
      continue;
    
    // Compute the water level at specific are using
    // interpolation. Check if the index is the last available
    // one. In this case use only one value.
    CA::Real level      = 0;
    if(index == _wles[i].wls.size() -1)
    {
      level = _wles[i].wls[index];
    }
    else
    {	
      CA::Real y0 = _wles[i].wls[index];
      CA::Real y1 = _wles[i].wls[index+1];
      CA::Real x0 = _wles[i].times[index];
      CA::Real x1 = _wles[i].times[index+1];
      level = y0 + (y1-y0) * ( ((t+period_time_dt) - x0)/(x1 - x0) );
    }
    
    // Retrieve the difference from the previous water level and use
    // diss difference to compute the critical velocity.
    CA::Real  water_diff = std::abs(level - _datas[i].last_level);
    potential_va = std::max(potential_va, std::sqrt( static_cast<CA::Real>(9.81)*( water_diff) ) );
    
    
  }
 
  // ATTENTION! This does not need to be precise but just give a rough estimation
  return potential_va;
}


CA::Real WaterLevelManager::endTime()
{
  CA::Real t_end = 0.0;

  // Loop through the water level event(s).
  for(size_t i = 0; i<_wles.size(); ++i)
  {
    t_end = std::max(t_end,_wles[i].times.back());
  }

  return t_end;
}


int WaterLevelManager::initData(const WLEvent& wle, Data& data)
{
  data.index = 0;

  // If the area vector does not contain four values, then the box
  // stay empty. NO WATER LEVEL.
  if(wle.area.size()==4)
  {
    // Compute the given box area.
    data.box_area = CA::Box::create(_grid, wle.area[0], wle.area[1], wle.area[2], wle.area[3]);
  }
  // This is new. The area could have been identified as a zone, i.e. x,y,w,h/
  if(wle.zone.size()==4)
  {
    // Compute the given box area from the zone.
    CA::Point     tl( CA::Point::create(_grid, wle.zone[0], wle.zone[1]) );
    CA::Unsigned  w = static_cast<CA::Unsigned>( std::ceil(wle.zone[2] / _grid.length()) );
    CA::Unsigned  h = static_cast<CA::Unsigned>( std::ceil(wle.zone[3] / _grid.length()) );
    
    data.box_area = CA::Box(tl.x(),tl.y(),w,h);
  }


  return 0;
}
