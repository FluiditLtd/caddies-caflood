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

//! \file Rain.cpp
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2013-03


#include"ca2D.hpp"
#include"ArgsData.hpp"
#include"Rain.hpp"
#include"Utilities.hpp"
#include<iostream>
#include<fstream>

// -------------------------//
// Include the CA 2D functions //
// -------------------------//
#include CA_2D_INCLUDE(computeArea)
#include CA_2D_INCLUDE(computeVolume)
#include CA_2D_INCLUDE(addRain)
#include CA_2D_INCLUDE(addRainGrid)


// Initialise the RainEvents structure usign a CSV file. 
int initRainEventFromCSV(const std::string& filename, RainEvent& re, const std::string& data_dir)
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
      
      re.name = CA::trimToken(str," \t\r");
    }

    if(CA::compareCaseInsensitive("Rain Intensity",tokens[0],true))
    {
      found_tok=true;
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(found_tok,value,tokens[i],tokens[0]);

	re.rains.push_back(value);
      }
    }

    if(CA::compareCaseInsensitive("Time Stop",tokens[0],true))
    {
      found_tok=true;
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(found_tok,value,tokens[i],tokens[0]);

	re.times.push_back(value);
      }
    }

    if(CA::compareCaseInsensitive("Grid",tokens[0],true))
    {
      found_tok=true;
      for (size_t i=1; i<tokens.size(); ++i)
      {
          re.grids.push_back(data_dir + tokens[i]);
      }
    }

    if(CA::compareCaseInsensitive("Area",tokens[0],true))
    {
      found_tok=true;
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(found_tok,value,tokens[i],tokens[0]);

	re.area.push_back(value);
      }
    }

    if(CA::compareCaseInsensitive("Zone",tokens[0],true))
    {
      found_tok=true;
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(found_tok,value,tokens[i],tokens[0]);

	re.zone.push_back(value);
      }
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


RainManager::RainManager(CA::Grid&  GRID, const std::vector<RainEvent>& res):
  _grid(GRID),
  _res(res),
  _datas(res.size())
{
  for(size_t i = 0; i<_res.size(); ++i)
  {
    initData(_res[i], _datas[i]);
  }
}


RainManager::~RainManager()
{
  
}


//! Add the computational domain of the rain events into the given domain.
void RainManager::addDomain(CA::BoxList& compdomain)
{
  for(size_t i = 0; i<_datas.size(); ++i)
  {
    if (_datas[i].grid == nullptr)
        compdomain.add(_datas[i].box_area);
    else
        compdomain.add(_grid.box());
  }
}


void RainManager::analyseArea(CA::CellBuffReal& TMP, CA::CellBuffState& MASK, CA::BoxList&  domain)
{
  for(size_t i = 0; i<_datas.size(); ++i)
  {
    if (_datas[i].grid != nullptr)
        continue;

    TMP.fill(domain, 0.0);
    CA::Execute::function(_datas[i].box_area, computeArea, _grid, TMP, MASK);    
    TMP.sequentialOp(_datas[i].box_area, _datas[i].grid_area, CA::Seq::Add);
  }
}


void RainManager::prepare(CA::Real t, CA::Real period_time_dt, CA::Real next_dt)
{
  // Loop through the rain event(s).
  for(size_t i = 0; i<_res.size(); ++i)
  {
    // Compute the difference of rain betwenn the expected and the
    // added. This rain should be added/subtracted as one off when it
    // reaches high value.
    _datas[i].one_off_rain = (_datas[i].expected_rain - _datas[i].total_rain);

    // Set the rain and volume to zero.
    _datas[i].rain   = 0.0;
    _datas[i].volume = 0.0;

    // Get the index.
    size_t index = _datas[i].index;

    // Check if the simulation time now is equal or higher than the
    // time when this rain intensity ends. If it is the case,
    // increase the index to the next rain intensity.
    bool reload = false;
    if(t >= _res[i].times[index]) {
        index++;
        reload = true;
    }

    // If the index is larger than the available rain/time, do
    // nothing.
    if(index >= _res[i].rains.size() && index >= _res[i].grids.size())
      continue;

    // Load the next raster in, if demanded
    if (reload || _datas[i].grid == nullptr) {
        if (index < _res[i].grids.size() && !_res[i].grids[index].empty()) {
            if (_datas[i].grid == nullptr)
                _datas[i].grid = new CA::CellBuffReal(_grid);

            std::cout << "Loading rain raster " << _res[i].grids[index] << " for time " << t << std::endl;

            CA::ESRI_ASCIIGrid<CA::Real> rainGrid;
            rainGrid.readAsciiGrid(_res[i].grids[index]);

            CA::Box      realbox(_grid.box().x()+1,_grid.box().y()+1,_grid.box().w()-2,_grid.box().h()-2);
            _datas[i].grid->insertData(realbox, &rainGrid.data[0], rainGrid.ncols, rainGrid.nrows);
        }
    }

    // In case there are actual intensities present, use them.
    if (index < _res[i].rains.size()) {
        // Get the rain (transformed in metres from mm) for each dt of the next period.
        _datas[i].rain = (CA::Real)((_res[i].rains[index] * 0.001) * (next_dt / 3600.0));

        // Get the rain (transformed in metres from mm) for the next period.
        CA::Real period_rain = (CA::Real)((_res[i].rains[index] * 0.001) * (period_time_dt / 3600.0));

        // The expected amount of rain for the next period.
        _datas[i].expected_rain = period_rain;

        // Add the volume.
        _datas[i].volume = period_rain * _datas[i].grid_area;
    }
    // In case of grid source, let's just leave everything zero
    else {
        _datas[i].rain = 0.0;
        _datas[i].expected_rain = 0.0;
        _datas[i].volume = 0.0;
    }

    // Reset the total amount of rain for the next period.
    _datas[i].total_rain = 0.0;
    
    // Update index.
    _datas[i].index = index;
  }
}


CA::Real RainManager::volume()
{
  CA::Real rain_volume = 0.0;
  
  // Loop through the rain event(s).
  for(size_t i = 0; i<_res.size(); ++i)
    rain_volume  += _datas[i].volume;

  return rain_volume;
}


void RainManager::add(CA::CellBuffReal& WD, CA::CellBuffState& MASK, CA::Real t, CA::Real next_dt)
{
  // Loop through the rain event(s).
  for(size_t i = 0; i<_res.size(); ++i)
  {
    if (_datas[i].grid != nullptr) {
        CA::Real volume = 0.0;
        _datas[i].grid->sequentialOp(_grid.box(), volume, CA::Seq::Add);
        volume = (float)(volume * 0.001 * (next_dt / 3600.0) * _grid.area());

        CA::Real rainmultiplier = (CA::Real)(0.001 * (next_dt / 3600.0));
        CA::Execute::function(_grid.box(), addRainGrid, _grid, WD, MASK, *_datas[i].grid, rainmultiplier);

        _datas[i].volume += volume;
        _datas[i].expected_rain += volume;
        _datas[i].total_rain += volume;
    }

    // Do not add the rain if it is zero.
    else if(_datas[i].rain>=SMALL_RAIN)
    {
      // The amount of rain to add.
      CA::Real rain = static_cast<double>(_datas[i].rain)+_datas[i].one_off_rain;
      _datas[i].one_off_rain = 0.0; // Reset the one of rain.
	
      CA::Execute::function(_datas[i].box_area, addRain, _grid, WD, MASK, rain); 

      // Increse the amount of rain added into the period.
      _datas[i].total_rain += _datas[i].rain;
    }
  }
}


CA::Real RainManager::potentialVA(CA::Real t, CA::Real period_time_dt)
{
  CA::Real potential_va = 0.0;
  
  // Use the maximum amount of rain that is going to fall in the next
  // period to calculate the possible water depth.  Use the critical
  // velocity with this water depth to compute the potential velocity.
  
  // Loop through the rain event(s).
  for(size_t i = 0; i<_res.size(); ++i)
  {
    // Get the index.
    size_t index = _datas[i].index;

    // Check if the simulation time now is equal or higher than the
    // time when this rain intensity ends. If it is the case,
    // increase the index to the next rain intensity.
    if(t >= _res[i].times[index])
	index++;

    // If the index is larger than the available rain/time, do
    // nothing.
    if(index >= _res[i].rains.size() )
      continue;

    CA::Real rain;
    if (_datas[i].grid == nullptr) {
        // Get the rain (transformed in metres from mm) for the next period.
        rain = (_res[i].rains[index] * 0.001) * (period_time_dt / 3600.0);
    }
    else {
        _datas[i].grid->sequentialOp(_grid.box(), rain, CA::Seq::Max);
        rain = (float)(rain * 0.001 * (period_time_dt / 3600.0));
    }
    
    // Compute the potential velocity.
    potential_va = std::max(potential_va, std::sqrt(static_cast<CA::Real>(9.81)*( rain )) );
  }

  return potential_va;
}


CA::Real RainManager::endTime()
{
  CA::Real t_end = 0.0;

  // Loop through the rain event(s).
  for(size_t i = 0; i<_res.size(); ++i)
  {
    // Loop through the time steps.
    for(size_t j = 0; j<_res[i].times.size(); j++)
    {
      // Check if it is still raining at this time. If yest then
      // updated end_t.
      if (j < _res[i].grids.size())
         t_end = std::max(t_end, _res[i].times[j]);
      else if(j < _res[i].rains.size() && _res[i].rains[j]>0.0)
	     t_end = std::max(t_end,_res[i].times[j]);
    }
  }

  return t_end;
}


int RainManager::initData(const RainEvent& re, Data& data)
{
  data.index = 0;

  // If the area vector does not contain four values, then it is all
  // the domain that contain data.
  if(re.area.size()!=4 && re.zone.size()!=4)
  {
    data.box_area = _grid.box();
  }
  else
  {
    // Compute the given box area.
    if(re.area.size()==4)
    {
      data.box_area = CA::Box::create(_grid, re.area[0], re.area[1], re.area[2], re.area[3]);
    }
    // This is new. The are could have been identified as a zone, i.e. x,y,w,h/
    if(re.zone.size()==4)
    {
      // Compute the given box area from the zone.
      CA::Point     tl( CA::Point::create(_grid, re.zone[0], re.zone[1]) );
      CA::Unsigned  w = static_cast<CA::Unsigned>( std::ceil(re.zone[2] / _grid.length()) );
      CA::Unsigned  h = static_cast<CA::Unsigned>( std::ceil(re.zone[3] / _grid.length()) );
      
      data.box_area = CA::Box(tl.x(),tl.y(),w,h);
    }
  }

  data.grid_area = 0.0;
  data.volume    = 0.0;

  return 0;
}
