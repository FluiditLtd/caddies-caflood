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
#include CA_2D_INCLUDE(addRain)


// Initialise the RainEvents structure usign a CSV file. 
int initRainEventFromCSV(const std::string& filename, RainEvent& re)
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

    if(CA::compareCaseInsensitive("Event Name",tokens[0],true))
    {
      std::string str;
      READ_TOKEN(str,tokens[1],tokens[0]);
      
      re.name = CA::trimToken(str," \t\r");
    }

    if(CA::compareCaseInsensitive("Rain Intensity",tokens[0],true))
    {
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(value,tokens[i],tokens[0]);

	re.rains.push_back(value);
      }
    }

    if(CA::compareCaseInsensitive("Time Stop",tokens[0],true))
    {
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(value,tokens[i],tokens[0]);

	re.times.push_back(value);
      }
    }

    if(CA::compareCaseInsensitive("Area",tokens[0],true))
    {
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(value,tokens[i],tokens[0]);

	re.area.push_back(value);
      }
    }

    if(CA::compareCaseInsensitive("Zone",tokens[0],true))
    {
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(value,tokens[i],tokens[0]);

	re.zone.push_back(value);
      }
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
    compdomain.add(_datas[i].box_area);
  }
}


void RainManager::analyseArea(CA::CellBuffReal& TMP, CA::CellBuffState& MASK, CA::BoxList&  domain)
{
  for(size_t i = 0; i<_datas.size(); ++i)
  {
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
    // Set the rain and volume to zero.
    _datas[i].rain   = 0.0;
    _datas[i].volume = 0.0;

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

    // Get the rain (transformed in metres from mm) for each dt of the next period.
    _datas[i].rain = (_res[i].rains[index] * 0.001) * (next_dt/3600.0);

    // Get the rain (transformed in metres from mm) for the next period.
    CA::Real period_rain = (_res[i].rains[index] * 0.001) * (period_time_dt/3600.0);

    // Add the volume.
    _datas[i].volume = period_rain * _datas[i].grid_area;
    
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


void RainManager::add(CA::CellBuffReal& WD, CA::CellBuffState& MASK)
{
  // Loop through the rain event(s).
  for(size_t i = 0; i<_res.size(); ++i)
  {
    // Do not add the rain if it is zero.
    if(_datas[i].rain>=SMALL_RAIN)
    {
      CA::Execute::function(_datas[i].box_area, addRain, _grid, WD, MASK, _datas[i].rain); 
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

    // Get the rain (transformed in metres from mm) for the next period.
    CA::Real rain = (_res[i].rains[index] * 0.001) * (period_time_dt/3600.0);
    
    // Compute the potential velocity.
    potential_va = std::max(potential_va, std::sqrt(static_cast<CA::Real>(9.81)*( rain )) );
  }

  return potential_va;
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
