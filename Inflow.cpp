//! \file Inflow.cpp
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2013-03


#include"ca2D.hpp"
#include"ArgsData.hpp"
#include"Inflow.hpp"
#include"Utilities.hpp"
#include<iostream>
#include<fstream>

// -------------------------//
// Include the CA 2D functions //
// -------------------------//
#include CA_2D_INCLUDE(computeArea)


int initIEventFromCSV(const std::string& filename, IEvent& ie)
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
      
      ie.name = CA::trimToken(str," \t\r");
    }

    if(CA::compareCaseInsensitive("Inflow",tokens[0],true))
    {
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(value,tokens[i],tokens[0]);

	ie.ins.push_back(value);
      }
    }

    if(CA::compareCaseInsensitive("Time",tokens[0],true))
    {
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(value,tokens[i],tokens[0]);

	ie.times.push_back(value);
      }
    }

    if(CA::compareCaseInsensitive("Area",tokens[0],true))
    {
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(value,tokens[i],tokens[0]);

	ie.area.push_back(value);
      }
    }

    if(CA::compareCaseInsensitive("Zone",tokens[0],true))
    {
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(value,tokens[i],tokens[0]);

	ie.zone.push_back(value);
      }
    }
  }
  
  return 0;
}


InflowManager::InflowManager(CA::Grid&  GRID, const std::vector<IEvent>& ies):
  _grid(GRID),
  _ies(ies),
  _datas(ies.size())
{
  for(size_t i = 0; i<_ies.size(); ++i)
  {
    initData(_ies[i], _datas[i]);
  }
}


InflowManager::~InflowManager()
{
  
}


//! Add the computational domain of the inflow events into the given domain.
void InflowManager::addDomain(CA::BoxList& compdomain)
{
  for(size_t i = 0; i<_datas.size(); ++i)
  {
    compdomain.add(_datas[i].box_area);
  }
}


void InflowManager::analyseArea(CA::CellBuffReal& TMP, CA::CellBuffState& MASK, CA::BoxList&  domain)
{
  for(size_t i = 0; i<_datas.size(); ++i)
  {
    TMP.fill(domain, 0.0);
    CA::Execute::function(_datas[i].box_area, computeArea, _grid, TMP, MASK);    
    TMP.sequentialOp(_datas[i].box_area, _datas[i].grid_area, CA::Seq::Add);
  }
}


void InflowManager::prepare(CA::Real t, CA::Real period_time_dt, CA::Real next_dt)
{

}


CA::Real InflowManager::volume()
{
  CA::Real inflow_volume = 0.0;
  
  // Loop through the inflow event(s).
  for(size_t i = 0; i<_ies.size(); ++i)
    inflow_volume  += _datas[i].volume;

  return inflow_volume;
}


void InflowManager::add(CA::CellBuffReal& WD, CA::CellBuffState& MASK)
{

}


CA::Real InflowManager::potentialVA(CA::Real t, CA::Real period_time_dt)
{
  CA::Real potential_va = 0.0;
  

  return potential_va;
}


int InflowManager::initData(const IEvent& ie, Data& data)
{
  data.index = 0;

  // If the area vector does not contain four values, then the box
  // stay empty. NO INFLOW.
  if(ie.area.size()==4)
  {
    // Compute the given box area.
    data.box_area = CA::Box::create(_grid, ie.area[0], ie.area[1], ie.area[2], ie.area[3]);
  }
  
  // This is new. The area could have been identified as a zone, i.e. x,y,w,h/
  if(ie.zone.size()==4)
  {
    // Compute the given box area from the zone.
    CA::Point     tl( CA::Point::create(_grid, ie.zone[0], ie.zone[1]) );
    CA::Unsigned  w = static_cast<CA::Unsigned>( std::ceil(ie.zone[2] / _grid.length()) );
    CA::Unsigned  h = static_cast<CA::Unsigned>( std::ceil(ie.zone[3] / _grid.length()) );
    
    data.box_area = CA::Box(tl.x(),tl.y(),w,h);
  }


  data.grid_area = 0.0;
  data.volume    = 0.0;


  return 0;
}
