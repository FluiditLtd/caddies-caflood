//! \file Events.cpp
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2013-03


#include"ArgsData.hpp"
#include"Events.hpp"
#include"Utilities.hpp"
#include<iostream>
#include<fstream>


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
    std::vector<std::string> tokens( CA::getLineTokens(ifile, ',') );
    
    // If the tokens vector is empty we reached the eof or an
    // empty line... continue.
    if(tokens.empty())
      continue;       

    if(CA::compareCaseInsensitive("Event Name",tokens[0],true))
    {
      std::string str;
      READ_TOKEN(str,tokens[1],tokens[0]);
      
      wle.name = CA::trimToken(str," \t\r");
    }

    if(CA::compareCaseInsensitive("Water Level",tokens[0],true))
    {
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(value,tokens[i],tokens[0]);

	wle.wls.push_back(value);
      }
    }

    if(CA::compareCaseInsensitive("Time",tokens[0],true))
    {
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(value,tokens[i],tokens[0]);

	wle.times.push_back(value);
      }
    }

    if(CA::compareCaseInsensitive("Area",tokens[0],true))
    {
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(value,tokens[i],tokens[0]);

	wle.area.push_back(value);
      }
    }

    if(CA::compareCaseInsensitive("Zone",tokens[0],true))
    {
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(value,tokens[i],tokens[0]);

	wle.zone.push_back(value);
      }
    }
  }
  
  return 0;
}

int initWLEData(const CA::Grid&  GRID, const WLEvent& wle, WLEData& wledata)
{
  wledata.index = 0;

  // If the area vector does not contain four values, then the box
  // stay empty. NO WATER LEVEL.
  if(wle.area.size()==4)
  {
    // Compute the given box area.
    wledata.box_area = CA::Box::create(GRID, wle.area[0], wle.area[1], wle.area[2], wle.area[3]);
  }
  // This is new. The area could have been identified as a zone, i.e. x,y,w,h/
  if(wle.zone.size()==4)
  {
    // Compute the given box area from the zone.
    CA::Point     tl( CA::Point::create(GRID, wle.zone[0], wle.zone[1]) );
    CA::Unsigned  w = static_cast<CA::Unsigned>( std::ceil(wle.zone[2] / GRID.length()) );
    CA::Unsigned  h = static_cast<CA::Unsigned>( std::ceil(wle.zone[3] / GRID.length()) );
    
    wledata.box_area = CA::Box(tl.x(),tl.y(),w,h);
  }

  return 0;
}


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

int initIEData(const CA::Grid&  GRID, const IEvent& ie, IEData& iedata)
{
  iedata.index = 0;

  // If the area vector does not contain four values, then the box
  // stay empty. NO INFLOW.
  if(ie.area.size()==4)
  {
    // Compute the given box area.
    iedata.box_area = CA::Box::create(GRID, ie.area[0], ie.area[1], ie.area[2], ie.area[3]);
  }
  
  // This is new. The area could have been identified as a zone, i.e. x,y,w,h/
  if(ie.zone.size()==4)
  {
    // Compute the given box area from the zone.
    CA::Point     tl( CA::Point::create(GRID, ie.zone[0], ie.zone[1]) );
    CA::Unsigned  w = static_cast<CA::Unsigned>( std::ceil(ie.zone[2] / GRID.length()) );
    CA::Unsigned  h = static_cast<CA::Unsigned>( std::ceil(ie.zone[3] / GRID.length()) );
    
    iedata.box_area = CA::Box(tl.x(),tl.y(),w,h);
  }


  iedata.grid_area = 0.0;
  iedata.volume    = 0.0;

  return 0;
}
