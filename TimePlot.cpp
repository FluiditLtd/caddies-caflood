//! \file TimePlot.cpp
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2013-03


#include"ca2D.hpp"
#include"ArgsData.hpp"
#include"TimePlot.hpp"
#include"Utilities.hpp"
#include<iostream>
#include<fstream>


// Initialise the TimePlot structure usign a CSV file. 
int initTimePlotFromCSV(const std::string& filename, TimePlot& tp)
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

    if(CA::compareCaseInsensitive("Time Plot Name",tokens[0],true))
    {
      std::string str;
      READ_TOKEN(str,tokens[1],tokens[0]);
      
      tp.name = CA::trimToken(str);
    }

    if(CA::compareCaseInsensitive("Physical Variable",tokens[0],true))
      READ_TOKEN(tp.pv,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Points Name",tokens[0],true))
    {
      for (size_t i=1; i<tokens.size(); ++i)
      {
	std::string str;
	READ_TOKEN(str,tokens[i],tokens[0]);

	tp.pnames.push_back(CA::trimToken(str));
      }
    }

    if(CA::compareCaseInsensitive("Points X Coo",tokens[0],true))
    {
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(value,tokens[i],tokens[0]);

	tp.xcoos.push_back(value);
      }
    }

    if(CA::compareCaseInsensitive("Points Y Coo",tokens[0],true))
    {
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Real value;
	READ_TOKEN(value,tokens[i],tokens[0]);

	tp.ycoos.push_back(value);
      }
    }
    
    if(CA::compareCaseInsensitive("Period",tokens[0],true))
      READ_TOKEN(tp.period,tokens[1],tokens[0]);

  }

  return 0;
}


int initTPData(const std::string& filename, const CA::Grid&  GRID, CA::CellBuffReal&  ELV,
	       const TimePlot& tp, TPData& tpdata)
{
  // Create file
  tpdata.filename = filename;
  tpdata.file.reset( new std::ofstream(filename.c_str()) );

  if(!tpdata.file->good())
    return 1;

  // Set  manipulators 
  tpdata.file->setf(std::ios::fixed, std::ios::floatfield);
  tpdata.file->precision(6); 

  // Write the header
  (*tpdata.file)<<"Iter, Time (s), ";

  // Write point name.
  switch(tp.pv)
  {
  case PV::WD:
  case PV::WL:
  case PV::VEL:    
    // Only once for water depth and water level.
    for(size_t p =0; p< tp.pnames.size(); p++)
    {	
      (*tpdata.file)<<tp.pnames[p]<<", ";
    }
    break;
  }
  (*tpdata.file)<<std::endl;

  
  // Loop through coordinates,
  for(size_t p =0; p< tp.pnames.size(); p++)
  {	
    tpdata.pl.add( CA::Point::create(GRID,tp.xcoos[p],tp.ycoos[p]) );
  }

  // Create buffer where to store the point data.
  tpdata.pvals.resize(tpdata.pl.size());

  // If we need to save the weater level.
  if( tp.pv == PV::WL )
  {
    // Create buffer where to store the elevation data.
    tpdata.pelvs.resize(tpdata.pl.size());
    // Retrieve the levation data.
    ELV.retrievePoints(tpdata.pl,&(tpdata.pelvs[0]),tpdata.pl.size());      
  }

  if(tp.period > 0.0)
    tpdata.time_next = tp.period;
  else
    tpdata.time_next = std::numeric_limits<CA::Real>::max();

  return 0;
}
