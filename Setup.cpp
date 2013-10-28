//! \file Setup.cpp
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2013-03


#include"ArgsData.hpp"
#include"Setup.hpp"
#include"Utilities.hpp"
#include<iostream>
#include<fstream>


// Initialise the setup structure usign a CSV file. 
int initSetupFromCSV(const std::string& filename, Setup& setup)
{
  // Default values
  setup.tolerance          = 0.001;
  setup.ignore_wd          = 0.001;

  setup.output_console     = false ;
  setup.output_period      = 300   ;
  setup.output_computation = false ;
  setup.check_vols         = false ;
  setup.remove_data        = true ;
  setup.remove_prec_data   = true ;
  setup.rast_vel_as_vect   = true ;
  setup.rast_wd_tol        = 0.01 ;
  setup.update_peak_dt     = false;
  setup.expand_domain      = false;
  
  // Read values
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

    if(CA::compareCaseInsensitive("Simulation Name",tokens[0],true))
    {
      std::string str;
      READ_TOKEN(str,tokens[1],tokens[0]);

      setup.sim_name = CA::trimToken(str);
    }

    if(CA::compareCaseInsensitive("Short Name",tokens[0],true))
    {
      std::string str;
      READ_TOKEN(str,tokens[1],tokens[0]);

      setup.short_name = CA::trimToken(str);
    }

    if(CA::compareCaseInsensitive("Version",tokens[0],true))
    {
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Unsigned value;
	READ_TOKEN(value,tokens[i],tokens[0]);

	setup.sim_version.push_back(value);
      }
    }

    if(CA::compareCaseInsensitive("Time Start",tokens[0],true))
      READ_TOKEN(setup.time_start,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Time End",tokens[0],true))
      READ_TOKEN(setup.time_end,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Max DT",tokens[0],true))
      READ_TOKEN(setup.time_maxdt,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Min DT",tokens[0],true))
      READ_TOKEN(setup.time_mindt,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Update DT",tokens[0],true))
      READ_TOKEN(setup.time_updatedt,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Alpha",tokens[0],true))
      READ_TOKEN(setup.time_alpha,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Max Iter",tokens[0],true))
      READ_TOKEN(setup.time_maxiters,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Roughness Global",tokens[0],true))
      READ_TOKEN(setup.roughness_global,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Tolerance",tokens[0],true))
      READ_TOKEN(setup.tolerance,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Ignore WD",tokens[0],true))
      READ_TOKEN(setup.ignore_wd,tokens[1],tokens[0]);


    if(CA::compareCaseInsensitive("Boundary Ele",tokens[0],true))
      READ_TOKEN(setup.boundary_elv,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Elevation ASCII",tokens[0],true))
    {
      std::string str;
      READ_TOKEN(str,tokens[1],tokens[0]);
      
      setup.elevation_ASCII = CA::trimToken(str);
    }

    if(CA::compareCaseInsensitive("Rain Event CSV",tokens[0],true))
    {
      for (size_t i=1; i<tokens.size(); ++i)
      {
	std::string str;
	READ_TOKEN(str,tokens[i],tokens[0]);

	setup.rainevent_files.push_back(CA::trimToken(str));
      }
    }

    if(CA::compareCaseInsensitive("Water Level Event CSV",tokens[0],true))
    {
      for (size_t i=1; i<tokens.size(); ++i)
      {
	std::string str;
	READ_TOKEN(str,tokens[i],tokens[0]);

	setup.wlevent_files.push_back(CA::trimToken(str));
      }
    }

    if(CA::compareCaseInsensitive("Inflow Event CSV",tokens[0],true))
    {
      for (size_t i=1; i<tokens.size(); ++i)
      {
	std::string str;
	READ_TOKEN(str,tokens[i],tokens[0]);

	setup.inflowevent_files.push_back(CA::trimToken(str));
      }
    }

    if(CA::compareCaseInsensitive("Time Plot CSV",tokens[0],true))
    {
      for (size_t i=1; i<tokens.size(); ++i)
      {
	std::string str;
	READ_TOKEN(str,tokens[i],tokens[0]);

	setup.timeplot_files.push_back(CA::trimToken(str));
      }
    }

    if(CA::compareCaseInsensitive("Raster Grid CSV",tokens[0],true))
    {
      for (size_t i=1; i<tokens.size(); ++i)
      {
	std::string str;
	READ_TOKEN(str,tokens[i],tokens[0]);

	setup.rastergrid_files.push_back(CA::trimToken(str));
      }
    }

    if(CA::compareCaseInsensitive("Output Console",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(setup.output_console,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Output Period",tokens[0],true))
      READ_TOKEN(setup.output_period,tokens[1],tokens[0]);
    
    if(CA::compareCaseInsensitive("Output Computation Time",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(setup.output_computation,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Check Volumes",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(setup.check_vols,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Remove Proc Data",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(setup.remove_data,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Remove Pre-Proc Data",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(setup.remove_prec_data,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Raster VEL Vector",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(setup.rast_vel_as_vect,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Raster WD Tol",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(setup.rast_wd_tol,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Update Peak Every DT",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(setup.update_peak_dt,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Expand Domain",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(setup.expand_domain,str,tokens[0]);
    }

  }

  // The updatedt cannot be more than 60 seconds.
  setup.time_updatedt = std::min(static_cast<CA::Real>(60.0),setup.time_updatedt);

  // The maximum dt cannot be more than update dt.
  setup.time_maxdt = std::min(setup.time_maxdt,setup.time_updatedt);

  return 0;
}
