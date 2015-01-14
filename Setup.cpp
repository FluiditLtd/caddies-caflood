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
  setup.sim_name           ="Sim";
  setup.short_name         ="sim";
  setup.preproc_name       = "";
  setup.model_type         = MODEL::WCA2Dv2; // default.
 
  setup.roughness_global   = 0.01;
  setup.infrate_global     = 0.0;

  setup.tolerance          = 0.001;
  setup.tol_slope          = 0.1;
  setup.ignore_wd          = 0.001;

  setup.output_console     = false ;
  setup.terrain_info       = false ;
  setup.ts_plot            = false ;
  setup.output_period      = 300   ;
  setup.output_computation = false ;
  setup.check_vols         = false ;
  setup.remove_data        = true ;
  setup.remove_prec_data   = true ;
  setup.rast_vel_as_vect   = true ;
  setup.rast_wd_tol        = 0.01 ;
  setup.rast_boundary      = false;
  setup.rast_places        = 6;
  setup.update_peak_dt     = false;
  setup.expand_domain      = false;
  setup.ignore_upstream    = false;
  setup.upstream_reduction = 1.0;

  
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
    // If true the token was identified;
    bool found_tok = false;

    std::vector<std::string> tokens( CA::getLineTokens(ifile, ',') );
    
    // If the tokens vector is empty we reached the eof or an
    // empty line... continue.
    if(tokens.empty())
      continue;

    if(CA::compareCaseInsensitive("Simulation Name",tokens[0],true))
    {
      std::string str;
      READ_TOKEN(found_tok,str,tokens[1],tokens[0]);

      setup.sim_name = CA::trimToken(str);
    }

    if(CA::compareCaseInsensitive("Short Name",tokens[0],true))
    {
      std::string str;
      READ_TOKEN(found_tok,str,tokens[1],tokens[0]);

      setup.short_name = CA::trimToken(str);
    }

    if(CA::compareCaseInsensitive("Pre-proc Name",tokens[0],true))
    {
      std::string str;
      READ_TOKEN(found_tok,str,tokens[1],tokens[0]);

      setup.preproc_name = CA::trimToken(str);
    }

    if(CA::compareCaseInsensitive("Version",tokens[0],true))
    {
      found_tok=true;
      for (size_t i=1; i<tokens.size(); ++i)
      {
	CA::Unsigned value;
	READ_TOKEN(found_tok,value,tokens[i],tokens[0]);

	setup.sim_version.push_back(value);
      }
    }

    if(CA::compareCaseInsensitive("Model Type",tokens[0],true))
      READ_TOKEN(found_tok,setup.model_type,tokens[1],tokens[0]);


    if(CA::compareCaseInsensitive("Time Start",tokens[0],true))
      READ_TOKEN(found_tok,setup.time_start,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Time End",tokens[0],true))
      READ_TOKEN(found_tok,setup.time_end,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Max DT",tokens[0],true))
      READ_TOKEN(found_tok,setup.time_maxdt,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Min DT",tokens[0],true))
      READ_TOKEN(found_tok,setup.time_mindt,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Update DT",tokens[0],true))
      READ_TOKEN(found_tok,setup.time_updatedt,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Alpha",tokens[0],true))
      READ_TOKEN(found_tok,setup.time_alpha,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Max Iter",tokens[0],true))
      READ_TOKEN(found_tok,setup.time_maxiters,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Roughness Global",tokens[0],true))
      READ_TOKEN(found_tok,setup.roughness_global,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Infiltration Global",tokens[0],true))
      READ_TOKEN(found_tok,setup.infrate_global,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Tolerance",tokens[0],true))
      READ_TOKEN(found_tok,setup.tolerance,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Slope Tol",tokens[0],true))
    {
      // Attention the tollerance must be passed in percentile.
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(found_tok,setup.tol_slope,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Ignore WD",tokens[0],true))
      READ_TOKEN(found_tok,setup.ignore_wd,tokens[1],tokens[0]);


    if(CA::compareCaseInsensitive("Boundary Ele",tokens[0],true))
      READ_TOKEN(found_tok,setup.boundary_elv,tokens[1],tokens[0]);


    if(CA::compareCaseInsensitive("Elevation ASCII",tokens[0],true))
    {
      std::string str;
      READ_TOKEN(found_tok,str,tokens[1],tokens[0]);
      
      setup.elevation_ASCII = CA::trimToken(str);
    }

    if(CA::compareCaseInsensitive("Rain Event CSV",tokens[0],true))
    {
      found_tok=true;
      for (size_t i=1; i<tokens.size(); ++i)
      {
	std::string str;
	READ_TOKEN(found_tok,str,tokens[i],tokens[0]);

	setup.rainevent_files.push_back(CA::trimToken(str));
      }
    }

    if(CA::compareCaseInsensitive("Water Level Event CSV",tokens[0],true))
    {
      found_tok=true;
      for (size_t i=1; i<tokens.size(); ++i)
      {
	std::string str;
	READ_TOKEN(found_tok,str,tokens[i],tokens[0]);

	setup.wlevent_files.push_back(CA::trimToken(str));
      }
    }

    if(CA::compareCaseInsensitive("Inflow Event CSV",tokens[0],true))
    {
      found_tok=true;
      for (size_t i=1; i<tokens.size(); ++i)
      {
	std::string str;
	READ_TOKEN(found_tok,str,tokens[i],tokens[0]);

	setup.inflowevent_files.push_back(CA::trimToken(str));
      }
    }

    if(CA::compareCaseInsensitive("Time Plot CSV",tokens[0],true))
    {
      found_tok=true;
      for (size_t i=1; i<tokens.size(); ++i)
      {
	std::string str;
	READ_TOKEN(found_tok,str,tokens[i],tokens[0]);

	setup.timeplot_files.push_back(CA::trimToken(str));
      }
    }

    if(CA::compareCaseInsensitive("Raster Grid CSV",tokens[0],true))
    {
      found_tok=true;
      for (size_t i=1; i<tokens.size(); ++i)
      {
	std::string str;
	READ_TOKEN(found_tok,str,tokens[i],tokens[0]);

	setup.rastergrid_files.push_back(CA::trimToken(str));
      }
    }

    if(CA::compareCaseInsensitive("Output Console",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(found_tok,setup.output_console,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Output Period",tokens[0],true))
      READ_TOKEN(found_tok,setup.output_period,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Terrain Info",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(found_tok,setup.terrain_info,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("TS Plot",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(found_tok,setup.ts_plot,str,tokens[0]);
    }
    
    if(CA::compareCaseInsensitive("Output Computation Time",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(found_tok,setup.output_computation,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Check Volumes",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(found_tok,setup.check_vols,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Remove Proc Data",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(found_tok,setup.remove_data,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Remove Pre-Proc Data",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(found_tok,setup.remove_prec_data,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Raster VEL Vector",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(found_tok,setup.rast_vel_as_vect,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Raster WD Tol",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(found_tok,setup.rast_wd_tol,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Raster Boundary",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(found_tok,setup.rast_boundary,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Raster Decimal Places",tokens[0],true))
      READ_TOKEN(found_tok,setup.rast_places,tokens[1],tokens[0]);

    if(CA::compareCaseInsensitive("Update Peak Every DT",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(found_tok,setup.update_peak_dt,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Expand Domain",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(found_tok,setup.expand_domain,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Ignore Upstream",tokens[0],true))
    {
      std::string str = CA::trimToken(tokens[1]);
      READ_TOKEN(found_tok,setup.ignore_upstream,str,tokens[0]);
    }

    if(CA::compareCaseInsensitive("Upstream Reduction",tokens[0],true))
      READ_TOKEN(found_tok,setup.upstream_reduction,tokens[1],tokens[0]);

    // If the token was not identified stop!
    if(!found_tok)
    {
      std::cerr<<"Element '"<<CA::trimToken(tokens[0])<<"' not identified"<<std::endl; \
      return 1;
    }
  }

  // The updatedt cannot be more than 60 seconds.
  setup.time_updatedt = std::min(static_cast<CA::Real>(60.0),setup.time_updatedt);

  // The maximum dt cannot be more than update dt.
  setup.time_maxdt = std::min(setup.time_maxdt,setup.time_updatedt);

  // The end time must be a multiple of updatedt.
  // Get the reminder
  CA::Real r = std::fmod(setup.time_end, setup.time_updatedt);
  if(r>0.0)
    setup.time_end += setup.time_updatedt - r;
  
  // If the preproc base name is empty use the simulation short name.
  if(setup.preproc_name.empty())
    setup.preproc_name = setup.short_name;

  return 0;
}
