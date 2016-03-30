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

//! \file postProc.cpp
//! Perform the post processing of the data for a CA 2D model.
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2013-07


#include"ca2D.hpp"
#include"Masks.hpp"
#include"ArgsData.hpp"
#include"Setup.hpp"
#include"TimePlot.hpp"
#include"RasterGrid.hpp"


// -------------------------//
// Include the CA 2D functions //
// -------------------------//
#include CA_2D_INCLUDE(zeroedWD)
#include CA_2D_INCLUDE(zeroedVA)


//! Inline function that given the actual time, the time of the
//! nearest action and the period of the possible new action, find the
//! possible NEW simulation time of possible new nearest action.
//! \attention An action is: the simulation finished, output to
//! console, raster, timeplot, or rain event change, etc..
inline CA::Real nextTimeNearestAction(CA::Real t, CA::Real t_nearest, CA::Real period)
{
  return std::min(t_nearest, t + period - std::fmod(t,period));
}


// -------------------------//
// Include the CA 2D functions //
// -------------------------//
#include CA_2D_INCLUDE(makeWL)



//! Perform the post processing of the data for a CA 2D model. 
int postProc(const ArgsData& ad, const Setup& setup, CA::AsciiGrid<CA::Real>& eg,
	     const std::vector<TimePlot>& tps, const std::vector<RasterGrid>& rgs)
{
  
  if(setup.output_computation)
  {
    std::cout<<"Post-processing : "<<setup.sim_name<< std::endl;
    std::cout<<"------------------------------------------" << std::endl; 
  }

  // ----  Timer ----
  
  // Get starting time.
  CA::Clock total_timer;

  // ----  CA GRID ----
  
  // Load the CA Grid from the DataDir. 
  // ATTENTION this should have an extra set of cells in each
  // direction.  The internal implementation could be different than a
  // square regular grid.
  CA::Grid  GRID(ad.data_dir,setup.preproc_name+"_Grid","0", ad.args.active());

  if(setup.output_console)
    std::cout<<"Loaded Grid data"<< std::endl;

  // Set if to print debug information on CA function.
  GRID.setCAPrint(false);

  // Create the full (extended) computational domain of CA grid. 
  CA::BoxList  fulldomain;
  CA::Box      fullbox = GRID.box();
  fulldomain.add(fullbox);

  // Create a borders object that contains all the borders and
  // corners of the grid.
  CA::Borders borders;
  
  borders.addSegment(CA::Top);
  borders.addSegment(CA::Bottom);
  borders.addSegment(CA::Right);
  borders.addSegment(CA::Left);
  
  borders.addCorner(CA::TopLeft);
  borders.addCorner(CA::TopRight);
  borders.addCorner(CA::BottomLeft);
  borders.addCorner(CA::BottomRight);

  // Create the real computational domain of CA grid, i.e. the
  // original DEM size not the extended one.
  CA::BoxList  realdomain;
  CA::Box      realbox(GRID.box().x()+1,GRID.box().y()+1,GRID.box().w()-2,GRID.box().h()-2);
  realdomain.add(realbox);

  // CREATE ASCII GRID TEMPORARY BUFFERS
  CA::AsciiGrid<CA::Real>& agtmp1 = eg;
  CA::AsciiGrid<CA::Real>  agtmp2 = eg;

  // -- INITIALISE ELEVATION ---

  // Create the elevation cell buffer.
  // It contains a "real" value in each cell of the grid.
  CA::CellBuffReal  ELV(GRID);

  // Create the water depth buffer.
  // It contains a "real" value in each cell of the grid.
  CA::CellBuffReal  WD(GRID);

  // Set the border of the elevation buffer to be no data. 
  ELV.bordersValue(borders,agtmp1.nodata);

  // Se the default value of the elevation to be nodata.
  ELV.fill(fulldomain, agtmp1.nodata);

  // Load the data not from the DEM file but from the pre-processed
  // file.
  if(!ELV.loadData(setup.preproc_name+"_ELV","0") )
  {
    std::cerr<<"Error while loading the Elevation pre-processed file"<<std::endl;
    return 1;
  }

  if(setup.output_console)
    std::cout<<"Loaded Elevation data"<< std::endl;

  // ----  CELL BUFFERS ----
    
  // Create two temporary Cell buffer
  CA::CellBuffReal TMP1(GRID);
  CA::CellBuffReal TMP2(GRID);

  // Create the MASK cell buffer. The mask is usefull to check wich
  // cell has data and nodata and which cell has neighbourhood with
  // data.
  CA::CellBuffState MASK(GRID);  
  
  // ----  SCALAR VALUES ----

  CA::Unsigned iter      = 0;		        // The actual iteration number.
  CA::Real     t         = setup.time_start;    // The actual time in seconds
  CA::Real     t_nearest = setup.time_end;      // The next time step
  CA::Real     nodata    = agtmp1.nodata;

  // -- CREATE FULL MASK ---
  
  CA::createCellMask(fulldomain,GRID,ELV,MASK,nodata);  

  // ----  INIT RASTER GRID ----

  // List of raster grid data
  std::vector<RGData> rgdatas(rgs.size());

  // Peak buffers.
  RGPeak rgpeak;

  // Allocate temporary Ascii file data to ouptut raster.
  if(rgs.size()>0)
  {
    agtmp1.data.resize(agtmp1.ncols * agtmp1.nrows, agtmp1.nodata);
    agtmp2.data.resize(agtmp2.ncols * agtmp2.nrows, agtmp2.nodata);
  }

  for(size_t i = 0; i<rgs.size(); ++i)
  {
    std::string filename = ad.output_dir+ad.sdir+setup.short_name+"_"+setup.rastergrid_files[i]; 
    initRGData(filename, GRID, nodata, rgs[i], rgdatas[i], rgpeak);
  }

  // ---- CONTAINER DATA TO REMOVE ---

  // Create the container with the data to remove.
  typedef std::vector< std::pair<std::string, std::string> > RemoveID;
  RemoveID removeIDsCB;
  RemoveID removeIDsEB;


  // Indicate that the water depth buffer was loaded.
  bool WDloaded  = false;

  // ------------------------- MAIN LOOP -------------------------------
  while(t<=setup.time_end)
  {
    // Need to realod WD for this time step.
    WDloaded  = false;

    // -------  OUTPUTS --------    

    // Retrieve the string of the time.
    std::string strtime;
    CA::toString(strtime,t);
    
    // ----  RASTER GRID ----

    for(size_t i = 0; i<rgdatas.size(); ++i)
    {
      // Check if it is time to process raster grid!
      // Or the final extend need to be processed
      if(t >= rgdatas[i].time_next || (rgs[i].final && t==setup.time_end) )
      {
	if(!WDloaded)
	{
	  // The water depth buffer is practically always needed (for WD/WL and VEL).
	  // Reset the buffer
	  WD.fill(fulldomain,agtmp1.nodata);
	  
	  // Load the water depth data.
	  if(!WD.loadData(setup.short_name+"_WD",strtime) )
	  {
	    std::cerr<<"Missing water depth data: "<<strtime<<std::endl;
	    continue;
	  }
	  
	  // Set the water depth to zero if it less than tollerance.
	  CA::State boundary = (setup.rast_boundary) ? 1 : 0;
	  CA::Execute::function(fulldomain, zeroedWD, GRID, WD, MASK, setup.rast_wd_tol, boundary);
	
	  // Add the ID to remove.
	  removeIDsCB.push_back(std::make_pair(setup.short_name+"_WD",strtime));

	  WDloaded=true;
	}
	// Perform the action depending on the type of variable.
	switch(rgs[i].pv)
	{
	case PV::VEL:
	  {	  
	    std::string filenameV;
	    std::string filenameA;
	    
	    // Create the name of the files if it is going to output a
	    // velocity field or simply rasters.
	    if(setup.rast_vel_as_vect == true)
	    {
	      filenameV = removeExtension(rgdatas[i].filename) + "_" + strtime + ".csv";
	      filenameA = removeExtension(rgdatas[i].filename) + "_" + strtime + ".csvt";
	    }
	    else
	    {
	      filenameV = removeExtension(rgdatas[i].filename) + "_V_" + strtime;
	      filenameA = removeExtension(rgdatas[i].filename) + "_A_" + strtime;
	    }
	    
	    // Reset the buffer
	    TMP1.fill(fulldomain,agtmp1.nodata);
	    TMP2.fill(fulldomain,agtmp1.nodata);
	    
	    // Load the velocity data on TMP1.
	    if(! TMP1.loadData(setup.short_name+"_V",strtime) )
	    {
	      std::cerr<<"Missing velocity data to create file: "<<filenameV<<std::endl;
	      continue;
	    }
	    
	    // Load the angle data on TMP2.
	    if(! TMP2.loadData(setup.short_name+"_A",strtime) )
	    {
	      std::cerr<<"Missing angle data to create file: "<<filenameA<<std::endl;
	      continue;
	    }

	    // Set the V and A to zero if water depth is less than tollerance.
	    CA::Execute::function(fulldomain, zeroedVA, GRID, TMP1, TMP2, WD, MASK, setup.rast_wd_tol);
	    
	    // Retrieve the velocity and angle data
	    TMP1.retrieveData(realbox, &agtmp1.data[0], agtmp1.ncols, agtmp1.nrows);	  
	    TMP2.retrieveData(realbox, &agtmp2.data[0], agtmp2.ncols, agtmp2.nrows);	  
	    
	    // check if we need to output a velocity field
	    if(setup.rast_vel_as_vect)
	    {
	      if(setup.output_console)
		std::cout<<"Write Raster Grid: "<<filenameV<<std::endl;
	      
	      // Create the CSV file
	      //std::ofstream file( filenameV.c_str() ); 
	      FILE* fout = fopen(filenameV.c_str(), "w");	      

	      // Set  manipulators 
	      //file.setf(std::ios::fixed, std::ios::floatfield);
	      //file.precision(6); 
	      
	      // Write the header
	      //file<<"X, Y, Speed, Angle_RAD, Angle_DEG, Angle_QGIS"<<std::endl;
	      fprintf(fout,"X, Y, Speed, Angle_RAD, Angle_DEG, Angle_QGIS\n");

	      // Loop thourhg the grid points.
	      for(CA::Unsigned j_reg=realbox.y(), j_mem=0; j_reg<realbox.h()+realbox.y(); ++j_reg, ++j_mem)
	      {
		for(CA::Unsigned i_reg=realbox.x(), i_mem=0; i_reg<realbox.w()+realbox.x(); ++i_reg, ++i_mem)
		{
		  // Create the point and find the coordinates.
		  CA::Point p(i_reg,j_reg);
		  p.setCoo(GRID);
		  
		  // Retrieve the speed and angle values (in radians),
		  // the compute the angle value in degrees.
		  CA::Real V  = agtmp1.data[j_mem * agtmp1.ncols + i_mem];
		  CA::Real AR = agtmp2.data[j_mem * agtmp2.ncols + i_mem];
		  CA::Real AD = AR*180/PI;
		  
		  // Compute the angle needed by QGis.
		  CA::Real AQ = -AR*180/PI + 90; 		
		  
		  // Write the results if the veclocity is more than zero.
		  if(V>0)
		  {
		    //file<<p.coo().x()<<","<<p.coo().y()<<","<<V<<","<<AR<<","<<AD<<","<<AQ<<","<<std::endl;
		    fprintf(fout,"%.12f,%.12f,%.6f,%.6f,%.6f,%.6f,\n", p.coo().x(),p.coo().y(),V,AR,AD,AQ);
		  }
		}
	      }        
	      
	      // Close the file.
	      //file.close();	      
	      fclose(fout);
	    }
	    // NOPE, simply output the rasters.
	    else
	    {
	      if(setup.output_console)
		std::cout<<"Write Raster Grid: "<<filenameV<<" "<<filenameA<<std::endl;
	      
	      // Write the  data.
		  agtmp1.writeAsciiGrid(filenameV,setup.rast_places);
		  agtmp2.writeAsciiGrid(filenameA,setup.rast_places);
	    }
	    
	    // Add the ID to remove.
	    removeIDsCB.push_back(std::make_pair(setup.short_name+"_V",strtime));
	    removeIDsCB.push_back(std::make_pair(setup.short_name+"_A",strtime));
	  }
	  break;
	  
	case PV::WL:
	  {
	    // Create the name of the file.
	    std::string filename = removeExtension(rgdatas[i].filename) + "_" + strtime;

	    // Reset the buffer
	    TMP1.fill(fulldomain,agtmp1.nodata);

	    // Make the water depth data and elevation into the water level.
	    CA::Execute::function(fulldomain, makeWL, GRID, TMP1, WD, ELV, MASK);
	    
	    // Retrieve the data
	    TMP1.retrieveData(realbox, &agtmp1.data[0], agtmp1.ncols, agtmp1.nrows);	  
	    	    
	    if(setup.output_console)
	      std::cout<<"Write Raster Grid: "<<filename<<std::endl;
	    
	    // Write the data.
		agtmp1.writeAsciiGrid(filename,setup.rast_places);
	  }
	  break;
	case PV::WD:
	  {
	    // Create the name of the file.
	    std::string filename = removeExtension(rgdatas[i].filename) + "_" + strtime;

	    // Water depth already loaded.
	    
	    // Retrieve the data
	    WD.retrieveData(realbox, &agtmp1.data[0], agtmp1.ncols, agtmp1.nrows);	  
	    	    
	    if(setup.output_console)
	      std::cout<<"Write Raster Grid: "<<filename<<std::endl;
	    
	    // Write the data.
		agtmp1.writeAsciiGrid(filename,setup.rast_places);
	  }
	  break;
	default:
	  break;
	}

	// Update the next time to process a raster grid.
	rgdatas[i].time_next += rgs[i].period;
      }

      // Its time to find the possible next time of interest.
      t_nearest = nextTimeNearestAction(t,t_nearest, rgs[i].period);

    }

    //Finish after t reach end
    if(t == setup.time_end)
      break;

    // Set the nearest important time as the enxt time and set the
    // possible nearest important time as the end of the simulation.
    t         = t_nearest;
    t_nearest = setup.time_end;   
  }

  // ----  PEAK RASTER GRID ----

  // Need to realod WD only once for PEAK.
  WDloaded  = false;
  
  for(size_t i = 0; i<rgdatas.size(); ++i)
  {
    // Check if the peak values need to be saved.
    if(rgs[i].peak == true)
    {

      if(!WDloaded)
      {
	// The water depth buffer is practically always needed (for WD/WL and VEL).
	// Reset the buffer
	WD.fill(fulldomain,agtmp1.nodata);
	
	// Load the water depth data.
	if(! WD.loadData(setup.short_name+"_WD","PEAK") )
	{
	  std::cerr<<"Missing water depth data: "<<"PEAK"<<std::endl;
	  continue;
	}

	// Set the water depth to zero if it less than tollerance.
	CA::State boundary = (setup.rast_boundary) ? 1 : 0;
	CA::Execute::function(fulldomain, zeroedWD, GRID, WD, MASK, setup.rast_wd_tol, boundary);

	// Add the ID to remove.
	removeIDsCB.push_back(std::make_pair(setup.short_name+"_WD","PEAK"));

	// Only once.
	WDloaded = true;
      }
	
      switch(rgs[i].pv)
      {
      case PV::VEL:
	{
	  // Create the name of the file.
	  std::string filenameV = removeExtension(rgdatas[i].filename) + "_V_PEAK";
	  
	  // Reset the buffer
	  TMP1.fill(fulldomain,agtmp1.nodata);
	  
	  // Load the data on TMP1.
	  if(! TMP1.loadData(setup.short_name+"_V","PEAK") )
	  {
	    std::cerr<<"Missing velocity data to create file: "<<filenameV<<std::endl;
	    continue;
	  }

	  // Set the V and A to zero if water depth is less than tollerance.
	  CA::Execute::function(fulldomain, zeroedVA, GRID, TMP1, TMP2, WD, MASK, setup.rast_wd_tol);
	  	  
	  // Retrieve the data
	  TMP1.retrieveData(realbox, &agtmp1.data[0], agtmp1.ncols, agtmp1.nrows);	  
	  
	  if(setup.output_console)
	    std::cout<<"Write Raster Grid: "<<filenameV<<std::endl;
	  
	  // Write the data.
	  agtmp1.writeAsciiGrid(filenameV,setup.rast_places);

	  // Add the ID to remove.
	  removeIDsCB.push_back(std::make_pair(setup.short_name+"_V","PEAK"));
	}
	break;

      case PV::WL:
	{
	  // Create the name of the file.
	  std::string filename = removeExtension(rgdatas[i].filename) + "_PEAK";
	  
	  // Reset the buffer
	  TMP1.fill(fulldomain,agtmp1.nodata);
	  	  

	  // Make the water depth data and elevation into the water level.
	  CA::Execute::function(fulldomain, makeWL, GRID, TMP1, WD, ELV, MASK);
	  
	  // Retrieve the data
	  TMP1.retrieveData(realbox, &agtmp1.data[0], agtmp1.ncols, agtmp1.nrows);	  
	  
	  if(setup.output_console)
	    std::cout<<"Write Raster Grid: "<<filename<<std::endl;
	  
	  // Write the data.
	  agtmp1.writeAsciiGrid(filename,setup.rast_places);
	}
	break;
	
      case PV::WD:
	{
	  // Create the name of the file.
	  std::string filename = removeExtension(rgdatas[i].filename) + "_PEAK";

	  // Water depth already loaded.
	  
	  // Retrieve the data
	  WD.retrieveData(realbox, &agtmp1.data[0], agtmp1.ncols, agtmp1.nrows);	  
	  
	  if(setup.output_console)
	    std::cout<<"Write Raster Grid: "<<filename<<std::endl;
	  
	  // Write the data.
	  agtmp1.writeAsciiGrid(filename,setup.rast_places);
	}
	break;
      default:
	break;
      }      
    }

  }

  if(setup.remove_data)
  {
    for(size_t i=0; i<removeIDsCB.size(); i++)
    {
      std::pair <std::string,std::string>& ID = removeIDsCB[i];; 
      CA::CellBuffReal::removeData(ad.data_dir,ID.first,ID.second);
    }    
    for(size_t i=0; i<removeIDsEB.size(); i++)
    {
      std::pair <std::string,std::string>& ID = removeIDsEB[i];; 
      CA::EdgeBuffReal::removeData(ad.data_dir,ID.first,ID.second);
    }    
  }

  if(setup.remove_prec_data)
  {
    // Remove Elevation data.
    CA::CellBuffReal::removeData(ad.data_dir,setup.preproc_name+"_ELV","0");
    
    // Remove Grid data.
    CA::Grid::remove(ad.data_dir,setup.preproc_name+"_Grid","0");
  }

  if(setup.output_console)
    std::cout<<"Cleaned data"<< std::endl;


  if(setup.output_computation)
  {
    std::cout<<"-----------------" << std::endl; 
    std::cout<<"Total run time taken (s) = " << total_timer.millisecond()/1000.0 << std::endl;
    std::cout<<"-----------------" << std::endl; 
  }

  
  return 0;
}
