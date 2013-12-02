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

//! \file CADDIES2D.cpp
//! Perform the CADDIES2D flood modelling 
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2013-03


#include"ca2D.hpp"
#include"Masks.hpp"
#include"ArgsData.hpp"
#include"Setup.hpp"
#include"Rain.hpp"
#include"Inflow.hpp"
#include"WaterLevel.hpp"
#include"TimePlot.hpp"
#include"RasterGrid.hpp"


// Base on:
// CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
// C                                                                           C                
// C     Paper name                                                            C
// C     Authors                                                               C
// C     Date                                                                  C
// C     CA2D Version:1.0                                                      C
// C                                                                           C 
// CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


/* Round a number to the given digits.*/
inline CA::Real fround(CA::Real n, unsigned d)
{
  return std::floor(n * std::pow(static_cast<CA::Real>(10.0), static_cast<CA::Real>(d)) 
		    + static_cast<CA::Real>(.5)) 
    / std::pow(static_cast<CA::Real>(10.0), static_cast<CA::Real>(d) );
}


/* Extend a domain */
inline CA::Box extendBox(CA::Box extent, CA::Box& fullbox, CA::Unsigned lines)
{
  // Check that the new box is not too big by using full.
  if(extent.x()>fullbox.x()+lines) 
    extent.setX(extent.x()-lines);
  else
    extent.setX(fullbox.x());
  if(extent.y()>fullbox.y()+lines) 
    extent.setY(extent.y()-lines);	
  else
    extent.setY(fullbox.y());
  extent.setW(std::min((extent.w()+lines*2),fullbox.w()));
  extent.setH(std::min((extent.h()+lines*2),fullbox.h()));
  
  return extent;
}

//! Output to console the information about the simulation.
void outputConsole(CA::Unsigned iter, CA::Unsigned oiter, CA::Real t, CA::Real dt, 
		   CA::Real avgodt, CA::Real minodt, CA::Real maxodt, 
		   CA::Real vamax, CA::Real upstr_elv,
		   const Setup& setup)
{
  std::cout<<"-----"<<std::endl;
  std::cout<<"Total iterations = "<<iter<<" Simulation time (MIN) = "<<t/60
	   <<" Last DT = "<<dt<<std::endl;
  std::cout<<"Last iterations  = "<<oiter <<" Average DT ="<<avgodt/oiter
	   <<" Min DT = "<<minodt<<" Max DT = "<<maxodt<<std::endl;	
  std::cout<<"UPSTRELV = "<<upstr_elv<<std::endl;
  std::cout<<"VAMAX    = "<<vamax<<std::endl;    

  std::cout<<"-----"<<std::endl;

}


//! Compute the next time step as a fraction of the period time step.
//! Check that dt is between  min max.
void computeDT(CA::Real& dt, CA::Unsigned& dtfrac, CA::Real dtn1, const Setup& setup)
{
  CA::Unsigned dtfracmax = static_cast<CA::Unsigned>(setup.time_maxdt / setup.time_mindt);
  
  // Find the fraction of time_maxdt second that is just less than dtn1.
  // If dt is smaller than dtn1 we need to decrease the
  // fraction and stop as soon as the result is higher than dtn1.
  // If dt is bigger than dtn1 we need to increase the
  // fraction and stop as soon as the restuts is lower tha dtn1.
  if(dt<=dtn1)
  {
    for(dtfrac; dtfrac>=1;dtfrac--)
    {
      CA::Real tmpdt = setup.time_maxdt/static_cast<CA::Real>(dtfrac);
      if(tmpdt>=dtn1)
	break;
    }
    if(dtfrac==1)
      dt = setup.time_maxdt;
    else
    {
      dtfrac+=1;
      dt = setup.time_maxdt/static_cast<CA::Real>(dtfrac);
    }
  }
  else
  {
    for(dtfrac; dtfrac<=dtfracmax;dtfrac++)
    {
      CA::Real tmpdt = setup.time_maxdt/static_cast<CA::Real>(dtfrac);
      if(tmpdt<=dtn1)
	break;
    }
    dt = setup.time_maxdt/static_cast<CA::Real>(dtfrac);
  }

  // Check that dt is between  min max.
  dt = std::min(std::max(dt,setup.time_mindt), setup.time_maxdt);  
}


// -------------------------//
// Include the CA 2D functions //
// -------------------------//
#include CA_2D_INCLUDE(setBoundaryEle)
#include CA_2D_INCLUDE(outflowWCA2Dv1)
#include CA_2D_INCLUDE(waterdepthWCA2Dv1)
#include CA_2D_INCLUDE(velocityWCA2Dv1)
#include CA_2D_INCLUDE(removeUpstr)
#include CA_2D_INCLUDE(updatePEAKC)
#include CA_2D_INCLUDE(updatePEAKE)


int CADDIES2D(const ArgsData& ad, const Setup& setup, const CA::AsciiGrid<CA::Real>& eg, 
	      const std::vector<RainEvent>& res, const std::vector<WLEvent>& wles, 
	      const std::vector<IEvent>& ies, 
	      const std::vector<TimePlot>& tps, const std::vector<RasterGrid>& rgs)
{

  if(setup.output_computation)
  {
    std::cout<<"Simulation : "<<setup.sim_name<< std::endl;
    std::cout<<"Model      : "<<ad.model<< std::endl;
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
  CA::Grid  GRID(ad.data_dir,setup.short_name+"_Grid","0", ad.args.active());

  if(setup.output_console)
    std::cout<<"Loaded Grid data"<< std::endl;

  // Print the Grid Information.
  if(setup.output_console)
  {
    std::cout<<"-----------------" << std::endl; 
    GRID.printInfo(std::cout);
    std::cout<<"-----------------" << std::endl; 
  }

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

  // Create the computational domain, i.e the area where the water is
  // present and moving.
  CA::BoxList  compdomain;
  
  // -- INITIALISE ELEVATION ---
  
  // Create the elevation cell buffer.
  // It contains a "real" value in each cell of the grid.
  CA::CellBuffReal  ELV(GRID);

  // Set the border of the elevation buffer to be no data. 
  ELV.bordersValue(borders,eg.nodata);

  // Se the default value of the elevation to be nodata.
  ELV.fill(fulldomain, eg.nodata);

  // Load the data not from the DEM file but from the pre-processed
  // file.
  if(!ELV.loadData(setup.short_name+"_ELV","0") )
  {
    std::cerr<<"Error while loading the Elevation pre-processed file"<<std::endl;
    return 1;
  }

  // Highest elevation
  CA::Real     high_elv = 90000;

  // Find highest elevation
  ELV.sequentialOp(fulldomain, high_elv, CA::Seq::Max);

  if(setup.output_console)
  {
    std::cout<<"Loaded Elevation data"<< std::endl;
    std::cout<<"Highest elevation = "<< high_elv<<std::endl;
  }
  // ----  CELL BUFFERS ----
    
  // Create  the water depth cell buffer.
  CA::CellBuffReal WD(GRID);
   
  // Create the MASK cell buffer. The mask is usefull to check wich
  // cell has data and nodata and which cell has neighbourhood with
  // data.
  CA::CellBuffState MASK(GRID);

  // Create the velocity cell buffer with the speed of the velocity
  CA::CellBuffReal V(GRID);

  // Create the velocity cell buffer with the angle of the velocity
  CA::CellBuffReal A(GRID);


  // ----  EDGES BUFFERS ----

  // Create the outflow edge buffer.
  CA::EdgeBuffReal OUTF(GRID);

  // During the computation this store:
  // WCA2Dv1 the total flow for  an edge for an update step.
  CA::EdgeBuffReal TOT(GRID);
  

  // ---- ALARMS ----

  // Create the alarm(s) checked during the outflux computation.
  // Alarm 1: indicates when there is going to be an outflux outside of the computational domain.
  CA::Alarms  OUTFALARMS(GRID,1);

  // Create the alarm(s) checked during the velocity computation
  // Alarm 1: indicates when there is still water movement over the elevation threshould.
  CA::Alarms  VELALARMS(GRID,1);


  // ----  SCALAR VALUES ----

  CA::Unsigned iter   = 0;		   // The actual iteration number.
  CA::Real     t      = setup.time_start;  // The actual time in seconds
  CA::Real     dt     = setup.time_maxdt;  // Starting delta time.
  CA::Real     dtn1   = 0.0;
  CA::Real     nodata = eg.nodata;

  // The simulation time when the events will not produce any further water.
  CA::Real     t_end_events = setup.time_start; 

  // The level of water that can be ignored.
  CA::Real     ignore_wd  = setup.ignore_wd;

  // The water difference between cell that can be ignored.
  CA::Real     tol_delwl = setup.tolerance;

  // The minimum lelvel of water that is used for computing the
  // velocity. This cannot be to small (less than a centimer)
  CA::Real     tol_va    = std::max(setup.ignore_wd,static_cast<CA::Real>(0.01));

  // This is the period of when the velocity is computed.
  CA::Real     period_time_dt = setup.time_updatedt;
  CA::Real     time_dt        = t + period_time_dt;          // The next simulation time when to update dt.

 // The fraction of time_maxdt second used to compute the next the dt
  CA::Unsigned dtfrac    = 1;		  
  
  // The parameter of the time step.
  CA::Real     alpha  = setup.time_alpha;

  // The inverse of the roughness
  CA::Real     irough = 1/setup.roughness_global;
   
  CA::Unsigned oiter  = 0;                 // Output step num of iteration.
  CA::Real     minodt = setup.time_maxdt;  // Output step minimum dt.
  CA::Real     maxodt = 0.0;               // Output step maximum dt.
  CA::Real     avgodt = 0.0;	           // Average dt;
  CA::Real     time_output  = t + setup.output_period;	   // The time of the next output.

  CA::Real     rain_volume   = 0.0;
  CA::Real     inflow_volume = 0.0;
  CA::Real     wd_volume     = 0.0;

  // Maximum velocity.
  CA::Real     vamax=0.0;

  // Upstream elevation value. This value is the elevation where the
  // water "probably cannot reach anymore".
  CA::Real     upstr_elv = high_elv;

  // The total number of cells.
  CA::Real    total_cells = 0.0;

  // The potential velocity that an event has created. This should be
  // used to create muliple loop of FLUX calculation.
  CA::Real     potential_va;

  // If true the outflow computation will check the box cell for a
  // possible expanding domain.
  bool         checkbox = setup.expand_domain;         

  // Variable which indicates if the PEAK values need to be updated.
  // The PEAK values can be updated on every update step or on every
  // steps.
  bool UpdatePEAK    = false;

  // Variable which indicates if the raster habe been saved in the lat
  // iteration.
  bool RGwritten     = false;
  

  // -- CREATE FULL MASK ---
  
  CA::createCellMask(fulldomain,GRID,ELV,MASK,nodata);

  //***********************//
  //* BOUNDARY DISCUSSION *//
  //***********************//

  // The MASK is used to check if a cell with nodata has a neighbour
  // with data (Bit 31 set to true in the MASK). This kind of cells
  // are called boundary cells and they are ignored by the computation
  // (not visited), i.e. the rain is not added into them, the outflow
  // are not computed. However, these boundary cells are used as
  // neighbour cell thus the bouandary cell can have inflow (outflow
  // from the data cell). Thus the elevation value in these bounday
  // cell is changed with the Boundary Elevation value. If this value
  // is very high the global baundary is a CLOSED one, if it is VERY
  // negative the global boundary is an OPEN one.

  // ATTENTION The water that finish in the boundary cell (open
  // boundary case) is not removed (it stay in the WD buffer). 

  // Set the boundary cell elevation to the given boundary_elv value.
  CA::Execute::function(fulldomain, setBoundaryEle, GRID, ELV, MASK, setup.boundary_elv);

  //CA_DUMP_BUFF(ELV,0);  

  // ----  INIT WATER LEVEL EVENT  ----

  // Initialise the object that manage the water level events.
  WaterLevelManager wl_manager(GRID,wles);
  
  // Add the area with water level in the computational domain.
  wl_manager.addDomain(compdomain);

  // Get the elevation information.
  wl_manager.getElevation(ELV);

  // Analyse the area where a water level event will happen. WD
  // is used as temporary buffer.
  if(setup.check_vols)
    wl_manager.analyseArea(WD,MASK,fulldomain);

  // ----  INIT RAIN EVENT  ----

  // Initialise the object that manage the rain.
  RainManager rain_manager(GRID,res);

  // Add the area with rain in the computational domain.
  rain_manager.addDomain(compdomain);

  // Analyse the area where it will rain to use for volume cheking. WD
  // is used as temporary buffer.
  if(setup.check_vols)
    rain_manager.analyseArea(WD,MASK,fulldomain);

  // ----  INIT INFLOW EVENT  ----

  // Initialise the object that manage the inflow.
  InflowManager inflow_manager(GRID,ies);

  // Add the area with the inflow in the computational domain.
  inflow_manager.addDomain(compdomain);

  // Analyse the area where it will inflow to use for volume
  // cheking. WD is used as temporary buffer.
  if(setup.check_vols)
    inflow_manager.analyseArea(WD,MASK,fulldomain);
  
  // ----  INIT TIME PLOTS AND RASTER GRID ----

  std::string basefilename = ad.output_dir+ad.sdir+setup.short_name;

  // Initialise the object that manage the time plots.
  TPManager tp_manager(GRID,ELV,tps,basefilename,setup.timeplot_files);
  
  // Initialise the object that manage the time plots.
  RGManager rg_manager(GRID,rgs,basefilename,setup.rastergrid_files);    

  // List of raster grid data
  std::vector<RGData> rgdatas(rgs.size());
 
  // Peak buffers
  RGPeak rgpeak;

  for(size_t i = 0; i<rgs.size(); ++i)
  {
    std::string filename = ad.output_dir+ad.sdir+setup.short_name+"_"+setup.rastergrid_files[i]; 
    initRGData(filename, GRID, nodata, rgs[i], rgdatas[i], rgpeak);
  }
 
  // -- INITIALISE  ---

  // Clear the outflow buffer to zero (borders included).
  OUTF.clear();
  TOT.clear();

  // Set the wather depth to be zero.
  WD.fill(fulldomain, 0.0);

  // If there is not request to expand domain. Set the computtational and extend domain to full domain.
  if(!setup.expand_domain)
  {
    compdomain.clear();
    compdomain.add(fullbox);
  }

  // -- CALCULATE POSSIBLE INITIAL DT ---

  // Find the possible velocity caused by the events.
  potential_va = 0.0;
  potential_va = std::max( potential_va, rain_manager  .potentialVA(t,period_time_dt) );
  potential_va = std::max( potential_va, inflow_manager.potentialVA(t,period_time_dt) );
  potential_va = std::max( potential_va, wl_manager    .potentialVA(t,period_time_dt) );

  // Compute the possible next time step using the critical velocity equations.
  dtn1 = std::min(setup.time_maxdt,alpha*GRID.length()/potential_va);
  
  // Compute the next time step as a fraction fo the period time step.
  // Check that dt is between  min max.
  computeDT(dt,dtfrac,dtn1,setup);

  // -- PREPARE THE EVENTS MANAGERS WITH THE NEW DT ---

  // Get the amount of events that would happen in each area for each dt
  // for the next period.
  rain_manager  .prepare(t,period_time_dt,dt);
  inflow_manager.prepare(t,period_time_dt,dt);
  wl_manager    .prepare(t,period_time_dt,dt);
  
  // -- PREMATURE END --

  // get the time when the events will not add any further water.
  t_end_events = std::max(t_end_events, rain_manager.endTime());
  t_end_events = std::max(t_end_events, inflow_manager.endTime());
  t_end_events = std::max(t_end_events, wl_manager.endTime());

  if(setup.output_computation)
  {
    std::cout<<"The events will end at "<<t_end_events<<" (s) simulation time"<< std::endl;
    std::cout<<"------------------------------------------" << std::endl; 
  }

  // ------------------------- MAIN LOOP -------------------------------
  while(iter<setup.time_maxiters && t<setup.time_end)
  {
    // Set this to false. This will be set to the righ value during an
    // update step or before the update itself.
    UpdatePEAK = false;

    // The raster have not been written this itertaion.
    RGwritten     = false;

    // If there is the request to expand the domain.
    // Deactivate Box alarm(s) and set them.
    if(setup.expand_domain)
    {
      OUTFALARMS.deactivateAll();
      OUTFALARMS.set();
    }

    // --- CONSOLE OUTPUT ---

    // Check if it is time to output to console.
    if(setup.output_console && t>= time_output)
    {
      outputConsole(iter,oiter,t,dt,avgodt,minodt,maxodt,vamax,upstr_elv,setup);

      oiter  = 0;
      avgodt = 0.0;
      minodt = setup.time_maxdt;  // Output step minimum dt.
      maxodt = 0.0;               // Output step maximum dt.
      
      // Compute the next output time.
      time_output += setup.output_period; 

      if(setup.check_vols == true)
      {
	// Compute the total volume of water that is in the water
	// depth (included the boundary cell).
	WD.sequentialOp(fulldomain, wd_volume, CA::Seq::Add);
	wd_volume *= GRID.length()*GRID.length();

	std::cout<<"Volume RAIN = "<<rain_volume<<" Volume INFLOW = "<<inflow_volume
		 <<" Volume WD = "<<wd_volume<<std::endl;	
	std::cout<<"-----------------" << std::endl; 
      }

      if(setup.output_computation)
      {
	std::cout<<"Partial run time taken (s) = " << total_timer.millisecond()/1000.0 << std::endl;
	std::cout<<"-----------------" << std::endl; 
      }
    }

    // --- SIMULATION TIME ---

    // Set the new time step.
    t+=dt;

    // Round the time step to be of 0.01 second precision and then
    // check if it is a multiple of 60 (with a 0.01 precision). If it
    // is the case we need to reset the time of the simulation. If we
    // don't do this the floating point error will creep into the time
    // of simulation.
    CA::Real tround = fround(t, 2);
    if(std::fmod(tround,period_time_dt) < static_cast<CA::Real>(0.01))
    {
      t = tround;
    }

    // Compute output step information.
    avgodt+=dt;
    
    if(dt>maxodt) maxodt=dt;
    if(dt<minodt) minodt=dt;    

    // --- COMPUTE OUTFLUX ---

    // Clear the outflow buffer to zero (borders included).
    OUTF.clear();
          
    // Compute outflow using WCA2Dv1.
    // Check if there is an outflow in the border of the box.
    CA::Execute::function(compdomain, outflowWCA2Dv1, GRID, OUTF, ELV, WD, MASK, OUTFALARMS,
			  ignore_wd, tol_delwl,dt, irough);

    // If there is a request to expand the domain.
    // Get the alarms states.
    if(setup.expand_domain)
    {
      OUTFALARMS.get();
    
      // Check if the box alarm is active, that mean there is some
      // outflux on the border of the computational domain.
      if(OUTFALARMS.isActivate(0))
      {
	// Set the computational domain to be the extend version and
	// create the new extended one.
	CA::Box extent(compdomain.extent());
	compdomain.clear();
	compdomain.add(extendBox(extent, fullbox,1));
      }  
    }
    
    // --- UPDATE WL AND WD  ---
    
    // Update the water depth with the outflux and store the total
    // amount of outflux for the WCA2Dv1 model. 
    CA::Execute::function(compdomain, waterdepthWCA2Dv1, GRID, WD, OUTF, TOT, MASK, dt, period_time_dt);

    // --- EXTRA LATERAL EVENT(s) ---

    // Add the eventual rain events.
    rain_manager.add(WD,MASK,t,dt);

    // Add the eventual inflow events.
    inflow_manager.add(WD,MASK,t,dt);

    // Add the eventual water level events.
    wl_manager.add(WD,ELV,MASK,t,dt);

    // --- COMPUTE NEXT DT, I.E. PERIOD STEP ---
    
    // Check if the dt need to be re-computed.
    if(t>=time_dt)
    {   
      if(setup.ignore_upstream)
      {
	// Deactivate the alarms checked during the velocity
	// calculation.
	VELALARMS.deactivateAll();
	VELALARMS.set();
      }

      // Lets make sure there are not any rounding errors.
      t = time_dt;

      // The peak value need to be updated
      UpdatePEAK = true;

      // Update the total volume from the events for the last period.
      rain_volume   += rain_manager.volume();
      inflow_volume += inflow_manager.volume();
      // NO water level at the moment.

      // --- UPDATE VA  ---
      
      // Compute the velocity using the total outflux.
      // Attention the tollerance is different here. 
      // Check if there is water movement over the upstream elevation threshould.
      CA::Execute::function(compdomain, velocityWCA2Dv1, GRID, V, A, WD, ELV, TOT, MASK, VELALARMS,
			    tol_va, period_time_dt, irough, upstr_elv);

      // CLear the total outflux.
      TOT.clear();
                  
      // Retrieve the maximum velocity 
      V.sequentialOp(compdomain, vamax,CA::Seq::MaxAbs);                

      // Find the maximum velocity
      CA::Real grid_max_va = vamax;

      // Find the possible velocity caused by the events.
      potential_va = 0.0;
      potential_va = std::max(potential_va, rain_manager  .potentialVA(t,period_time_dt) );
      potential_va = std::max(potential_va, inflow_manager.potentialVA(t,period_time_dt) );
      potential_va = std::max(potential_va, wl_manager    .potentialVA(t,period_time_dt) );

      // Compute the possible next dt from the grid velocity and from
      // the potential velocity fron an event.
      dtn1 = setup.time_maxdt;
      dtn1 = std::min(dtn1,alpha*GRID.length()/grid_max_va);
      dtn1 = std::min(dtn1,alpha*GRID.length()/potential_va);
            	     
      // Compute the next time step as a fraction fo the period time step.
      // Check that dt is between  min max.
      computeDT(dt,dtfrac,dtn1,setup);
      
      // When the dt need to be recomputed.
      time_dt += period_time_dt;
    
      // Prepare the Events manager for the next update .
      rain_manager  .prepare(t,period_time_dt,dt);
      inflow_manager.prepare(t,period_time_dt,dt);
      wl_manager    .prepare(t,period_time_dt,dt);

      if(setup.ignore_upstream)
      {
	// Check the ALARMS computed during the velocity step.
	VELALARMS.get();
	
	// If the alarms one is not active, then there was not any cell
	// which had the water level over the upstream elevation
	// threshold and some flux at the same time. So we can remove
	// the cell from the computation and then lower the upstream
	// threshold.
	// ATTENTION This action is performed only when all the events
	// finished to add water to the domain.
	if(!VELALARMS.isActivate(0) && t > t_end_events)
	{
	  CA::Execute::function(fulldomain, removeUpstr, GRID, MASK, ELV, upstr_elv);
	  upstr_elv -= setup.upstream_reduction;
	}
      }
    } // COMPUTE NEXT DT.

    // -------  OUTPUTS --------
         
    // Output time plots.
    tp_manager.output(t, iter, WD, V, setup.output_console);

    // Check if we need to update the peak at every time step.
    if(setup.update_peak_dt)
      UpdatePEAK = true;
    
    // Update the peak
    if(UpdatePEAK)
      rg_manager.updatePeak(compdomain,WD,V,MASK);

    // Output raster grid. Keep track if the raster have been written.
    RGwritten = rg_manager.output(t, WD, V, A, setup.short_name, setup.output_console);
	      
    // ---- END OF ITERATION ----

    // Increase time step (and output time step)
    iter++;    
    oiter++;
  }  

  // --- END OF MAIN LOOP ---

  // --- OUTPUT PEAK ---
  
  // Check if the raster have not been written in the last
  // iteration. If not, we need to make sure that we save the PEAK.
  if(!RGwritten)
  {
    // Make sure to output the last peack value.  
    rg_manager.updatePeak(compdomain,WD,V,MASK);
    rg_manager.outputPeak(t,WD, V, setup.short_name, setup.output_console);
  }

  // --- CONSOLE OUTPUT ---
  
  // Check if it is time to output to console.
  if(setup.output_console && t>= time_output)
  {    
    outputConsole(iter,oiter,t,dt,avgodt,minodt,maxodt,vamax,upstr_elv,setup);

    if(setup.check_vols == true)
    {
      // Compute the total volume of water that is in the water
      // depth (included the boundary cell).
      WD.sequentialOp(fulldomain, wd_volume, CA::Seq::Add);
      wd_volume *= GRID.length()*GRID.length();
      
      std::cout<<"Volume RAIN = "<<rain_volume<<" Volume INFLOW = "<<inflow_volume
	       <<" Volume WD = "<<wd_volume<<std::endl;	
      std::cout<<"-----------------" << std::endl;       
    }
  }

  // ---- TIME OUTPUT ----

  if(setup.output_computation)
  {
    std::cout<<"-----------------" << std::endl; 
    std::cout<<"Total run time taken (s) = " << total_timer.millisecond()/1000.0 << std::endl;
    std::cout<<"-----------------" << std::endl; 
  }

  
  return 0;
}
