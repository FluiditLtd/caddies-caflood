//! \file WCA2D.cpp
//! Perform the Weighted CA2D flood modelling 
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2013-03


#include"ca2D.hpp"
#include"Masks.hpp"
#include"ArgsData.hpp"
#include"Setup.hpp"
#include"Rain.hpp"
#include"Events.hpp"
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
		   CA::Real vamax, 
		   const Setup& setup)
{
  std::cout<<"-----"<<std::endl;
  std::cout<<"Total iterations = "<<iter<<" Simulation time (MIN) = "<<t/60
	   <<" Last DT = "<<dt<<std::endl;
  std::cout<<"Last iterations  = "<<oiter <<" Average DT ="<<avgodt/oiter
	   <<" Min DT = "<<minodt<<" Max DT = "<<maxodt<<std::endl;	
  std::cout<<"VAMAX= "<<vamax<<std::endl;    

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
#include CA_2D_INCLUDE(computeArea)
#include CA_2D_INCLUDE(computeCells)
#include CA_2D_INCLUDE(setBoundaryEle)
#include CA_2D_INCLUDE(addRain)
#include CA_2D_INCLUDE(addRaise)
#include CA_2D_INCLUDE(addInflow)
#include CA_2D_INCLUDE(outflowWeightedWDv2)
#include CA_2D_INCLUDE(computeWDv2)
#include CA_2D_INCLUDE(computeVelocityVAv1)
#include CA_2D_INCLUDE(updatePEAKC)
#include CA_2D_INCLUDE(updatePEAKE)


int WCA2D(const ArgsData& ad, const Setup& setup, const CA::AsciiGrid<CA::Real>& eg, 
	  const std::vector<RainEvent>& res, const std::vector<WLEvent>& wles, const std::vector<IEvent>& ies, 
	  const std::vector<TimePlot>& tps, const std::vector<RasterGrid>& rgs)
{

  if(setup.output_computation)
  {
    std::cout<<"Simulation : "<<setup.sim_name<< std::endl;
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

  if(setup.output_console)
    std::cout<<"Loaded Elevation data"<< std::endl;

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

  // During the computation this store the average amount of flow of
  // an edge for an update step.
  CA::EdgeBuffReal AVGOUTF(GRID);
  
  // ---- LOOKUP TABLE ----

  // Create the lookup table for the pow(R,2/3) used in the manning
  // equation.
  //CA::LookupTableReal POWR23(GRID,0,4,0.01);


  // ---- ALARMS ----

  // Create the alarm(s).
  // Alarm 1: indicates when there is going to be an outflux outside of the computational domain.
  CA::Alarms  ALARMS(GRID,1);


  // ----  SCALAR VALUES ----

  CA::Unsigned iter   = 0;		   // The actual iteration number.
  CA::Real     t      = setup.time_start;  // The actual time in seconds
  CA::Real     dt     = setup.time_maxdt;  // Starting delta time.
  CA::Real     dtn1   = 0.0;
  CA::Real     nodata = eg.nodata;

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

  // The area of the base queare grid cell.
  CA::Real     cell_area  = GRID.length()*GRID.length();

  // The total number of cells.
  CA::Real    total_cells = 0.0;

  // The potential velocity that an event has created. This should be
  // used to create muliple loop of FLUX calculation.
  CA::Real     potential_va;

  // If true the outflow computation will check the box cell for a
  // possible expanding domain.
  bool         checkbox = setup.expand_domain;         

  // These variables are used to indicates if the output to console
  // happen in the case of time plot or raster grid output.
  bool TPoutputed    = false;
  bool RGoutputed    = false;

  // Since WL variable does not exist. We need to save WD and then use
  // ELV to compute WL. The following variable is used to save WD only
  // once if both WD and WL are requested. The same is the case for
  // VEL, which needs WD.
  bool VAsaved       = false;
  bool VAPEAKsaved   = false;
  bool VAPEAKupdated = false;
  bool WDsaved       = false;
  bool WDPEAKsaved   = false;
  bool WDPEAKupdated = false;


  // Update the peak only every update time.
  bool UpdatePEAK    = false;
  
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

  // ---- NUMBER OF CELLS. ---

  // Find the total number of cells.
  // WD is used as temporary buffer here.
  WD.fill(fulldomain, 0.0);
  CA::Execute::function(fulldomain, computeCells, GRID, WD, MASK);    
  WD.sequentialOp(fulldomain, total_cells, CA::Seq::Add);	


  // ----  INIT WATER LEVEL EVENT  ----

  // List of rain event data.
  std::vector<WLEData> wledatas(wles.size());

  for(size_t i = 0; i<wles.size(); ++i)
  {
    initWLEData(GRID, wles[i], wledatas[i]);

    // Retrieve the minimum elevation of the given area.
    ELV.sequentialOp(wledatas[i].box_area, wledatas[i].min_elv, CA::Seq::MinAbs);	

    // Add the area with water level rise in the computational domain.
    compdomain.add(wledatas[i].box_area);
  }

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

  // List of inflow event data.
  std::vector<IEData> iedatas(ies.size());

  for(size_t i = 0; i<ies.size(); ++i)
  {
    initIEData(GRID, ies[i], iedatas[i]);

    // Compute area to use for volume cheking.
    if(setup.check_vols)
    {
      // WD is used as temporary buffer here.
      WD.fill(fulldomain, 0.0);
      CA::Execute::function(iedatas[i].box_area, computeArea, GRID, WD, MASK);    
      WD.sequentialOp(iedatas[i].box_area, iedatas[i].grid_area, CA::Seq::Add);	
    }

    // Add the area with inflow in the computational domain.
    compdomain.add(iedatas[i].box_area);
  }
  
  // ----  INIT TIME PLOTS ----

  // List of time plots data.
  std::vector<TPData> tpdatas(tps.size());
  
  for(size_t i = 0; i<tps.size(); ++i)
  {
    std::string filename = ad.output_dir+ad.sdir+setup.short_name+"_"+setup.timeplot_files[i];
    initTPData(filename, GRID, ELV, tps[i], tpdatas[i]);
  }

  // ----  INIT RASTER GRID ----

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
  AVGOUTF.clear();

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
  potential_va = std::max( potential_va, rain_manager.potentialVA(t,period_time_dt) );

  // Compute the possible next time step using the critical velocity equations.
  dtn1 = std::min(setup.time_maxdt,alpha*GRID.length()/potential_va);
  
  // Compute the next time step as a fraction fo the period time step.
  // Check that dt is between  min max.
  computeDT(dt,dtfrac,dtn1,setup);

  // -- PREPARE THE EVENTS MANAGERS WITH THE NEW DT ---

  // Get the amount of rain tahts houdl fall in each area for each dt
  // for the next period.
  rain_manager.prepare(t,period_time_dt,dt);
  

  // ------------------------- MAIN LOOP -------------------------------
  while(iter<setup.time_maxiters && t<setup.time_end)
  {
    // If there is tehe request to expand the domain.
    // Deactivate Box alarm(s) and set them.
    if(setup.expand_domain)
    {
      ALARMS.deactivateAll();
      ALARMS.set();
    }

    // --- CONSOLE OUTPUT ---

    // Check if it is time to output to console.
    if(setup.output_console && t>= time_output)
    {
      outputConsole(iter,oiter,t,dt,avgodt,minodt,maxodt,vamax,setup);

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
	wd_volume *= cell_area;

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

    // The potential velocity that an event can create is set to
    // zero. This value is used to create eventual muliple loop of
    // FLUX calculation. ATTENTION, at the moment is used only to
    // compute potential_dt.
    potential_va = 0;

    // --- INFLOW EVENT(s) ---

    // Loop through the inflow event(s).
    for(size_t i = 0; i<ies.size(); ++i)
    {
      size_t index = iedatas[i].index;

      // If the index is larger than the available ins/time, do
      // nothing.
      if(index >= ies[i].ins.size() )
	continue;
      
      // Compute the inflow volume at specific time using
      // interpolation. Check if the index is the last available
      // one, then there is no inflow.
      CA::Real volume = 0;
      if(index != ies[i].ins.size() -1)
      {
	CA::Real y0 = ies[i].ins[index];
	CA::Real y1 = ies[i].ins[index+1];
	CA::Real x0 = ies[i].times[index];
	CA::Real x1 = ies[i].times[index+1];
	CA::Real t0 = t - dt;
	CA::Real t1 = t;
	CA::Real yt0= y0 + (y1-y0) * ( (t0 - x0)/(x1 - x0) );
	CA::Real yt1= y0 + (y1-y0) * ( (t1 - x0)/(x1 - x0) );
	volume = 0.5*(t1-t0)*(yt1-yt0)+(t1-t0)*(yt0);
      }

      // If it is requested to check the volumes, compute the total
      // volume of inflow 
      if(setup.check_vols == true)
      {
	iedatas[i].volume += volume;
	inflow_volume     += volume;	
      }

      // ATTENTION The volume is the total volume, it need to be
      // divided by the number of cells that are going to receive the
      // inflow. 
      volume = volume/(iedatas[i].grid_area/cell_area);

      // Compute the potential velocity using the amount of extra
      // water level added.
      potential_va = std::max(potential_va, std::sqrt(static_cast<CA::Real>(9.81)*(volume/cell_area) ) );
      
      // Add (or subtract) the given volume into the water detph of the
      // given area.
      CA::Execute::function(iedatas[i].box_area, addInflow, GRID, WD, MASK, volume);           
      
      // Check if the simulation time now is equal or higher than the
      // time of the NEXT index.
      if(t >= ies[i].times[index+1])
	index++;
      
      // Update index.
      iedatas[i].index = index;
    }

    // --- WATER LEVEL EVENT(s) ---

    // Loop through the water level event(s).
    for(size_t i = 0; i<wles.size(); ++i)
    {
      size_t index = wledatas[i].index;

      // If the index is larger than the available rain/time, do
      // nothing.
      if(index >= wles[i].wls.size() )
	continue;
      
      // Compute the water level at specific are using
      // interpolation. Check if the index is the last available
      // one. In this case use only one value.
      CA::Real level      = 0;
      if(index == wles[i].wls.size() -1)
      {
	level = wles[i].wls[index];
      }
      else
      {	
	CA::Real y0 = wles[i].wls[index];
	CA::Real y1 = wles[i].wls[index+1];
	CA::Real x0 = wles[i].times[index];
	CA::Real x1 = wles[i].times[index+1];
	level = y0 + (y1-y0) * ( (t - x0)/(x1 - x0) );
      }
      
      // Compute the potential velocity using the maximum water depth
      // in the area, i.e. the level minus the minimum elevation
      CA::Real  water_depth = level - wledatas[i].min_elv;
      potential_va = std::max(potential_va, std::sqrt( static_cast<CA::Real>(9.81)*( water_depth) ) );

      // Given the way the CA2D model work, we need to set the water
      // depth instead of the water level. Thus the water depth value
      // at specific location is the value of the water level event
      // minus the elevation.
      CA::Execute::function(wledatas[i].box_area, addRaise, GRID, WD, ELV, MASK, level);     

      // Check if the simulation time now is equal or higher than the
      // time of the NEXT index.
      if(t >= wles[i].times[index+1])
	index++;
      
      // Update index.
      wledatas[i].index = index;
    }
    
    // --- COMPUTE OUTFLUX ---

    // Clear the outflow buffer to zero (borders included).
    OUTF.clear();
          
    // Compute outflow using weighted method version 2.
    // This version check if there is an outflow in the border of the box.
    CA::Execute::function(compdomain, outflowWeightedWDv2, GRID, OUTF, ELV, WD, MASK, ALARMS,
			  ignore_wd, tol_delwl,dt, irough);

    // If there is a request to expand the domain.
    // Get the alarms states.
    if(setup.expand_domain)
    {
      ALARMS.get();
    
      // Check if the box alarm is active, that mean there is some
      // outflux on the border of the computational domain.
      if(ALARMS.isActivate(0))
      {
	// Set the computational domain to be the extend version and
	// create the new extended one.
	CA::Box extent(compdomain.extent());
	compdomain.clear();
	compdomain.add(extendBox(extent, fullbox,1));
      }  
    }
    
    // --- UPDATE WL AND WD  ---
    
    // Update the water depth with the outflux and store the average
    // amount of outflux. 
    CA::Execute::function(compdomain, computeWDv2, GRID, WD, OUTF, AVGOUTF, MASK, dt, period_time_dt);

    // --- EXTRA LATERAL EVENT(s) ---

    // Add the eventual rain events.
    rain_manager.add(WD,MASK);


    // --- COMPUTE NEXT DT, I.E. PERIOD STEP ---
    
    // Check if the dt need to be re-computed.
    if(t>=time_dt)
    {   
      // Lets make sure there are not any rounding errors.
      t = time_dt;

      // The peak value need to be updated
      UpdatePEAK = true;

      // --- UPDATE VA  ---
      
      // Compute the velocity using the total outflux.
      // Attention the tollerance is different here. 
      CA::Execute::function(compdomain, computeVelocityVAv1, GRID, V, A, WD, ELV, AVGOUTF, MASK, 
			    tol_va, period_time_dt, irough);

      // CLear the average outflux. This improve the computed
      // velocity.
      AVGOUTF.clear();
                  
      // Retrieve the maximum velocity 
      V.sequentialOp(compdomain, vamax,CA::Seq::MaxAbs);                

      // Find the maximum velocity
      CA::Real grid_max_va = vamax;

      // Find the possible velocity caused by the events.
      potential_va = 0.0;
      potential_va = std::max(potential_va, rain_manager.potentialVA(t,period_time_dt) );

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
      rain_manager.prepare(t,period_time_dt,dt);

    } // COMPUTE NEXT DT.

    // -------  OUTPUTS --------
    
    TPoutputed    = false;
    RGoutputed    = false;

    // ----  TIME PLOTS ----
    
    for(size_t i = 0; i<tpdatas.size(); ++i)
    {
      // Check if it is time to plot and the file is good.
      if(t >= tpdatas[i].time_next && tpdatas[i].file->good())
      {
	if(!TPoutputed && setup.output_console)
	{
	  std::cout<<"Update Time Plot :";
	  TPoutputed = true;
	}
	
	switch(tps[i].pv)
	{
	case PV::VEL:
	  {
	    if(setup.output_console)
	      std::cout<<" VEL";

	    // Retrieve the speed
	    V.retrievePoints(tpdatas[i].pl,&(tpdatas[i].pvals[0]),tpdatas[i].pl.size());      

	    (*tpdatas[i].file)<<iter<<", "<<t/60.0<<", ";	  
	    // Write the speed
	    for(CA::Unsigned p =0; p< tpdatas[i].pl.size(); p++)
	    {	
	      (*tpdatas[i].file)<<tpdatas[i].pvals[p]<<", ";
	    }
	    (*tpdatas[i].file)<<std::endl;
	  }
	  break;
	case PV::WL:
	  {
	    if(setup.output_console)
	      std::cout<<" WL";

	    // Retrieve the water depth
	    WD.retrievePoints(tpdatas[i].pl,&(tpdatas[i].pvals[0]),tpdatas[i].pl.size());      

	    (*tpdatas[i].file)<<iter<<", "<<t/60.0<<", ";	  
	    // Write the water level by adding the previously saved elevation.
	    for(CA::Unsigned p =0; p< tpdatas[i].pl.size(); p++)
	    {	
	      (*tpdatas[i].file)<<tpdatas[i].pelvs[p] + tpdatas[i].pvals[p]<<", ";
	    }
	    (*tpdatas[i].file)<<std::endl;
	  }
	  break;
	case PV::WD:
	  {
	    if(setup.output_console)
	      std::cout<<" WD";

	    WD.retrievePoints(tpdatas[i].pl,&(tpdatas[i].pvals[0]),tpdatas[i].pl.size());      

	    (*tpdatas[i].file)<<iter<<", "<<t/60.0<<", ";	  
	    for(CA::Unsigned p =0; p< tpdatas[i].pl.size(); p++)
	    {	
	      (*tpdatas[i].file)<<tpdatas[i].pvals[p]<<", ";
	    }
	    (*tpdatas[i].file)<<std::endl;
	  }
	  break;
	}
	
	// Update the next time to plot.
	tpdatas[i].time_next += tps[i].period;
      }
    }

    if(TPoutputed && setup.output_console)
      std::cout<<std::endl;

    // ----  RASTER GRID ----

    VAsaved       = false;
    VAPEAKsaved   = false;
    VAPEAKupdated = false;
    WDsaved       = false;
    WDPEAKsaved   = false;
    WDPEAKupdated = false;
      
    // ----  UPDATE PEAK ----

    // Update the peak only one there is a update period step.
    if(UpdatePEAK || setup.update_peak_dt)
    {
      UpdatePEAK = false;

      for(size_t i = 0; i<rgdatas.size(); ++i)
      {
	// Check if the peak values need to be updated.
	if(rgs[i].peak == true)
	{
	  switch(rgs[i].pv)
	  {
	  case PV::VEL:
	    // Update the absolute maximum velocity.
	    // ATTENTION need to be tested.
	    if(!VAPEAKupdated)
	    {
	      CA::Execute::function(compdomain, updatePEAKC, GRID, (*rgpeak.V), V, MASK);
	      VAPEAKupdated = true;
	    }
	    // ATTENTION! The break is removed since in order to
	    // post-process VA we need WD. 
	    // break;	  
	  case PV::WL:
	  case PV::WD:	  
	    // Update the absolute maximum water depth only once.
	    if(!WDPEAKupdated)
	    {
	      CA::Execute::function(compdomain, updatePEAKC, GRID, (*rgpeak.WD), WD, MASK);
	      WDPEAKupdated = true;
	    }
	    break;
	  }
	}
      }
    }

    // ----  WRITE GRID ----

    for(size_t i = 0; i<rgdatas.size(); ++i)
    {
      // Check if it is time to plot!
      if(t >= rgdatas[i].time_next)
      {	
	if(!RGoutputed && setup.output_console)
	{
	  std::cout<<"Write Raster Grid: ";
	  RGoutputed = true;
	}

	// Retrieve the string of the time.
	std::string strtime;
	CA::toString(strtime,std::floor(t+0.5));

	// Save the buffer using direct I/O where the main ID is the
	// buffer name and the subID is the timestep.
	switch(rgs[i].pv)
	{
	case PV::VEL:	  
	  if(!VAsaved)
	  {
	    if(setup.output_console)
	      std::cout<<" VA";	  
	    V.saveData(setup.short_name+"_V",strtime);
	    A.saveData(setup.short_name+"_A",strtime);
	    VAsaved = true;
	  }
	  // ATTENTION! The break is removed since in order to
	  // post-process VA we need WD. 
	  //break;
	case PV::WL:
	case PV::WD:
	  if(!WDsaved)
	  {
	    if(setup.output_console)
	      std::cout<<" WD";
	  
	    WD.saveData(setup.short_name+"_WD",strtime);
	    WDsaved = true;
	  }
	  break;
	}
	
	// Check if the peak values need to be saved.
	if(rgs[i].peak == true)
	{
	  // Non velocity raster grid.
	  switch(rgs[i].pv)
	  {
	  case PV::VEL:	  
	    if(!VAPEAKsaved)
	    {
	      if(setup.output_console)
		std::cout<<" VAPEAK";
	  
	      rgpeak.V->saveData(setup.short_name+"_V","PEAK");
	      VAPEAKsaved = true;
	    }
	    // ATTENTION! The break is removed since in order to
	    // post-process VA we need WD. 
	    // break;
	  case PV::WL:
	  case PV::WD:
	    if(!WDPEAKsaved)
	    {
	      if(setup.output_console)
		std::cout<<" WDPEAK";
	  
	      rgpeak.WD->saveData(setup.short_name+"_WD","PEAK");
	      WDPEAKsaved = true;
	    }
	    break;
	  }
	  
	}
	// Update the next time to save a raster grid.
	rgdatas[i].time_next += rgs[i].period;
      }
    }

    if(RGoutputed && setup.output_console)
      std::cout<<std::endl;
	      
    // ---- END OF ITERATION ----

    // Increase time step (and output time step)
    iter++;    
    oiter++;
  }  

  // ----  PEAK RASTER GRID ----

  RGoutputed    = false;

  for(size_t i = 0; i<rgdatas.size(); ++i)
  {
    // Check if the peak values need to be saved.
    if(rgs[i].peak == true)
    {

      if(!RGoutputed && setup.output_console)
      {
	std::cout<<"Write Raster Grid: ";
	RGoutputed = true;
      }

      // Non velocity raster grid.
      switch(rgs[i].pv)
      {
      case PV::VEL:	  
	if(!VAPEAKsaved)
	{
	  if(setup.output_console)
	    std::cout<<" VAPEAK";
	  
	  rgpeak.V->saveData(setup.short_name+"_V","PEAK");
	  VAPEAKsaved = true;
	}
	// ATTENTION! The break is removed since in order to
	// post-process VA we need WD. 
	// break;
      case PV::WL:
      case PV::WD:
	if(!WDPEAKsaved)
	{
	  if(setup.output_console)
	    std::cout<<" WDPEAK";
	  
	  rgpeak.WD->saveData(setup.short_name+"_WD","PEAK");
	  WDPEAKsaved = true;
	}
	break;
      }
    }
  }

  if(RGoutputed && setup.output_console)
    std::cout<<std::endl;
	  
  // --- CONSOLE OUTPUT ---
  
  // Check if it is time to output to console.
  if(setup.output_console && t>= time_output)
  {    
    outputConsole(iter,oiter,t,dt,avgodt,minodt,maxodt,vamax,setup);

    if(setup.check_vols == true)
    {
      // Compute the total volume of water that is in the water
      // depth (included the boundary cell).
      WD.sequentialOp(fulldomain, wd_volume, CA::Seq::Add);
      wd_volume *= cell_area;
      
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
