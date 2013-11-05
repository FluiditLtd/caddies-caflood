//! \file WCA2D.cpp
//! Perform the Weighted CA2D flood modelling 
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
#include CA_2D_INCLUDE(setBoundaryEle)
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
  
  // ----  INIT TIME PLOTS ----

  
  // Initialise the object that manage the time plots.
  std::string basefilename = ad.output_dir+ad.sdir+setup.short_name;
  TPManager tp_manager(GRID,ELV,tps,basefilename,setup.timeplot_files);
  
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
    rain_manager.add(WD,MASK,t,dt);

    // Add the eventual inflow events.
    inflow_manager.add(WD,MASK,t,dt);

    // Add the eventual water level events.
    wl_manager.add(WD,ELV,MASK,t,dt);

    // --- COMPUTE NEXT DT, I.E. PERIOD STEP ---
    
    // Check if the dt need to be re-computed.
    if(t>=time_dt)
    {   
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

    } // COMPUTE NEXT DT.

    // -------  OUTPUTS --------
    
    // Output time plots.
    tp_manager.output(t, iter, WD, V, setup.output_console);


    // ----  RASTER GRID ----

    RGoutputed    = false;

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
