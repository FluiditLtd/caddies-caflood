#ifndef _SETUP_HPP_
#define _SETUP_HPP_


//! \file Setup.hpp
//!  Contains the structure with the main setup for the CA2D model.
//! \author Michele Guidolin, University of Exeter, 
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2013-03

#include"BaseTypes.hpp"
#include<string>
#include<vector>


//! Structure with the data that define the main variables used to
//! execute the CA2D model. 
struct Setup
{
  //  --- SIMULATION  --- 
  std::vector<CA::Unsigned> sim_version; 
  std::string    sim_name;	//!< Name of the simulation. 
  std::string    short_name;	//!< Short name used as base for output files. 

  //  --- TIME VALUES  ---  
  CA::Real      time_start;	//!< Starting time of the simulation (seconds). 
  CA::Real      time_end;	//!< Ending time of the simulation (seconds). 
  CA::Real      time_maxdt;	//!< The initial and maximum time step.
  CA::Real      time_mindt;	//!< The minimum time step.
  CA::Real      time_updatedt;	//!< The time step to update the time step.
  CA::Real      time_alpha;	//!< The proportion of adaptive time step taken. 
  CA::Unsigned  time_maxiters;	//!< The maximum number of iterations. 

  //  --- SIMULATION PARAMETERS ---
  CA::Real      roughness_global;  //!< Global default value.  
  CA::Real      tolerance;	   //!< The water difference between cell that can be ignored.
  CA::Real      ignore_wd;	   //!< The water depth that can be ignored.

  //! The elevation of boundary cell, a high value represent a CLOSED
  //! boundary a low value represents an OPEN boundary.
  CA::Real      boundary_elv;	   

  //  --- ELEVATION  ---  
  //! ARC/INFO ASCII GRID format file with the specific
  //! elavation value for each cell and the no data value.
  std::string   elevation_ASCII; 
  
  //  --- RAIN EVENT  ---  
  //! CSV file(s) with the configuration of the rain event(s) to add
  //! into the grid.
  std::vector<std::string> rainevent_files;

  //  --- WATER LEVEL EVENT  ---  
  //! CSV file(s) with the configuration of the water level event(s)
  //! to add into the grid.
  std::vector<std::string> wlevent_files;

  //  --- INFLOW EVENT  ---  
  //! CSV file(s) with the configuration of the inflow event(s)
  //! to add into the grid.
  std::vector<std::string> inflowevent_files;

  
  //  --- TIMEPLOT  ---  
  //! CSV file(s) with the configuration of the time plot(s) to
  //! produce, i.e. the value of a physical variable (WD,WL,VEL) at
  //! specific time and location.
  std::vector<std::string> timeplot_files;

  //  --- RASTERS  ---  
  //! CSV file(s) with the configuration of the raster grid(s) to
  //! produce, i.e. the grid value of a physical variable (WD,WL,VEL) at
  //! specific time.
  std::vector<std::string> rastergrid_files;

  
  //  --- OUTPUT  ---  
  bool output_console;		//!< If true output info into console.
  CA::Real output_period;	//!< The period to ouput in seconds.
  bool output_computation;	//!< If true ouput computation time.
  
  //  --- CHEKS  ---
  bool check_vols;		//!< If true compute the various input/ouput volumes. 

  //  --- REMOVE  ---
  bool remove_data;		//!< If true remove the data created by the model (no pre-proc).
  bool remove_prec_data;	//!< If true remove the data created by the pre-proces

  // ---  RASTER OUTPUT ---
  bool     rast_vel_as_vect;	//!< If true output the velocity raster as a vector field.  
  CA::Real rast_wd_tol;		//!< If true output the velocity raster as a vector field.  

  // ---  PEAK UPDATE ---
  bool  update_peak_dt;	        //!< If true update the peak at every time step. Default false.

  // --- EXPAND DOMAIN ---
  bool expand_domain;		//!< If true expand the computational domain when needed.
};


//! Initialise the setup structure usign a CSV file. 
//! Each row represents a new "variable" where the 
//! first column is the name of the element 
//! and the following columns have the multiple/single values.
//! \attention The order of elements is not important.
//! \param[in]  filename This is the file where the data is read.
//! \param[out] setup    The structure containing the read data.
//! \return A non zero value if there was an error.
int initSetupFromCSV(const std::string& filename, Setup& setup);


#endif
