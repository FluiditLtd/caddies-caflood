#ifndef _EVENTS_HPP_
#define _EVENTS_HPP_


//! \file Events.hpp
//!  Contains the structure(s) of the various event(s) of the CA2D model.
//! \author Michele Guidolin, University of Exeter, 
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2013-03


#include"ca2D.hpp"
#include"BaseTypes.hpp"
#include"Box.hpp"
#include"Inflow.hpp"
#include<string>
#include<vector>


//! Structure with the configuration value that define a water level
//! event in the CA2D model. A water level event is the level of water
//! in a specific area at specific time. The time is in seconds and
//! the level is in meters. The specific level at a specific time is
//! given by linear interpolation between the previous and next
//! values.
struct WLEvent
{
  std::string    name;		//!< Name of the event. 

  std::vector<CA::Real> wls;	//!< The list of water level in meters.
  std::vector<CA::Real> times;	//!< The times in seconds.
  std::vector<CA::Real> area;	//!< The area where the water level will be set.  
  std::vector<CA::Real> zone;	//!< The zone (x,y,w,h) where the water level will heppen.
};


//! The structure used during the model computation to store the water
//! level event data.
struct WLEData
{
  size_t   index;	        //!< The index of the water level data (wls/times).
  CA::Box  box_area;	 	//!< The box of the area where the water level is set.
  CA::Real min_elv;		//!< The minimum elevation of the area.

  WLEData():
    index(0), box_area(CA::Box::Empty()), min_elv(0)
  {}

  ~WLEData()
  {}
};

//! Initialise the water level event structure usign a CSV file. 
//! Each row represents a new "variable" where the 
//! first column is the name of the element 
//! and the following columns have the multiple/single values.
//! \attention The order of elements is not important.
//! \param[in]  filename This is the file where the data is read.
//! \param[out] setup    The structure containing the read data.
//! \return A non zero value if there was an error.
int initWLEventFromCSV(const std::string& filename, WLEvent& wle);


//! Initialise the water level event data that is used during the
//! computation from the rain event configuration.
//! \return A non zero value if there was an error.
int initWLEData(const CA::Grid&  GRID, const WLEvent& wle, WLEData& wledata);



//! The structure used during the model computation to store the inflow
//! event data.
struct IEData
{
  size_t   index;	        //!< The index of the inflow data (ins/times).
  CA::Box  box_area;	 	//!< The box of the area where the inflow is set.
  CA::Real grid_area;		//!< Compute the exact grid area, it used for volume checking. 
  CA::Real volume;		//!< Compute the total volume of water added by rain. 


  IEData():
    index(0), box_area(CA::Box::Empty())
  {}

  ~IEData()
  {}
};


//! Initialise the inflow  event data that is used during the
//! computation from the rain event configuration.
//! \return A non zero value if there was an error.
int initIEData(const CA::Grid&  GRID, const IEvent& ie, IEData& iedata);



#endif
