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
#include<string>
#include<vector>


//! Structure with the configuration value that define a rain event in
//! the CA2D model. A rain event is the amount of rain intensity
//! falling over time. The time is represented the end onf the
//! intensity from t=0 of the rain in second while the rain intensity
//! value is represented in mm/hr.
struct RainEvent
{
  std::string    name;		//!< Name of the event. 

  std::vector<CA::Real> rains;	//!< The list of rain intensities in mm/hr 
  std::vector<CA::Real> times;	//!< The times when the rain intensities stop in seconds.
  std::vector<CA::Real> area;	//!< The area where the rain will fall (if empty all domain).  
  std::vector<CA::Real> zone;	//!< The zone (x,y,w,h) where the rain will heppen (if empty all domain).
};


//! The structure used during the model computation to store the rain
//! event data.
struct REData
{
  size_t   index;	        //!< The index of the rain data (rains/times).
  CA::Box  box_area;	 	//!< The box of the area where the rain fall.
  CA::Real grid_area;		//!< Compute the exact grid area, it used for volume checking. 
  CA::Real volume;		//!< Compute the total volume of water added by rain. 

  REData():
    index(0), box_area(CA::Box::Empty())
  {}

  ~REData()
  {}
};


//! Initialise the rain event structure usign a CSV file. 
//! Each row represents a new "variable" where the 
//! first column is the name of the element 
//! and the following columns have the multiple/single values.
//! \attention The order of elements is not important.
//! \param[in]  filename This is the file where the data is read.
//! \param[out] setup    The structure containing the read data.
//! \return A non zero value if there was an error.
int initRainEventFromCSV(const std::string& filename, RainEvent& re);


//! Initialise the rain event data that is used during the
//! computation from the rain event configuration.
//! \return A non zero value if there was an error.
int initREData(const CA::Grid&  GRID, const RainEvent& re, REData& redata);


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


//! Structure with the configuration value that define an inflow of
//! water event in the CA2D model. The time is in seconds and the
//! inflwo is in cubic meters per second. The specific inflow at a
//! specific time is given by linear interpolation between the
//! previous and next values.
struct IEvent
{
  std::string    name;		//!< Name of the event. 

  std::vector<CA::Real> ins;	//!< The list of inflows in cubic meters per second.
  std::vector<CA::Real> times;	//!< The times in seconds.
  std::vector<CA::Real> area;	//!< The area (tl,tr,bl,br) where the inflow will heppen in the border.
  std::vector<CA::Real> zone;	//!< The zone (x,y,w,h) where the inflow will heppen in the border.
};


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


//! Initialise the inflow event structure usign a CSV file. 
//! Each row represents a new "variable" where the 
//! first column is the name of the element 
//! and the following columns have the multiple/single values.
//! \attention The order of elements is not important.
//! \param[in]  filename This is the file where the data is read.
//! \param[out] setup    The structure containing the read data.
//! \return A non zero value if there was an error.
int initIEventFromCSV(const std::string& filename, IEvent& ie);


//! Initialise the inflow  event data that is used during the
//! computation from the rain event configuration.
//! \return A non zero value if there was an error.
int initIEData(const CA::Grid&  GRID, const IEvent& ie, IEData& iedata);



#endif
