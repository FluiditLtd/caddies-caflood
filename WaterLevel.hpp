#ifndef _WATERLEVEL_HPP_
#define _WATERLEVEL_HPP_


//! \file WaterLevel.hpp
//!  Contains the structure(s) and classe(s) that are used to manage a WaterLevel event.
//! \author Michele Guidolin, University of Exeter, 
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2013-03


#include"ca2D.hpp"
#include"BaseTypes.hpp"
#include"Box.hpp"
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

//! Initialise the water level event structure usign a CSV file. 
//! Each row represents a new "variable" where the 
//! first column is the name of the element 
//! and the following columns have the multiple/single values.
//! \attention The order of elements is not important.
//! \param[in]  filename This is the file where the data is read.
//! \param[out] setup    The structure containing the read data.
//! \return A non zero value if there was an error.
int initWLEventFromCSV(const std::string& filename, WLEvent& wle);


//! Class that manage all WaterLevel events
class WaterLevelManager
{
private:
  
  //! The structure used during the model computation to store the WaterLevel
  //! event data.
  struct Data
  {
    size_t   index;	        //!< The index of the water level data (wls/times).
    CA::Box  box_area;	 	//!< The box of the area where the water level is set.
    CA::Real grid_area;		//!< Compute the exact grid area, it used for volume checking. 

    CA::Real volume;	        //!< Compute the total volume of the last update period.
    CA::Real last_level;	//!< The last known level set in the cell. This is the elevation at beginning.
    
    Data():
      index(0), box_area(CA::Box::Empty()), grid_area(0.0), volume(0.0), last_level(0.0)
    {}
    
    ~Data()
  {}
  };

public:

  //! Construct a WaterLevel manager
  WaterLevelManager(CA::Grid&  GRID, const std::vector<WLEvent>& wes);

  //! Destroy a WaterLevel Manager.
  ~WaterLevelManager();

  //! Add the computational domain of the WaterLevel events into the given domain.
  void addDomain(CA::BoxList& compdomain);

  //! Analyse the area where the various WaterLevel event will happen.
  void analyseArea(CA::CellBuffReal& TMP, CA::CellBuffState& MASK, CA::BoxList&  domain);

  //! Retrieve the minimum elevation of the given area.
  void getElevation(CA::CellBuffReal& Ele);

  //! Prepare the WaterLevel events for the next update step considering the
  //! simulation time, the lenght of the update step and the next tim
  //! step.
  void prepare(CA::Real t, CA::Real period_time_dt, CA::Real next_dt);

  //! Return the volume of WaterLevel of the last period_time_dt. 
  //! \attention This is the PERIOD volume.
  CA::Real volume();

  //! Add the amount of WaterLevel 
  void add(CA::CellBuffReal& WD, CA::CellBuffReal& ELV, CA::CellBuffState& MASK, CA::Real t, CA::Real next_dt);

  //! Compute the potential velocity that could happen in the next
  //! update/period step.
  //! This is used to limit the time step.
  CA::Real potentialVA(CA::Real t, CA::Real period_time_dt); 

  
protected:

  //! Initialise a single WaterLevel event data that is used during the
  //! computation from the WaterLevel event configuration.
  int initData(const WLEvent& wle, Data& data);

private:

  //! Reference to the grid.
  CA::Grid& _grid;

  //! Reference to the List of WaterLevel events
  const std::vector<WLEvent>& _wles;

  //! List of WaterLevel event data.
  std::vector<Data> _datas;
};


#endif
