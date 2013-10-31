#ifndef _RAIN_HPP_
#define _RAIN_HPP_


//! \file Rain.hpp
//!  Contains the structure(s) and classe(s) that are used to manage a rain event.
//! \author Michele Guidolin, University of Exeter, 
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2013-03


#include"ca2D.hpp"
#include"BaseTypes.hpp"
#include"Box.hpp"
#include<string>
#include<vector>

//! Define a small rain to ignore.
#define SMALL_RAIN  1.E-10

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

//! Initialise the rain event structure usign a CSV file. 
//! Each row represents a new "variable" where the 
//! first column is the name of the element 
//! and the following columns have the multiple/single values.
//! \attention The order of elements is not important.
//! \param[in]  filename This is the file where the data is read.
//! \param[out] setup    The structure containing the read data.
//! \return A non zero value if there was an error.
int initRainEventFromCSV(const std::string& filename, RainEvent& re);


//! Class that manage all Rain events
class RainManager
{
private:
  
  //! The structure used during the model computation to store the rain
  //! event data.
  struct Data
  {
    size_t   index;	        //!< The index of the rain data (rains/times).
    CA::Box  box_area;	 	//!< The box of the area where the rain fall.
    CA::Real grid_area;		//!< Compute the exact grid area, it used for volume checking. 
    CA::Real volume;		//!< Compute the total volume of water added by rain. 

    CA::Real rain;		//!< The amount of rain for each dt of the next period. 
    
    Data():
      index(0), box_area(CA::Box::Empty())
    {}
    
    ~Data()
  {}
  };

public:

  //! Construct a Rain manager
  RainManager(CA::Grid&  GRID, const std::vector<RainEvent>& res);

  //! Destroy a Rain Manager.
  ~RainManager();

  //! Add the computational domain of the rain events into the given domain.
  void addDomain(CA::BoxList& compdomain);

  //! Analyse the area where the various rain event will fall.
  void analyseArea(CA::CellBuffReal& TMP, CA::CellBuffState& MASK, CA::BoxList&  domain);

  //! Prepare the rain events for the next update step considering the
  //! simulation time, the lenght of the update step and the next tim
  //! step.
  void prepare(CA::Real t, CA::Real period_time_dt, CA::Real next_dt);

  //! Add the amount of rain that was previously prepared into the
  //! water depth
  void add(CA::CellBuffReal& WD, CA::CellBuffState& MASK);

  //! Compute the potential velocity that could happen in the next
  //! update/period step.
  //! This is used to limit the time step.
  CA::Real potentialVA(CA::Real t, CA::Real period_time_dt); 

protected:

  //! Initialise a single rain event data that is used during the
  //! computation from the rain event configuration.
  int initData(const RainEvent& re, Data& redata);

private:

  //! Reference to the grid.
  CA::Grid& _grid;

  //! Reference to the List of rain events
  const std::vector<RainEvent>& _res;

  //! List of rain event data.
  std::vector<Data> _datas;
};

#endif
