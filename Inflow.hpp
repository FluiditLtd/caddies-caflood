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

#ifndef _INFLOW_HPP_
#define _INFLOW_HPP_


//! \file Inflow.hpp
//!  Contains the structure(s) and classe(s) that are used to manage an inflow event.
//! \author Michele Guidolin, University of Exeter, 
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2013-03


#include"ca2D.hpp"
#include"BaseTypes.hpp"
#include"Box.hpp"
#include<string>
#include<vector>

//! Define a small inflow to ignore.
#define SMALL_INFLOW  1.E-10


//! Structure with the configuration value that define an inflow of
//! water event in the CA2D model. The time is in seconds and the
//! inflow is in cubic meters per second. The specific inflow at a
//! specific time is given by linear interpolation between the
//! previous and next values.
struct IEvent
{
  std::string    name;		//!< Name of the event. 

  std::vector<CA::Real> ins;	//!< The list of inflows in cubic meters per second.
  std::vector<CA::Real> times;	//!< The times in seconds.
  std::vector<CA::Real> area;	//!< The area (tl,tr,bl,br) where the inflow will heppen in the border.
  std::vector<CA::Real> zone;	//!< The zone (x,y,w,h) where the inflow will heppen in the border.
  CA::Real              u;	//!< Used to compute the Analytical solution.
  CA::Real              n;	//!< Used to compute the Analytical solution.
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


//! Class that manage all Inflow events
class InflowManager
{
private:
  
  //! The structure used during the model computation to store the inflow
  //! event data.
  struct Data
  {
    size_t   index;	        //!< The index of the inflow data (ins/times).
    CA::Box  box_area;	 	//!< The box of the area where the inflow is set.
    CA::Real grid_area;		//!< Compute the exact grid area, it used for volume checking. 

    CA::Real volume;	        //!< Compute the total volume of inflow of the last update period.

    // These next variable are used to solve the problem of different
    // volume when using float type for Real. The idea is to compute
    // the amount of missing/extra inflow  in comparison to the
    // one expected and add/subtract this `one-off` inflow. 

    double   total_inflow;	//!< The total inflow added during a update period.
    double   expected_inflow;	//!< The expected inflow added during a update period.

    double   one_off_inflow;	//!< The inflow to add/subtract.

    
    Data():
      index(0), box_area(CA::Box::Empty()), grid_area(0.0), volume(0.0),
      total_inflow(0.0), expected_inflow(0.0), one_off_inflow(0.0)
    {}
    
    ~Data()
  {}
  };

public:

  //! Construct a Inflow manager
  InflowManager(CA::Grid&  GRID, const std::vector<IEvent>& ies);

  //! Destroy a Inflow Manager.
  ~InflowManager();

  //! Add the computational domain of the inflow events into the given domain.
  void addDomain(CA::BoxList& compdomain);

  //! Analyse the area where the various inflow event will heppen.
  void analyseArea(CA::CellBuffReal& TMP, CA::CellBuffState& MASK, CA::BoxList&  domain);

  //! Prepare the inflow events for the next update step considering the
  //! simulation time, the lenght of the update step and the next time
  //! step.
  void prepare(CA::Real t, CA::Real period_time_dt, CA::Real next_dt);

  //! Return the volume of inflow of the last period_time_dt. 
  //! \attention This is the PERIOD volume.
  CA::Real volume();

  //! Add the amount of inflow 
  void add(CA::CellBuffReal& WD, CA::CellBuffState& MASK, CA::Real t, CA::Real next_dt);

  //! Compute the potential velocity that could happen in the next
  //! update/period step.
  //! This is used to limit the time step.
  CA::Real potentialVA(CA::Real t, CA::Real period_time_dt); 

  //! Return the simulation time when the events will not add any
  //! further water.
  CA::Real endTime();
  
protected:

  //! Initialise a single inflow event data that is used during the
  //! computation from the inflow event configuration.
  int initData(const IEvent& ie, Data& iedata);

private:

  //! Reference to the grid.
  CA::Grid& _grid;

  //! Reference to the List of inflow events
  const std::vector<IEvent>& _ies;

  //! List of inflow event data.
  std::vector<Data> _datas;

};

#endif
