/*
 * CouplingManager.hpp
 *
 *  Created on: Jun 15, 2021
 *      Author: sunelma
 */

#ifndef _COUPLINGMANAGER_HPP_
#define _COUPLINGMANAGER_HPP_

#include"ca2D.hpp"
#include"BaseTypes.hpp"
#include"Box.hpp"
#include<string>
#include<vector>

//! Structure with the configuration value that define an inflow of
//! water event in the CA2D model. The time is in seconds and the
//! inflow is in cubic meters per second. The specific inflow at a
//! specific time is given by linear interpolation between the
//! previous and next values.
struct ICoupling
{
public:

  std::string    name;      //!< Coupling component.
  CA::Real       x;         //!< X coordinate of the coupling point
  CA::Real       y;         //!< Y coordinate of the coupling point
  CA::Real       elv;       //!< Elevation of the coupling point
  
  CA::Real       head;      //!< CAFLOOD simulated head
  CA::Real       flow;      //!< Net flow calculated by the hydraulic simulator (+ to surface, - to network)
  CA::Real       prevFlow;  //!< Net flow calculated by the hydraulic simulator (+ to surface, - to network)
  CA::Real       actualFlow;//!< The realized flow calculated by the 2D simulator (+ to surface, - to network)
  CA::Box        box_area;  //!< The box of the area where the flow is set.

  ICoupling():
      box_area(CA::Box::Empty()) {
  }
};



//! Initialise the inflow event structure usign a CSV file. 
//! Each row represents a new "variable" where the 
//! first column is the name of the element 
//! and the following columns have the multiple/single values.
//! \attention The order of elements is not important.
//! \param[in]  filename This is the file where the data is read.
//! \param[out] setup    The structure containing the read data.
//! \return A non zero value if there was an error.
int initICouplingsFromCSV(const std::string& filename, std::vector<ICoupling>& couplings);


class CouplingManager {
private:

  //! Reference to the grid.
  CA::Grid& grid;

  //! Reference to the list of coupling locations events
  std::vector<ICoupling>& coupling;

  //! PointList for the coupling locations for reading the water depths efficiently.
  CA::PointList points;

  //! Reference to the water depth memory for the coupled points
  CA::Real *waterDepthBuffer;

  int port;
  int sockfd;
  CA::Real time_start;
  CA::Real time_end;
  CA::Real readValuesUntil;
  CA::Real previousValuesUntil;
  CA::Real networkWaitingUntil;
  double coupledVolume = 0.0;
  bool inputEnded = false;
  bool stopped = false;

public:
    CouplingManager(CA::Grid&  GRID, CA::CellBuffReal& ELV, std::vector<ICoupling>& aCoupling, CA::Real time_start, CA::Real time_end, int port);
    ~CouplingManager();

    void addDomain(CA::BoxList& compdomain);
    inline bool isStopped() { return stopped; }
    void input(CA::Real t);
    void output(CA::Real time, CA::CellBuffReal& WD, CA::CellBuffReal& ELV);
    void add(CA::CellBuffReal& WD, CA::CellBuffState& MASK, CA::Real area, CA::Real t, CA::Real dt);
    void end();
    void close();

    CA::Real potentialVA(CA::Real t, CA::Real period_time_dt);

    CA::Real volume() {
        return getAndResetCoupledVolume();
    }

    CA::Real endTime();

    CA::Real getAndResetCoupledVolume() {
        CA::Real ret = (CA::Real)coupledVolume;
        coupledVolume = 0.0;
        return ret;
    }

private:
    void createBoxes();
    void readElevations(CA::CellBuffReal& ELV);
    void write(std::string line);
    std::string read();
};

#endif
