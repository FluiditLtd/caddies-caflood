/*
 * CouplingManager.cpp
 *
 *  Created on: Jun 15, 2021
 *      Author: sunelma
 */

#include "ca2D.hpp"
#include "ArgsData.hpp"
#include "CouplingManager.hpp"
#include "Utilities.hpp"
#include <iostream>
#include <fstream>
#include CA_2D_INCLUDE(addInflow)

//! Define a small inflow to ignore.
#define SMALL_INFLOW  1.E-10


int initICouplingsFromCSV(const std::string& filename, std::vector<ICoupling>& couplings)
{
    std::ifstream ifile(filename.c_str());
    
    if(!ifile)
    {
        std::cerr<<"Error opening couplings CSV file: "<<filename<<std::endl;
        return 1;
    }
    
    // Parse the file line by line until the end of file 
    // and retrieve the tokens of each line.
    while(!ifile.eof())
    {   
      bool found_tok = false;
        std::vector<std::string> tokens( CA::getLineTokens(ifile, ',') );

        // If the tokens vector is empty we reached the eof or an
        // empty line... continue.
        if(tokens.empty())
            continue;

        if (tokens.size() < 3)
        {
            std::cerr<<"Not enough tokens for a coupling"<<std::endl;
            return 1;
        }

        ICoupling coupling;
        {
            std::string str;
            READ_TOKEN(found_tok,str,tokens[0],tokens[0]);

            coupling.name = CA::trimToken(str," \t\r");
        }

        {
            CA::Real value;
            READ_TOKEN(found_tok, value, tokens[1], tokens[0]);
            coupling.x = value;
            READ_TOKEN(found_tok, value, tokens[2], tokens[0]);
            coupling.y = value;
        }

        coupling.head = 0;
        coupling.flow = 0;
        couplings.push_back(coupling);
    }
    
    return 0;
}


CouplingManager::CouplingManager(CA::Grid&  GRID, std::vector<ICoupling>& aCoupling, CA::Real aTime_start, CA::Real aTime_end):
  grid(GRID),
  coupling(aCoupling),
  time_start(aTime_start),
  time_end(aTime_end),
  readValuesUntil(0),
  previousValuesUntil(0),
  networkWaitingUntil(-1) {
}

CouplingManager::~CouplingManager() {
	
}


void CouplingManager::input(CA::Real time) {
    // Nothing to do as there are enough data
    if (readValuesUntil >= time || inputEnded || stopped || coupling.empty())
        return;

    // There won't be anything input, if the network simulator
    // is waiting for values after the current time.
    if (networkWaitingUntil >= time)
        return;

    // Inform the other end, that we are actually waiting to get some values in
    std::cout << "WAITING," << time << std::endl;
    std::cout.flush();

    while (readValuesUntil < time && !inputEnded && std::cin.good() && !std::cin.eof()) {
        std::vector<std::string> tokens( CA::getLineTokens(std::cin, ',') );
        if (tokens.size() > 0) {
            if (tokens[0] == "END") {
                inputEnded = true;
                break;
            }
            else if (tokens[0] == "STOP") {
                inputEnded = true;
                stopped = true;
                break;
            }
            else if (tokens[0] == "WAITING") {
                CA::Real newTime;
                if (!CA::fromString(newTime, tokens[1]))
                    break;

                networkWaitingUntil = newTime;
                if (newTime >= time)
                    break;
            }
            else if (tokens[0] == "FLOW") {
                CA::Real newTime;
                previousValuesUntil = readValuesUntil;
                if (!CA::fromString(newTime, tokens[1]))
                        continue;

                readValuesUntil = newTime;

                for (int index = 0; index < tokens.size() - 2 && index < coupling.size(); index++) {
                    CA::Real flow;
                    if (!CA::fromString(flow, tokens[index + 2]))
                        continue;

                    coupling[index].prevFlow = coupling[index].flow;
                    coupling[index].flow = flow;
                }
            }
        }
    }
}


void CouplingManager::output(CA::Real time, CA::CellBuffReal& WD, CA::CellBuffReal& ELV) {
    if (coupling.empty() || stopped || inputEnded)
        return;

    // No need to output values if nobody is going to use them
    if (time < networkWaitingUntil)
        return;

    std::cout << "HEAD," << time;
  
    for (auto iter = coupling.begin(); iter != coupling.end(); iter++) {
        auto& point = *iter;
        CA::Real depth;
        CA::Real elevation;

        // Retrieve the data from the CellBUff into the temporary buffer.
        WD.retrieveData(point.box_area, &depth, 1, 1);
        ELV.retrieveData(point.box_area, &elevation, 1, 1);
        std::cout << "," << depth << "," << (depth + elevation);
    }
    std::cout << std::endl;
    std::cout.flush();
}

void CouplingManager::add(CA::CellBuffReal& WD, CA::CellBuffState& MASK, CA::Real t, CA::Real dt)
{
  // Loop through the couplings
  for(size_t i = 0; i < coupling.size(); ++i)
  {
    // Compute the inflow volume at specific time using
    // interpolation. Check if the index is the last available
    // one, then there is no inflow.
    CA::Real volume;
    if (t < readValuesUntil)
        volume = coupling[i].prevFlow * dt;
    else
        volume = coupling[i].flow * dt;

    // Add (or subtract) the given volume into the water detph of the
    // given area.
    // Do not add it if it is zero.
    if(std::abs(volume)>=SMALL_INFLOW)      
    {
      CA::Execute::function(coupling[i].box_area, addInflow, grid, WD, MASK, volume);           
    }
  }
}

void CouplingManager::createBoxes()
{
    std::cout<<"---- COUPLING POINTS ----"<<std::endl;
    for (auto iter = coupling.begin(); iter != coupling.end(); iter++)
    {
        auto& item = *iter;
        CA::Point     tl( CA::Point::create(grid, item.x, item.y) );
        item.box_area = CA::Box(tl.x(), tl.y(), 1, 1);
        std::cout<<item.name<< " " << item.x << " " << item.y << " box: " << item.box_area.x() << " " << item.box_area.y() << " "  << item.box_area.w() << " " << item.box_area.h() << std::endl;
    }
}

void CouplingManager::end() {
    if (coupling.size() > 0 && !inputEnded && !stopped) {
        std::cout << "END" << std::endl;
        std::cout.flush();
    }
}

CA::Real CouplingManager::potentialVA(CA::Real t, CA::Real period_time_dt)
{
    CA::Real potential_va = 0.0;

    for (auto iter = coupling.begin(); iter != coupling.end(); iter++)
    {
        auto& item = *iter;
        CA::Real volume = item.flow * period_time_dt;
        CA::Real wd = volume / grid.length() / grid.length();

        // Compute the potential velocity.
        potential_va = std::max(potential_va, std::sqrt(wd * static_cast<CA::Real>(9.81)));
    }

    // ATTENTION! This does not need to be precise but just give a rough estimation
    return potential_va;
}


CA::Real CouplingManager::endTime()
{
    if (coupling.size() == 0)
        return time_start;
    else
        return time_end;
}