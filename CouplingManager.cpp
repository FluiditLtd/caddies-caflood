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

#ifdef _WIN32
#include <winsock.h>
#pragma comment(lib, "Ws2_32.lib")
#else
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#endif

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


CouplingManager::CouplingManager(CA::Grid&  GRID, CA::CellBuffReal& ELV, std::vector<ICoupling>& aCoupling, CA::Real aTime_start, CA::Real aTime_end, int aPort):
  grid(GRID),
  coupling(aCoupling),
  time_start(aTime_start),
  time_end(aTime_end),
  port(aPort),
  readValuesUntil(0),
  previousValuesUntil(0),
  networkWaitingUntil(-1),
  sockfd(0),
  points() {

    if (port > 0) {
#ifdef _WIN32
        WSADATA wsaData;
        if (WSAStartup(MAKEWORD(2, 2), &wsaData) != NO_ERROR) {
            std::cerr << "Error in WSAStartup" << std::endl;
            throw;
        }
#endif

        sockfd = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
        sockaddr_in addr;
        addr.sin_family = AF_INET;
        addr.sin_port = htons(port);
        addr.sin_addr.s_addr = inet_addr("127.0.0.1");
        if (connect(sockfd, reinterpret_cast<const sockaddr *>(&addr), sizeof(addr)) != 0) {
            std::cerr << "Cannot connect to the coupling host at 127.0.0.1:" << port << std::endl;
#ifdef _WIN32
        closesocket(sockfd);
#endif
            throw;
        }
    }

    waterDepthBuffer = new CA::Real[coupling.size()];
    createBoxes();
    readElevations(ELV);
}

CouplingManager::~CouplingManager() {
	delete waterDepthBuffer;
}

//! Add the computational domain of the coupled points into the given domain.
void CouplingManager::addDomain(CA::BoxList& compdomain)
{
    for (auto iter = coupling.begin(); iter != coupling.end(); iter++) {
        compdomain.add((*iter).box_area);
    }
}

void CouplingManager::write(std::string line) {
    if (port <= 0) {
        std::cout << line;
        std::cout.flush();
    }
    else {
        const char *data = line.c_str();
        int length = line.length();
        int pos = 0;
        while (length > 0) {
            int bytes = send(sockfd, &data[pos], length, 0);
            if (bytes >= 0) {
                pos += bytes;
                length -= bytes;
            }
            else
                break;
        }
    }
}

std::string CouplingManager::read() {
    if (port <= 0) {
        if (!std::cin.good() || std::cin.eof())
            return "";

        std::string line;
        std::getline(std::cin, line, '\n');
        return line;
    }
    else {
        std::stringstream line("");
        char c;
        int len;
        len = recv(sockfd, &c, 1, 0);
        while (len > 0 && c != '\n') {
            line << c;
            len = recv(sockfd, &c, 1, 0);
        }
        return line.str();
    }
}

void CouplingManager::input(CA::Real time) {
    // Nothing to do as there are enough data
    if (readValuesUntil >= time - 0.0001 || inputEnded || stopped || coupling.empty())
        return;

    // There won't be anything input, if the network simulator
    // is waiting for values after the current time.
    if (networkWaitingUntil >= time - 0.001)
        return;

    // Inform the other end, that we are actually waiting to get some values in
    std::stringstream line("");

    line << "WAITING," << time << "\n";
    write(line.str());

    std::string readLine = "";
    do {
        readLine = read();
        std::stringstream readLineStream(readLine);
        std::vector<std::string> tokens(CA::getLineTokens(readLineStream, ','));
        if (tokens.size() > 0) {
            if (tokens[0] == "END") {
                inputEnded = true;
                break;
            } else if (tokens[0] == "STOP") {
                inputEnded = true;
                stopped = true;
                break;
            } else if (tokens[0] == "WAITING") {
                CA::Real newTime;
                if (!CA::fromString(newTime, tokens[1]))
                    break;

                networkWaitingUntil = newTime;
                if (newTime > readValuesUntil)
                    readValuesUntil = newTime;

                if (newTime >= time - 0.0001)
                    break;
            } else if (tokens[0] == "FLOW") {
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
    } while (readValuesUntil < time - 0.0001 && !inputEnded && !readLine.empty());
}


void CouplingManager::output(CA::Real time, CA::CellBuffReal& WD, CA::CellBuffReal& ELV) {
    if (coupling.empty() || stopped || inputEnded)
        return;

    // No need to output values if nobody is going to use them
    if (time < networkWaitingUntil)
        return;

    std::stringstream line("");
    line << "HEAD," << time;

    // Read the water depth values into the water depth buffer
    WD.retrievePoints(points, waterDepthBuffer, coupling.size());

    int index = 0;
    for (auto iter = coupling.begin(); iter != coupling.end(); iter++, index++) {
        auto& point = *iter;

        // Retrieve the data from the CellBUff into the temporary buffer.
        CA::Real depth = waterDepthBuffer[index];
        line << "," << depth << "," << (depth + point.elv);
    }
    line << "\n";
    write(line.str());
}

void CouplingManager::add(CA::CellBuffReal& WD, CA::CellBuffState& MASK, CA::Real area, CA::Real t, CA::Real dt)
{
  CA::PointList points;
  std::vector<CA::Real> volumes;

  // Loop through the couplings
  for (auto iter = coupling.begin(); iter != coupling.end(); iter++)
  {
    // Compute the inflow volume at specific time using
    // interpolation. Check if the index is the last available
    // one, then there is no inflow.
    CA::Real volume;
    ICoupling &point = *iter;
    if (t < readValuesUntil)
        volume = point.prevFlow * dt;
    else
        volume = point.flow * dt;

    // Add (or subtract) the given volume into the water detph of the
    // given area.
    // Do not add it if it is zero.
    if(std::abs(volume)>=SMALL_INFLOW)      
    {
        points.add(point.box_area.topLeft());
        volumes.push_back(volume);
    }
  }

  unsigned long size = points.size();
  CA::Real *buffer = new CA::Real[size];
  WD.retrievePoints(points, buffer, size);

  for (size_t i = 0; i < size; i++)
      buffer[i] = std::max(static_cast<CA::Real>(0), buffer[i] + volumes[i] / area);

  WD.insertPoints(points, buffer, size);
  delete buffer;
}

void CouplingManager::createBoxes()
{
    std::cout<<"---- COUPLING POINTS ----"<<std::endl;
    for (auto iter = coupling.begin(); iter != coupling.end(); iter++)
    {
        auto& item = *iter;
        CA::Point     tl( CA::Point::create(grid, item.x, item.y) );
        points.add(tl);

        item.box_area = CA::Box(tl.x(), tl.y(), 1, 1);
        //std::cout<<item.name<< " " << item.x << " " << item.y << " box: " << item.box_area.x() << " " << item.box_area.y() << " "  << item.box_area.w() << " " << item.box_area.h() << std::endl;
    }
}

void CouplingManager::readElevations(CA::CellBuffReal& ELV) {
    ELV.retrievePoints(points, waterDepthBuffer, coupling.size());

    int index = 0;
    for (auto iter = coupling.begin(); iter != coupling.end(); iter++, index++) {
        auto &item = *iter;
        item.elv = waterDepthBuffer[index];
    }
}

void CouplingManager::end() {
    if (coupling.size() > 0) {
        write("END\n");
    }
}

void CouplingManager::close() {
#ifdef _WIN32
    closesocket(sockfd);
#else
    ::close(sockfd);
#endif
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
