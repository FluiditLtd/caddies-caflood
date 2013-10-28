//! \file ArgsData.cpp
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2013-03


#include"ArgsData.hpp"
#include"Utilities.hpp"
#include<iostream>
#include<fstream>


// Transform an input string into a physical variable enum.
std::istream& operator>>(std::istream& in, PV::Type& pv)
{
  CA::String tmp;
  in>>tmp;

  if( CA::compareCaseInsensitive(tmp, "WD" ) )
  {
    pv = PV::WD; 
    return in;
  }
  if( CA::compareCaseInsensitive(tmp, "WL" ) )
  {
    pv = PV::WL; 
    return in;
  }
  if( CA::compareCaseInsensitive(tmp, "VEL" ) )
  {
    pv = PV::VEL; 
    return in;
  }
  pv = PV::UNKNOWN;
  in.setstate(std::ios::failbit);
  return in;
}


// Transform a physical variable enum into an output string.
std::ostream& operator<<(std::ostream& out, PV::Type& pv)
{
  switch(pv)
  {
  case PV::WD :
    out<<"WD";
    break;
  case PV::WL :
    out<<"WL";
    break;
  case PV::VEL :
    out<<"VEL";
    break;
  default:
    out<<"UNKNOWN";
  }

  return out;
}


int initVELBuff(CA::Grid& GRID, double nodata, VELBuff& vb)
{
  vb.init = true;
  vb.V.reset( new CA::CellBuffReal(GRID) ); 
  vb.U.reset( new CA::CellBuffReal(GRID) ); 
  vb.V->fill(GRID.box(),nodata);
  vb.U->fill(GRID.box(),nodata);
  return 0;
}
