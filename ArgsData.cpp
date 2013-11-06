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
