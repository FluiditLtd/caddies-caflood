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

//! \file TSPlot.cpp
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2014-07


#include"ca2D.hpp"
#include"ArgsData.hpp"
#include"TSPlot.hpp"
#include"Utilities.hpp"
#include<iostream>
#include<fstream>


TSPlot::TSPlot(std::string name, bool plot):
_file()
{
  if(plot)
  {
    // Create file
    _file.reset( new std::ofstream(name.c_str()) );
    
    if(_file->good())
    {
      
      // Set  manipulators 
      _file->setf(std::ios::fixed, std::ios::floatfield);
      _file->precision(6); 
      
      // Write the header
      (*_file)<<"t (s),  dt (s)"<<std::endl;
    }
  }
}
  

TSPlot::~TSPlot()
{

}


void TSPlot::output(CA::Real t, CA::Real dt)
{
  if(_file && _file->good())
  {
    // Write line
    (*_file)<<t<<", "<<dt<<std::endl;
  }
}
