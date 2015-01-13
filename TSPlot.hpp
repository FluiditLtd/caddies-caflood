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

#ifndef _TSPLOT_HPP_
#define _TSPLOT_HPP_


//! \file TSPlot.hpp
//! Contains the manager that plot the time step(s)
//! \author Michele Guidolin, University of Exeter, 
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2014-07


#include"ca2D.hpp"
#include"BaseTypes.hpp"
#include"ArgsData.hpp"
#include"PointList.hpp"
#include<string>
#include<vector>
#include<fstream>


//! Class that manage all Time Steps Plots outputs
class TSPlot
{
public:

  //! Construct a TimeSteps Plot manager
  //! \param name This is the name of the file.
  TSPlot(std::string name,bool plot);

  //! Destroy a Time Steps Plot Manager.
  ~TSPlot();
  

  //! Output the time plots.
  //! \params t       The simulation time.
  //! \param  dt      The last dt
  void output(CA::Real t, CA::Real dt);

private:

  cpp11::shared_ptr<std::ofstream> _file;   //!< The file where to output the time plot data.

};


#endif
