//! \file main.cpp
//! Program/demo that perform a generic (selected) flood modelling 
//! The ca algorithm is developed in collaboration with BIDUR GHIMIRE & ALBERT CHEN   
//! contact: m.guidolin [at] exeter.ac.uk
//! \date 2013-03


#include"ca2D.hpp"
#include"Arguments.hpp"
#include"Options.hpp"
#include"ArgsData.hpp"
#include"Setup.hpp"
#include"Rain.hpp"
#include"Events.hpp"
#include"TimePlot.hpp"
#include"RasterGrid.hpp"


//! Perform the pre-processing using the given parameters (see pre_proc.cpp).
//! \param[in] ad       The arguments data.
//! \param[in] setup    The setup of the simulation.
//! \param[in] ele_file The elevation file.
//! \return A non zero value if there was an error.
int preProc(const ArgsData& ad, const Setup& setup, const std::string& ele_file); 


//! Perform the post processing of the data for a CA 2D model. 
//! \param[in] ad       The arguments data.
//! \param[in] setup    The setup of the simulation.
//! \param[in] eg       The elevation grid.
//! \param[in] tps      The list of time plot outputs.
//! \param[in] rgs      The list of raster grid outputs.
//! \return A non zero value if there was an error.
int postProc(const ArgsData& ad, const Setup& setup, CA::AsciiGrid<CA::Real>& eg,
	     const std::vector<TimePlot>& tps, const std::vector<RasterGrid>& rgs);


//! Perform the WCA2D flood modelling algorithm using the given
//! parameters (see WCA2D.cpp).
//! \param[in] ad       The arguments data.
//! \param[in] setup    The setup of the simulation.
//! \param[in] eg       The elevation grid.
//! \param[in] res      The list of rain event inputs.
//! \param[in] wles     The list of water level event inputs.
//! \param[in] ies      The list of inflow event inputs.
//! \param[in] tps      The list of time plot outputs.
//! \param[in] rgs      The list of raster grid outputs.
//! \return A non zero value if there was an error.
int WCA2D(const ArgsData& ad, const Setup& setup, const CA::AsciiGrid<CA::Real>& eg, 
	  const std::vector<RainEvent>& res, const std::vector<WLEvent>& wles, const std::vector<IEvent>& ies, 
	  const std::vector<TimePlot>& tps, const std::vector<RasterGrid>& rgs);


int main(int argc, char* argv[])
{
  // Initialise the 2D caAPI.
  CA::init2D(&argc,&argv);

  // Create the argument structure.
  ArgsData ad;

  // The number of argument.
  CA::Unsigned na = 0;

  // There is only three compulsory arguments which are the directory
  // where the configuration files and data is found, the initial
  // setup file of the given CA algorithm, and the directory where the
  // output data is saved. The setup file should be in the working
  // directory.
  ad.args.add(na++,"working-dir", "The working dir which contain the setup file and the data","",  false);
  ad.args.add(na++,"setup-file",  "The setup file of the ca algorithm","",  false);
  ad.args.add(na++,"output-dir",  "The output dir where output files are saved","",  false);

  // Add optional arguments.
  ad.args.add(na++,"help",      "Display the help and exit.", "",  true, false, true);
  ad.args.add(na++,"info",      "Print information about the configuration on console",  "",  true, false);
  ad.args.add(na++,"pre-proc",  "Perform the pre-processing of the data (only)",         "",  true, false);
  ad.args.add(na++,"WCA2D",     "Perform the WCA2D flood model",                        "",  true, false);

  // Add the options from the CA implementation
  ad.args.addList(CA::options());

  // Parse the argument list, if the parse fail (it returns false) then exit!
  if(ad.args.parse(argc,argv,std::cout) == false )
  {
#if defined _WIN32 || defined __CYGWIN__   
    std::cout<<"Use the /help options to show help"<<std::endl;
    std::cout<<"Press 'Return' to continue"<<std::endl;
    std::cin.get();
#else
    std::cout<<"Use the -help options to show help"<<std::endl;
    std::cout<<"Press 'Return' to continue"<<std::endl;
    std::cin.get();
#endif
    return EXIT_FAILURE;    
  }
  
  // Cycle throught the arguments that were activated
  for(CA::Arguments::List::const_iterator i = ad.args.active().begin(); i != ad.args.active().end();  ++i)
  {
    if((*i)->tag == 0)		// Working dir.
    {
      ad.working_dir = (*i)->value;
    }

    if((*i)->tag == 1)		// setup file.
    {
      ad.setup_file = (*i)->value;
    }

    if((*i)->tag == 2)		// Ouptut dir.
    {
      ad.output_dir = (*i)->value;
    }

    if((*i)->name == "help")
    {
      ad.args.help(std::cout,true);
      return -1;
    }
    
    if((*i)->name == "info")
      ad.info = true;

    if((*i)->name == "pre-proc")
      ad.pre_proc = true;

    if((*i)->name == "WCA2D")
    {
      ad.pre_proc = true;
      ad.WCA2D = true;
      ad.post_proc = true;
    }

  }

  // Set the data directory.
  ad.data_dir = ad.output_dir + ad.sdir;
  
  if(ad.info)
  {
    std::cout<<"caflood 1.0"<<std::endl;
    std::cout<<"CA API Version       : "<<caVersion<<std::endl;
    std::cout<<"       Impl Name     : "<<caImplName<<std::endl;
    std::cout<<"       Impl Short    : "<<caImplShortName<<std::endl;
    std::cout<<"       Impl Version  : "<<caImplVersion<<std::endl;
    std::cout<<"       Impl Precision: "<<caImplPrecision<<std::endl;
    std::cout<<"Working Dir          : "<<ad.working_dir<<std::endl;
    std::cout<<"Setup File           : "<<ad.setup_file<<std::endl;
    std::cout<<"Output Dir           : "<<ad.output_dir<<std::endl;
    std::cout<<"Data   Dir           : "<<ad.data_dir<<std::endl;
  }

  if(ad.info)
    std::cout<<std::endl<<"Load simulation setup"<<std::endl;

  // Create the setup structure.
  Setup setup;
  
  // Load setup information from the CA2D.csv file in the working dir.  
  std::string setup_file = ad.working_dir+ad.sdir+ad.setup_file;

  if(initSetupFromCSV(setup_file, setup)!=0)
  {
    std::cerr<<"Error reading Setup CSV file: "<<setup_file<<std::endl;
    return EXIT_FAILURE;    
  }
  
  if(ad.info)
  {
    std::cout<<"Simulation              : "<<setup.sim_name<<std::endl;
    std::cout<<"Short Name              : "<<setup.short_name<<std::endl;
    std::cout<<"Version                 : ";
    for(size_t i = 0; i< setup.sim_version.size(); ++i)
      std::cout<<setup.sim_version[i]<<" ";
    std::cout<<std::endl;
    std::cout<<"Time Start              : "<<setup.time_start<<std::endl;
    std::cout<<"Time End                : "<<setup.time_end<<std::endl;
    std::cout<<"Max DT                  : "<<setup.time_maxdt<<std::endl;
    std::cout<<"Min DT                  : "<<setup.time_mindt<<std::endl;
    std::cout<<"Update DT               : "<<setup.time_updatedt<<std::endl;
    std::cout<<"Alpha                   : "<<setup.time_alpha<<std::endl;
    std::cout<<"Max Iterations          : "<<setup.time_maxiters<<std::endl;
    std::cout<<"Roughness Global        : "<<setup.roughness_global<<std::endl;
    std::cout<<"Ignore WD               : "<<setup.ignore_wd<<std::endl;
    std::cout<<"Tolerance               : "<<setup.tolerance<<std::endl;
    std::cout<<"Boundary Ele            : "<<setup.boundary_elv<<std::endl;
    std::cout<<"Elevation ASCII         : "<<setup.elevation_ASCII<<std::endl;
    std::cout<<"Rain Event CSV          : ";
    for(size_t i = 0; i< setup.rainevent_files.size(); ++i)
      std::cout<<setup.rainevent_files[i]<<" ";
    std::cout<<std::endl;
    std::cout<<"Water Level Event CSV   : ";
    for(size_t i = 0; i< setup.wlevent_files.size(); ++i)
      std::cout<<setup.wlevent_files[i]<<" ";
    std::cout<<std::endl;
    std::cout<<"Inflow Event CSV        : ";
    for(size_t i = 0; i< setup.wlevent_files.size(); ++i)
      std::cout<<setup.wlevent_files[i]<<" ";
    std::cout<<std::endl;
    std::cout<<"Time Plot CSV           : ";
    for(size_t i = 0; i< setup.timeplot_files.size(); ++i)
      std::cout<<setup.timeplot_files[i]<<" ";
    std::cout<<std::endl;
    std::cout<<"Raster CSV              : ";
    for(size_t i = 0; i< setup.rastergrid_files.size(); ++i)
      std::cout<<setup.rastergrid_files[i]<<" ";
    std::cout<<std::endl;
    std::cout<<"Output Console          : "<<setup.output_console<<std::endl;
    std::cout<<"Output Period           : "<<setup.output_period<<std::endl;
    std::cout<<"Output Computation Time : "<<setup.output_computation<<std::endl;
    std::cout<<"Check Volumes           : "<<setup.check_vols<<std::endl;
    std::cout<<"Remove Proc Data        : "<<setup.remove_data<<std::endl;
    std::cout<<"Remove Pre-Proc Data    : "<<setup.remove_prec_data<<std::endl;
    std::cout<<"Raster VEL Vector       : "<<setup.rast_vel_as_vect<<std::endl;
    std::cout<<"Raster WD Tol           : "<<setup.rast_wd_tol<<std::endl;
    std::cout<<"Update Peak Every DT    : "<<setup.update_peak_dt<<std::endl;
    std::cout<<"Expand Domain           : "<<setup.expand_domain<<std::endl;
  }

  if(ad.info)
    std::cout<<std::endl<<"Load elevation data "<<std::endl;

  //! Load the eventual elevation.
  std::string ele_file = ad.working_dir+ad.sdir+setup.elevation_ASCII; 
  CA::AsciiGrid<CA::Real> eg;

  // ATTENTION. Load only the header here .. not the actual data
  CA::readAsciiGridHeader(eg,ele_file);

  if(ad.info)
  {
    std::cout<<"Elevation ASCII    : "<<setup.elevation_ASCII<<std::endl;
    std::cout<<"ncols              : "<<eg.ncols<<std::endl;
    std::cout<<"nrows              : "<<eg.nrows<<std::endl;
    std::cout<<"xllcorner          : "<<eg.xllcorner<<std::endl;
    std::cout<<"yllcorner          : "<<eg.yllcorner<<std::endl;
    std::cout<<"cellsize           : "<<eg.cellsize<<std::endl;
    std::cout<<"nodata             : "<<eg.nodata<<std::endl;
  }

  if(ad.info)
    std::cout<<std::endl<<"Load rain event inputs configuration "<<std::endl;

  // Load any eventual rain event.
  std::vector<RainEvent> res; 
  for(size_t i = 0; i< setup.rainevent_files.size(); ++i)
  {
    std::string file = ad.working_dir+ad.sdir+setup.rainevent_files[i]; 

    RainEvent re;
    if(initRainEventFromCSV(file, re)!=0)
    {
      std::cerr<<"Error reading Rain Event CSV file: "<<file<<std::endl;
      return EXIT_FAILURE;    
    }

    if(ad.info)
    {
      std::cout<<"Rain Event         : "<<re.name<<std::endl;
      std::cout<<"Rain Intensity     : ";
      for(size_t i = 0; i< re.rains.size(); ++i)
	std::cout<<re.rains[i]<<" ";
      std::cout<<std::endl;
      std::cout<<"Time Stop          : ";
      for(size_t i = 0; i< re.times.size(); ++i)
	std::cout<<re.times[i]<<" ";
      std::cout<<std::endl;
      std::cout<<"Area               : ";
      for(size_t i = 0; i< re.area.size(); ++i)
	std::cout<<re.area[i]<<" ";
      std::cout<<std::endl;
    }

    res.push_back(re);    
  } 

  if(ad.info)
    std::cout<<std::endl<<"Load water level event inputs configuration "<<std::endl;

  // Load any eventual water level events
  std::vector<WLEvent> wles; 
  for(size_t i = 0; i< setup.wlevent_files.size(); ++i)
  {
    std::string file = ad.working_dir+ad.sdir+setup.wlevent_files[i]; 

    WLEvent wle;
    if(initWLEventFromCSV(file, wle)!=0)
    {
      std::cerr<<"Error reading Water Level Event CSV file: "<<file<<std::endl;
      return EXIT_FAILURE;    
    }

    if(ad.info)
    {
      std::cout<<"Water Level Event  : "<<wle.name<<std::endl;
      std::cout<<"Water Level        : ";
      for(size_t i = 0; i< wle.wls.size(); ++i)
	std::cout<<wle.wls[i]<<" ";
      std::cout<<std::endl;
      std::cout<<"Time               : ";
      for(size_t i = 0; i< wle.times.size(); ++i)
	std::cout<<wle.times[i]<<" ";
      std::cout<<std::endl;
      std::cout<<"Area               : ";
      for(size_t i = 0; i< wle.area.size(); ++i)
	std::cout<<wle.area[i]<<" ";
      std::cout<<std::endl;
    }

    wles.push_back(wle);    
  } 

  if(ad.info)
    std::cout<<std::endl<<"Load inflow event inputs configuration "<<std::endl;

  // Load any eventual inflow events
  std::vector<IEvent> ies; 
  for(size_t i = 0; i< setup.inflowevent_files.size(); ++i)
  {
    std::string file = ad.working_dir+ad.sdir+setup.inflowevent_files[i]; 

    IEvent ie;

    if(initIEventFromCSV(file, ie)!=0)
    {
      std::cerr<<"Error reading Inflow Event CSV file: "<<file<<std::endl;
      return EXIT_FAILURE;    
    }

    if(ad.info)
    {
      std::cout<<"Inflow Event       : "<<ie.name<<std::endl;
      std::cout<<"Inflow (cumecs)    : ";
      for(size_t i = 0; i< ie.ins.size(); ++i)
	std::cout<<ie.ins[i]<<" ";
      std::cout<<std::endl;
      std::cout<<"Time               : ";
      for(size_t i = 0; i< ie.times.size(); ++i)
	std::cout<<ie.times[i]<<" ";
      std::cout<<std::endl;
      std::cout<<"Area               : ";
      for(size_t i = 0; i< ie.area.size(); ++i)
	std::cout<<ie.area[i]<<" ";
      std::cout<<std::endl;
    }

    ies.push_back(ie);    
  } 

  if(ad.info)
    std::cout<<std::endl<<"Load time plot outputs configuration "<<std::endl;

  // Load any eventual time plots.
  std::vector<TimePlot> tps; 
  for(size_t i = 0; i< setup.timeplot_files.size(); ++i)
  {
    std::string file = ad.working_dir+ad.sdir+setup.timeplot_files[i]; 

    TimePlot tp;
    if(initTimePlotFromCSV(file, tp)!=0)
    {
      std::cerr<<"Error reading Time Plot CSV file: "<<file<<std::endl;
      return EXIT_FAILURE;    
    }

    if(ad.info)
    {
      std::cout<<"Time Plot Name     : "<<tp.name<<std::endl;
      std::cout<<"Physical Variable  : "<<tp.pv<<std::endl;
      std::cout<<"Points X Name      : ";
      for(size_t i = 0; i< tp.xcoos.size(); ++i)
	std::cout<<tp.pnames[i]<<" ";
      std::cout<<std::endl;
      std::cout<<"Points X Coo       : ";
      for(size_t i = 0; i< tp.xcoos.size(); ++i)
	std::cout<<tp.xcoos[i]<<" ";
      std::cout<<std::endl;
      std::cout<<"Points Y Coo       : ";
      for(size_t i = 0; i< tp.ycoos.size(); ++i)
	std::cout<<tp.ycoos[i]<<" ";
      std::cout<<std::endl;  
      std::cout<<"Period             : "<<tp.period<<std::endl;
    }

    tps.push_back(tp);    
  } 

  if(ad.info)
    std::cout<<std::endl<<"Load raster grid outputs configuration "<<std::endl;

  // Load any eventual raster grids.
  std::vector<RasterGrid> rgs; 
  for(size_t i = 0; i< setup.rastergrid_files.size(); ++i)
  {
    std::string file = ad.working_dir+ad.sdir+setup.rastergrid_files[i]; 

    RasterGrid rg;
    if(initRasterGridFromCSV(file, rg)!=0)
    {
      std::cerr<<"Error reading Time Plot CSV file: "<<file<<std::endl;
      return EXIT_FAILURE;    
    }

    if(ad.info)
    {
      std::cout<<"Raster Grid Name   : "<<rg.name<<std::endl;
      std::cout<<"Physical Variable  : "<<rg.pv<<std::endl;
      std::cout<<"Peak               : "<<rg.peak<<std::endl;
      std::cout<<"Period             : "<<rg.period<<std::endl;
    }

    rgs.push_back(rg);    
  } 

  //! Now perform the pre-processing
  if(ad.pre_proc)
  {
    if(ad.info)
      std::cout<<std::endl<<"Starting pre-processing data "<<std::endl;

    if(preProc(ad,setup,ele_file)!=0)
    {
      std::cerr<<"Error while performing pre-processing"<<std::endl;
      return EXIT_FAILURE;    
    }

    if(ad.info)
    {
      std::cout<<std::endl;  
      std::cout<<"Ending pre-processing data "<<std::endl;
    }
  }

  //! Now perform the CA2D flood modelling if requested.   
  if(ad.WCA2D)
  {
    if(ad.info)
      std::cout<<std::endl<<"Starting WCA2D flood modelling "<<std::endl;

    if(WCA2D(ad,setup,eg,res,wles,ies,tps,rgs)!=0)
    {
      std::cerr<<"Error while performing WCA2D flood modelling"<<std::endl;
      return EXIT_FAILURE;    
    }
    
    if(ad.info)
    {
      std::cout<<std::endl;  
      std::cout<<"Ending WCA2D flood modelling "<<std::endl;
    }
  }

  //! Now perform the post-processing
  if(ad.post_proc)
  {
    if(ad.info)
      std::cout<<std::endl<<"Starting post-processing data "<<std::endl;

    if(postProc(ad,setup,eg,tps,rgs)!=0)
    {
      std::cerr<<"Error while performing post-processing"<<std::endl;
      return EXIT_FAILURE;    
    }

    if(ad.info)
    {
      std::cout<<std::endl;  
      std::cout<<"Ending post-processing data "<<std::endl;
    }
  }

  // Close the caAPI.
  CA::finalise2D();

  return EXIT_SUCCESS;
}





