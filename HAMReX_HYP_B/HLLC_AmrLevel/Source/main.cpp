
#include <new>
#include <iostream>
#include <iomanip>

#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_EB2.H>

//#include "structdefs.H"
//#include "funcdefs.H"

using namespace amrex;

//void initialize_EB2 (const Geometry& geom, const int required_level, const int max_level);
void makeEmbeddedBoundary(Geometry&, const int&);

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    Real dRunTime1 = amrex::second();
    
    /* -----------------------------------------------------
     * Global simulation parameters required for main 
     * coarseTimeStep loop, pulled from inputs file 
     * through ParmParse
     * -----------------------------------------------------*/
    
    ParmParse pp;

    int  max_step;
    pp.get("nsteps_max", max_step);
    Real strt_time;
    pp.get("startT", strt_time);
    int NGROW;
    pp.get("NGROW", NGROW);
    
    std::string Test;
    pp.get("testcase", Test);
    
    ParmParse ppTest(Test);
    Real stop_time; 
    ppTest.get("finalT", stop_time);


    if (strt_time < 0.0) {
        amrex::Abort("MUST SPECIFY a non-negative strt_time");
    }

    if (max_step < 0 && stop_time < 0.0) {
	amrex::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

    {
	
	Amr amr;
		AmrLevel::SetEBSupportLevel(EBSupport::full);
        //AmrLevel::SetEBMaxGrowCells(NGROW,4,2); //Q: 4 and 2 ???  
		makeEmbeddedBoundary(amr.Geom(amr.maxLevel()), amr.maxLevel());
		amrex::Print() << "max refinement level = " << amr.maxLevel() << "\n";
	
	amr.init(strt_time,stop_time);

	while ( amr.okToContinue() &&
  	       (amr.levelSteps(0) < max_step || max_step < 0) &&
	       (amr.cumTime() < stop_time || stop_time < 0.0) )

	{
	    //
	    // Do a coarse timestep.  Recursively calls timeStep()
	    //
	    amr.coarseTimeStep(stop_time);
	}

	// Write final checkpoint and plotfile
	if (amr.stepOfLastCheckPoint() < amr.levelSteps(0)) {
	    amr.checkPoint();
	}

	if (amr.stepOfLastPlotFile() < amr.levelSteps(0)) {
	    amr.writePlotFile();
	}

    }

    Real dRunTime2 = amrex::second() - dRunTime1;

    ParallelDescriptor::ReduceRealMax(dRunTime2, ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Run time = " << dRunTime2 << std::endl;

    amrex::Finalize();

    return 0;
}
