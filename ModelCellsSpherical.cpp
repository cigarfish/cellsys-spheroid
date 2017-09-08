//#include <cstdio>
#include<iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cfloat>
#include <fstream>
#include <time.h>
#include <ctime>
#include <functional>
#include <random>
//
#include <malloc.h>
#include <stdio.h>

using namespace std::placeholders;
using namespace std;

#include "../utils/float_utils.h"
#include <iomanip>
#include <fstream>
#include <iostream>
#include <limits>

#include "ModelCellsSpherical.h"
#include "../../../tools/model/CSModelTools.h"
#include "../../Core.h"
#include "../../../tools/parameters/CSParameterContext.h"
#include "../../../tools/parameters/CSParameterContextTemporary.h"
#include "../../../tools/model/BoundingBoxList.h"
#include "../../Interactions/CSInteractionHertz.h"
#include "../../Interactions/CSInteractionJKR.h"
#include "../../Interactions/CSInteractionHertzEnergy.h"
#include "../../Interactions/CSInteractionJKREnergy.h"
#include "../../Interactions/CSInteractionFrictionMatrix.h"
#include "../../../gui/QCSSimulationThread.h"
#include "../../Elements/ModelElementBarrierTriangle.h"
//#include "../../Elements/ModelElementVesselSphere.h"
#include "../../../tools/dataIO/vtp/CSVTPWriter.h"
#include "../../../tools/dataIO/xml/CSXMLWriter.h"
#include "../../utils/macros.h"

#include "../../../tools/math/CGP.h"
#include "../../../tools/math/SimpleSolver.h"

#include "../../BasicDatatypes/GraphSphere.h"

#include <QtCore>

#include <H5Cpp.h>
#include "../../../tools/dataIO/hdf5/CSHDF5Helpers.h"
#include <time.h>
#include "../../../tools/math/mathematics.h"
#include <locale.h>

const std::string ModelCellsSpherical::xmlType = "ModelCellsSpherical";

const std::string ModelCellsSpherical::mContactModelNames[] =
{
	"Hertz Model",
	"JKR Model",
	""
};

const std::string ModelCellsSpherical::mLobuleShapeNames[] =
{
	"None",
	"Hexagonal",
	"Quadric",
	""
};

// Preliminary:  fixed CellSphericalPolar 'region angle'.
//   This value corresponds to an adhesive area of 33% of the cells surface
#define POLAR_REGION_ANGLE .841

/*
* Constructor
*/

ModelCellsSpherical::ModelCellsSpherical()
	: CSModel(3),
	mScenario(ScenarioSingleCell),
	mContactModel(ContactModelHertz),
	mUseDynamicalTimeSteps(true),
	mBloodVesselNetwork(false),
	mTimeScalingFactor(2),
	mpInteractionHertz(nullptr),
	mpInteractionJKR(nullptr),
	mpInteractionFrictionMatrix(nullptr),
	mpJKRElementDone(nullptr),
	mJKRElementDoneSize(0),
	mpVelocities(nullptr),
	mProblemAllocationSize(0),
	mpParameters(nullptr),
	mpFrictionMatrices(nullptr),
	mpGraphBloodVesselNetwork(nullptr),
	mUsePolarCells(false),
	observeDivision(false),
	initial_radius(-1),
	nof_cells(-1),
	minimum_cell_distance(0),
	scen_two_cells_dist(1),
	mCellsAtLastUpdate(0),
	mUseCCFriction(true),
	mMaximumDisplacementCellRadius(0.1),
	mMaximumDisplacementSquared(0.0),
	mStateChoice(0),
	mSolver(nullptr)
{
	// Console output
	core->tools->output->consoleModelCellsSpherical_Simulation << "Init model.\n";
	if ((*core->tools->output).debugMonolayer) (*core->tools->output).logfile
		<< "ModelMonolayer::ModelMonolayer(). Create model\n";

	// Unit conversions, Biological -> Dimensionless values for modeling etc.
	biolink = new BiologyLink(0); // Init for ModelCellsSpherical model

								  // register all parameters in our parameter context
	RegisterParameters();
	lastColormodeCells = 0;

	// for test
	is2D = true;
	//is2D = false;
}

ModelCellsSpherical::~ModelCellsSpherical()
{
	if (mpFrictionMatrices)
		free(mpFrictionMatrices);

	if (mpJKRElementDone)
		free(mpJKRElementDone);

	if (mSolver)
		delete mSolver;
}

/*
 *  region Setup
 */

void ModelCellsSpherical::ResetDerived()
{
	no_cells = 0;
	// Preliminary:  hardcode if to use dynamical step size adaptation:
	//    mUseDynamicalTimeSteps = true;//true;
	mTimeScalingFactor = 2;
	//timeStep = 1;
	//Decide if use cell-cell friction or not
	//mUseCCFriction=true;//default is true;
	lastTimeConsoleUpdate = 0.;
	step = 0;

	UpdateParametersFromBioLink();

	if (mSolver) delete mSolver;

	if ((*core->tools->output).debugMonolayer) (*core->tools->output).logfile
		<< "ModelMonolayer::Reset(). Reset monolayer model\n";

	int dims = 3;
	if (is2D)
	{
		core->tools->output->consoleModelCellsSpherical_Simulation
			<< "Reset model (strictly monolayer).\n";
		dims = 2;
	}
	else
	{
		core->tools->output->consoleModelCellsSpherical_Simulation
			<< "Reset model (tumor spheroid).\n";
	}

	//mRandom.Init();

	if (mpInteractionFrictionMatrix)
		delete mpInteractionFrictionMatrix;
	mpInteractionFrictionMatrix = nullptr;

	if (mpInteractionHertz)
		delete mpInteractionHertz;
	mpInteractionHertz = nullptr;

	if (mpInteractionJKR)
		delete mpInteractionJKR;
	mpInteractionJKR = nullptr;

	if (mpFrictionMatrices)
	{
		free(mpFrictionMatrices);
		mpFrictionMatrices = nullptr;
	}

	if (mpInteractionFrictionMatrix)
	{
		free(mpInteractionFrictionMatrix);
		mpInteractionFrictionMatrix = nullptr;
	}

	std::cout << "removed former interactions" << std::endl;
	mProblemAllocationSize = 0;

	std::cout << "kill blood vessels" << std::endl;
	if (mpGraphBloodVesselNetwork != nullptr)
		delete mpGraphBloodVesselNetwork;
	mpGraphBloodVesselNetwork = nullptr;

	#pragma region Remove existing cells
	if (cells.size() > 0)
	{
		std::vector<CellSpherical *>::iterator cellIt;
		for (cellIt = cells.begin(); cellIt != cells.end(); ++cellIt)
		{
			delete (*cellIt);
		}
		cells.erase(cells.begin(), cells.end());
		cells.clear();
	}

	// Remove cells from arena
	std::cout << "kill arena" << std::endl;
	mpArena->clear();

	std::cout << "============================================\n";
	std::cout << "Intialize simulation using scenario " << (int)mScenario << std::endl;

	#pragma region Init without reading cells ...
	if (!this->mReadCells)
	{
		if (mScenario == ScenarioSingleCell)
		{
			if (mUsePolarCells)
			{
				AddPolarCell(0, 0, 0);
			}
			else
			{
				AddCell(0, 0, 0);
				//AddCell(0, 0, 1);
				// test cell Divide()
				//CellSpherical * newCell = static_cast<CellSpherical *>(cells[0]->Divide());
				//AddCell(newCell);
				//double nc1x = newCell->position.x + 1.4;
				//double nc1y = newCell->position.y + 1.3;
				//double nc1z = newCell->position.z - 1.4;
				//AddCell(nc1x, nc1y, nc1z);
			}
		}
		else if (mScenario == ScenarioEmbeddingMedium)
		{
			std::cout << "Embedding Medium" << endl;
			double populationRadiusQ = mPopulationRadius * mPopulationRadius;

			ModelElementHollowSphere * hull = new ModelElementHollowSphere(0, 0, 0);
			hull->mRadius = mPopulationRadius;
			hull->mStatic = true;
			cells2->add(hull);
			mpArena->addObject(hull->GLObject());

			// Init a spherically arranged population of cells
			for (double x = -mPopulationRadius; x < mPopulationRadius; x = x + mPopulationInitialDistance)
			{
				for (double y = -mPopulationRadius; y < mPopulationRadius; y = y + mPopulationInitialDistance)
				{
					if (is2D)
					{
						if (x * x + y * y <= populationRadiusQ)
						{
							AddCell(x, y, 0);
							cells.back()->setState(Cell::StateQuiescent);
						}
					}
					else
					{
						for (double z = -mPopulationRadius; z < mPopulationRadius; z = z + mPopulationInitialDistance)
						{
							if (x*x + y*y + z*z <= populationRadiusQ)
							{
								AddCell(x, y, z);
								cells.back()->setState(Cell::StateQuiescent);
							}
						}
					}
				}
			}
			double minDistQ = 10e10;
			int minCellIndex = 0;

			for (unsigned int i = 0; i < this->cells.size(); i++)
			{
				this->cells[i]->SetCellSubType(SubType::ETERNAL_QUIESCENT);
				double distQ = this->cells[i]->position.x*this->cells[i]->position.x +
					this->cells[i]->position.y*this->cells[i]->position.y +
					this->cells[i]->position.z*this->cells[i]->position.z;

				if (distQ < minDistQ)
				{
					minDistQ = distQ;
					minCellIndex = i;
				}
			}

			this->cells[minCellIndex]->SetCellSubType(SubType::NORMAL);
		}
		else if (mScenario == ScenarioSpherePacking)
		{
			double radius_capsule = 40e-5;
			this->spherePacking(radius_capsule / this->biolink->length_scale);
			for (unsigned int i = 0; i < this->cells.size(); i++)
			{
				this->cells[i]->setState(Cell::StateQuiescent);
			}
		}
		else
		{
			std::cout << "Scenario not recognized!" << endl;
		}
		mNextKillTime = 0.0;
	}
	else
	{
		// read in 3D
		cells2->setDimensions(3);
		is2D = false;

		fprintf(stderr, "%s \n", mBloodVesselNetworkPath.c_str());
		this->AddReadCells(mBloodVesselNetworkPath);

		for (unsigned int i = 0; i < this->cells.size(); i++) 
		{
			double dist = std::sqrt(this->cells[i]->position.x * this->cells[i]->position.x +
									this->cells[i]->position.y * this->cells[i]->position.y);
			if (dist < 4.0)
			{
				this->mCellsToDelete.push_back(this->cells[i]);
			}
		}
		mKillInterval = 300.0;
		mKillPeriod = 86400.0;
		mKillEnd = mKillPeriod;
		mKillRate = (mKillInterval / mKillPeriod) * mCellsToDelete.size();

		mNextKillTime = mKillInterval;
	}

	// note is overwritten by read data
	nextProliferationUpdate = 0.;

#pragma region Vessel network init (Prelim: NOT HERE!!)

	//vessel network
	//
	// TODO MOVE TO A GENERATOR!!
	//
	//
	if (mBloodVesselNetwork)
	{
		AddBloodVesselNetwork(mBloodVesselNetworkPath, 2);

		// Blood vessel networks always in 3D
		is2D = false;

		for (unsigned int i = 0; i < mpGraphBloodVesselNetwork->mvNode.size(); i++)
		{
			mpGraphBloodVesselNetwork->mvNode[i]->color.red = 0;
			mpGraphBloodVesselNetwork->mvNode[i]->color.green = 0;
			mpGraphBloodVesselNetwork->mvNode[i]->color.blue = 1;
		}
	}

	//
	// SAME HERE MOVE THIS SOMEWHERE SEPARATE
	//
	mLobuleShape = (LobuleShape)((CSParameterChoice *)mpParameters->
		findParameter("Lobule Shape")->dataPointer())->currentIndex();
	if (mLobuleShape != 0)
	{
		AddShape(mlobule_radius, mlobule_height, 1, 1, 0, 0.5, 1);
	}

	mContactModel
		= (ContactModel)((CSParameterChoice *)mpParameters->
			findParameter("Contact Model")->dataPointer())->currentIndex();

	enableSimulation = false;

	//core->tools->output->consoleModelCellsSpherical_Simulation
	//	<< "Rest is done!\n";
}

void ModelCellsSpherical::SetupSimulation()
{
	// the maximum allowed displacement for dynamic step size adaptation:
	// set to 10% of a cell's initial radius, but squared.
	mMaximumDisplacementSquared = mMaximumDisplacementCellRadius
		* mMaximumDisplacementCellRadius
		* defaultInitialCellRadius
		* defaultInitialCellRadius;

	// Remember simulation start for progress bar
	//
	// some of the commmon stuff could be put into CSModel XXX
	//
	simulateFromDays = biolink->getTimeInDays(time);

	std::cout << "sim from = " << simulateFromDays << " to " << simulateUntilDays << std::endl;
	// Target simulation time has already been reached -> No initial observation
	if (simulateFromDays >= simulateUntilDays) return;

	mpSimulationThread->setUpdateInterval(50);

	if (mBloodVesselNetwork)
	{
		mpGraphBloodVesselNetwork->setInitialOverlaps();
		mpGraphBloodVesselNetwork->setLength0();
		mpGraphBloodVesselNetwork->setStatic(ModelElementVesselSphere::CentralVein);
		mpGraphBloodVesselNetwork->setStatic(ModelElementVesselSphere::PortalVein);
	}

	mpSimulationThread->setUpdateInterval(50);

	enableSimulation = true;

	if (enableObservation)
	{
		core->tools->output->logfile << " In SetupSimulation: " << name << "\n";;

		// Set next observation time
		nextObservationTime = simulateFromDays + observeEveryDays;
	}

	//
	// Initializing interactions:
	//
	mpInteractionFrictionMatrix = new CSInteractionFrictionMatrix();
	std::cout << "CREATING NEW FRICTION MATRIX" << std::endl;
	
	// Connecting model parameters to mpInteractionFrictionMatrix
	std::cout << "Simulation Parameters" << std::endl;
	std::cout << "Friction Matrix" << std::endl;
	std::cout << "gamma Cells || " << gammaCellsParallel << ", |_ " << gammaCellsPerpendicular << std::endl;
	mpInteractionFrictionMatrix->mpGammaCellsParallel = &gammaCellsParallel;
	mpInteractionFrictionMatrix->mpGammaCellsPerpendicular = &gammaCellsPerpendicular;

	if (mContactModel == ContactModelHertz)
	{
		mpInteractionHertz = new CSInteractionHertz();

		// Connecting model parameters to mpInteractionHertz
		mpInteractionHertz->mIs2D = is2D;
		mpInteractionHertz->mpSingleBondEnergy = &singleBondEnergy;
		mpInteractionHertz->mpAdhesionDensity = &adhesionDensity;

		// Connecting 'output' of mpInteractionHertz to 'input' of
		// mpInteractionFrictionMatrix from the Hertz force calculation:
		mpInteractionFrictionMatrix->mpDistance = &mpInteractionHertz->mDistance;
		mpInteractionFrictionMatrix->mpContactArea = &mpInteractionHertz->mContactArea;
	}
	else
	{
		mpInteractionJKR = new CSInteractionJKR();

		// Connecting model parameters to mpInteractionJKR
		mpInteractionJKR->mIs2D = is2D;
		mpInteractionJKR->mpSingleBondEnergy = &singleBondEnergy;
		mpInteractionJKR->mpAdhesionDensity = &adhesionDensity;

		// Connecting 'output' of mpInteractionHertz to 'input' of
		// mpInteractionFrictionMatrix from the JKR force calculation:
		mpInteractionFrictionMatrix->mpDistance = &mpInteractionJKR->mDistance;
		mpInteractionFrictionMatrix->mpContactArea = &mpInteractionJKR->mContactArea;
	}

	UpdateParametersForAllCells();
	printParameters();

	// MOVE THE SOLVER TO A TYPE OF PLUGIN THAT IS JUST SIMPLY HANDLED IN
	// CSMODEL
	if (mUseCCFriction)
	{
		std::cout << "Using CGP:" << std::endl;
		mSolver = new CGP();
		if (is2D) mSolver->set2D();
	}
	else
	{
		std::cout << "Using SimpleSolver:" << std::endl;
		mSolver = new SimpleSolver();
	}

	mSolver->setBoundingBoxList(cells2);

	//std::stringstream ss;
	//mpParameters->dump(ss);
	//std::cout << "Parameters+" << ss.rdbuf() << std::endl;
}
// end Setup

/*
* Simulation
*/
//! Calculates the next model state t + dt

void ModelCellsSpherical::SimulateDiffusionTest()
{
	simulateCommon();

	assert(cells2->size() == 1);
	assert(cells.size() == 1);

	std::ofstream F;
	std::string outputFileName = mOutputPath + mOutputPrefix + "_testdiffusion.dat";
	F.open(outputFileName.c_str(), std::ios::out | std::ios::app);
	auto cell = cells[0];
	F << "\n" << biolink->getTimeInDays(time)
		<< "\t" << time
		<< "\t" << cell->mRadius;
	F << "\t" << cell->maxoverlap;
	F << "\t" << cell->accumulatedForceAbsolute;
	F << "\t" << biolink->getForceInNanoNewton(cell->directedForce.x);
	F << "\t" << biolink->getLengthInMicrometers(cell->position.x)
		<< "\t" << biolink->getLengthInMicrometers(cell->position.y)
		<< "\t" << biolink->getLengthInMicrometers(cell->position.z);
	F << "\t" << biolink->getForceInNanoNewton(cell->mLangevinForce.x)
		<< "\t" << biolink->getForceInNanoNewton(cell->mLangevinForce.y)
		<< "\t" << biolink->getForceInNanoNewton(cell->mLangevinForce.z);

	F.close();

	time += timeStep;
}

void ModelCellsSpherical::SimulateTwoCells()
{
	simulateCommon();

	assert(cells2->size() == 2);
	assert(cells.size() == 2);

	// TODO move this crap to helpers seriously this is done several times
	// all over the place!!!
	std::ofstream F;
	std::string outputFileName = mOutputPath + mOutputPrefix + "_twocells.dat";
	F.open(outputFileName.c_str(), std::ios::out | std::ios::app);
	double dist = (cells[0]->position - cells[1]->position).Norm();
	F << "\n" << time << "\t" << dist << "\t" << cells[0]->mRadius << "\t" << cells[1]->mRadius;
	F << "\t" << cells[0]->maxoverlap << "\t" << cells[1]->maxoverlap;
	F << "\t" << cells[0]->accumulatedForceAbsolute << "\t" << cells[1]->accumulatedForceAbsolute;
	F << "\t" << biolink->getForceInNanoNewton(cells[0]->directedForce.x)
		<< "\t" << biolink->getForceInNanoNewton(cells[1]->directedForce.x);
	F << "\t" << cells[0]->position.x << "\t" << cells[0]->position.y << "\t"
		<< cells[0]->position.z;
	F << "\t" << cells[1]->position.x << "\t" << cells[1]->position.y << "\t"
		<< cells[1]->position.z;
	F << "\t" << CSModelTools::GetDistance3D(cells[0]->position, cells[1]->position);

	F.close();
	std::cout << "Overlap Threshold ";
	for (CellSpherical* cell : cells)
		std::cout << " c= " << overlap_threshold * cell->mRadius;
	std::cout << std::endl;
	for (CellSpherical * cell : cells)
		cell->maxoverlap = 100;
	//mpParaview->exec(this);

	time += timeStep;
}

void ModelCellsSpherical::simulateWithNutrients()
{}

void ModelCellsSpherical::SimulateGeneric()
{
	int stepstop = 7100;
	if (step <= stepstop)
	{
		// bug file
		//std::string outputFileName = "debugFile.txt";
		//std::ofstream FO;
		//FO.open(outputFileName.c_str(), std::ios::out | std::ios::app);
		//FO << "ModelCellsSpherical.cpp / SimulateGeneric()	step " << step << ", lastcolormode: " << lastColormodeCells << endl;
		//
		simulateCommon();

		double mCellularBeforeToUpdate = 0.05;
		if (cells.size() - mCellsAtLastUpdate > mCellularBeforeToUpdate * (cells.size())
			|| biolink->getTimeInDays(time - lastTimeConsoleUpdate) > 0.5)
		{
			mCellsAtLastUpdate = cells.size();
			lastTimeConsoleUpdate = time;
		}

		GrowAndDivide();

		// test: add bars representing cell-cell interaction
		bool takebar = false;
		if (step == stepstop) { takebar = true; }
		if (takebar)
		{
			BoundingBoxList::iterator cellIterator = cells2->begin();
			while (cellIterator != cells2->end())
			{
				ModelElement *object_1 = *cellIterator;
				CSListContainer< unsigned long > & contactsElements = object_1->mIntersectingList;//mContacts;
				unsigned long * iterator1;
				for (iterator1 = contactsElements.begin();
					iterator1 != contactsElements.end();
					++iterator1)
				{
					ModelElement * cell_2e = cells2->element(*iterator1);
					ARGBColor color_be; color_be.alpha = 1.0; color_be.red = 1.0; color_be.green = 0.1; color_be.blue = 0.1;
					CSGLBar *neighbor_be = new CSGLBar(&object_1->position, &cell_2e->position, color_be, 0.5);
					mpArena->addObject(neighbor_be);
				}
				++cellIterator;
			}
		}

		//FO << "ModelCellsSpherical.cpp / SimulateGeneric()	after grow&divide. Number of cells: " << cells.size() << endl;
		//if (cells.size() > 0)
		//{
		//	for (int ii = 0; ii < cells.size(); ii++)
		//	{
		//		double cx = cells[ii]->position.x;
		//		double cy = cells[ii]->position.y;
		//		double cz = cells[ii]->position.z;
		//		double cr = cells[ii]->mRadius;
		//		FO << "ModelCellsSpherical.cpp / SimulateGeneric()	cell " << ii+1 << ": " << cx << ", " << cy << ", " << cz << ", radius: " << cr << endl;
		//	}
		//}

		//FO << "ModelCellsSpherical.cpp / SimulateGeneric()	time: " << time << endl;
		//FO.close();
	}

	//std:cout << "step:" << step << std::endl;
	//if (cells.size() > 8)
	//    assert(false);

	time += timeStep;
}

void ModelCellsSpherical::simulateCommon()
{
	// bug file
	//std::string outputFileName = "debugFile.txt";
	//std::ofstream FO;
	//FO.open(outputFileName.c_str(), std::ios::out | std::ios::app);
	//FO << "ModelCellsSpherical.cpp / simulateCommon()	starts" << endl;
	//
	for (CellSpherical* cell : cells)
	{
		if (abs(cell->position) > 1e2)
		{
			std::cout << "Cell:" << cell->SimObjId << " pos:" << cell->position
				<< " force:" << cell->directedForce
				<< " daughter:" << std::boolalpha << cell->mDaughterCell
				<< " necrotic:" << std::boolalpha << cell->isNecrotic()
				<< " Q:" << std::boolalpha << cell->isQuiescent()
				<< " D:" << std::boolalpha << cell->isDividing()
				<< std::endl;
			assert(false);
		}
	}

	InitForces();
	UpdateInteractions();

	//FO << "ModelCellsSpherical.cpp / simulateCommon()	after InitForce & UpdateInteractions." << endl;

	if (mUseDumbbell)
	{
		DumbbellMetropolis();
		//FO << "ModelCellsSpherical.cpp / simulateCommon()	after dumbbellMetropolis." << endl;
	}
	else
	{
		//FO << "ModelCellsSpherical.cpp / simulateCommon()	no dumbbellMetropolis." << endl;
	}

	UpdateSingleCellEffects();

	//FO << "ModelCellsSpherical.cpp / simulateCommon()	after updateSingleCellEffects." << endl;

	if (mpGraphBloodVesselNetwork != nullptr)
	{
		mpGraphBloodVesselNetwork->calcSpringForce(1000.0);
		//FO << "ModelCellsSpherical.cpp / simulateCommon()	after calcSpringForce." << endl;
	}
	else
	{
		//FO << "ModelCellsSpherical.cpp / simulateCommon()	no calcSpringForce." << endl;
	}

	setFrictionMatrix();

	//FO << "ModelCellsSpherical.cpp / simulateCommon()	after setFrictionMatrix." << endl;

	//std::cout << "SIMULATE COMMON " << biolink->getTimeInDays(time);
	//std::cout << " timeStep: " << biolink->getTimeInSeconds(timeStep);
	//std::cout << "(s) max(v^2): " << mVelocityMaxSquared;
	//std::cout << " d_m: " << mMaximumDisplacementSquared;
	//std::cout << " t_m: " << biolink->getTimeInSeconds(mDynamicTimeStepMin);

	// Compute new position based on forces

	if (mUseDynamicalTimeSteps)
	{
		double timeStepFactor = 1;
		// a flag to register the scale direction of a possible leading scaling
		// to prevent infinite loops (down-up-down... scalings
		int scaled = 0;
		unsigned int scaleSquared = mTimeScalingFactor*mTimeScalingFactor;
		std::size_t no_attempts{ 0 };

		while (true)
		{
			no_attempts++;
			mSolver->solve();
			mVelocityMaxSquared = mSolver->getVelocityMaxSquared();
			double displacementSquared = mVelocityMaxSquared * timeStep * timeStep;

			if (displacementSquared <= mMaximumDisplacementSquared)
			{
				//std::cout << " scaled D=" << scaleSquared*displacementSquared;
				if (scaleSquared*displacementSquared <= mMaximumDisplacementSquared)
				{
					if (scaled < 0)
						break;

					// while ( 4*displacementSquared <= mMaximumDisplacementSquared )
					//check range of dynamic timestep
					if (timeStep < mDynamicTimeStepMax)
					{
						displacementSquared *= scaleSquared;
						timeStepFactor *= mTimeScalingFactor;
						timeStep *= mTimeScalingFactor;
						//std::cout << " doubling timeStep to "
						//    << biolink->getTimeInSeconds(timeStep) << "(s)"
						//    << std::endl;

						rescaleLangevinContrib(mTimeScalingFactor);
					}
					else
						break;
					scaled = 1;

					continue;
				}

				break;
			}

			//check range of dynamic timestep
			if (timeStep > mDynamicTimeStepMin)
			{
				//std::cout << "timeStep: " << timeStep << " Dyn:" << mDynamicTimeStepMin << std::endl;
				displacementSquared /= scaleSquared;
				timeStepFactor /= mTimeScalingFactor;
				timeStep /= mTimeScalingFactor;
				//std::cout <<" reducing timeStep to: "
				//    << biolink->getTimeInSeconds(timeStep) << "(s)";

				rescaleLangevinContrib(mTimeScalingFactor);
			}
			else
				break;

			// prevent infinite loop.
			if (scaled > 0)
				break;

			scaled = -1;
			//std::cout << "reducing timeStep " << timeStep << std::endl;

			if (no_attempts > 3)
				std::cerr << "WARNING(" << time << "): simulate common number of attempts larger 3! Scaled:"
				<< scaled << " new time step is "
				<< biolink->getTimeInSeconds(timeStep) << "(s)."
				<< " max: " << std::sqrt(mMaximumDisplacementSquared)
				<< " max^2: " << mMaximumDisplacementSquared
				<< std::endl;
		}

	}
	else
	{
		mSolver->solve();
	}

	//FO << "ModelCellsSpherical.cpp / simulateCommon()	after solveMatrix." << endl;
	//FO.close();

	updatePosition();
}

void ModelCellsSpherical::updatePosition()
{
	// bug file
	//std::string outputFileName = "debugFile.txt";
	//std::ofstream FO;
	//FO.open(outputFileName.c_str(), std::ios::out | std::ios::app);
	//FO << "ModelCellsSpherical.cpp / updatePosition()	starts" << endl;
	//update position of all elements
	for (unsigned int i = 0; i < cells2->size(); ++i)
	{
		ModelElement * element = cells2->element(i);

		if (!element)
			continue;

		if (element->mStatic)
		{
			//FO << "ModelCellsSpherical.cpp / updatePosition()		element " << i << " is static!" << endl;
			continue; // skip rest of for loop
		}


					  // ToDo:  Code pertaining dumbbell should be handled by Cells.
					  // ToDo:  Hierarchically handled bitmap for types, so that derived
					  //        types will not be needed for this and similar tests.
		if (mUseDumbbell)
		{
			//cout << "DUMBBELL";
			if (element->mType == ModelElement::TypeCellSpherical ||
				element->mType == ModelElement::TypeCellSphericalPolar)
			{
				CellSpherical * cell = static_cast<CellSpherical *>(element);
				if (cell->mpDivisionPartner)
				{
					double *daughterVelocity = mpVelocities + 3 * element->mGlobalIndex;
					double *motherVelocity = mpVelocities + 3 * cell->mpDivisionPartner->mGlobalIndex;

					daughterVelocity[0] = 0.5 * (daughterVelocity[0] + motherVelocity[0]);
					daughterVelocity[1] = 0.5 * (daughterVelocity[1] + motherVelocity[1]);
					daughterVelocity[2] = 0.5 * (daughterVelocity[2] + motherVelocity[2]);

					motherVelocity[0] = daughterVelocity[0];
					motherVelocity[1] = daughterVelocity[1];
					motherVelocity[2] = daughterVelocity[2];
				}
			}
		}

		//const double MaxPositionChange =
		//    defaultInitialCellRadius * mMaximumDisplacementCellRadius;

		unsigned int j = 3 * i;

		//if (std::abs(mpVelocities[j] * timeStep) > 1.25*MaxPositionChange
		//    || std::abs(mpVelocities[j+1] * timeStep) > 1.25*MaxPositionChange
		//    || std::abs(mpVelocities[j+2] * timeStep) > 1.25*MaxPositionChange)
		//{
		//    std::cout << "EXIT(" << time << ") timeStep=" << timeStep
		//        << " update x=" << std::abs(mpVelocities[j] * timeStep)
		//        << " update y=" << std::abs(mpVelocities[j+1] * timeStep)
		//        << " update z=" << std::abs(mpVelocities[j+2] * timeStep)
		//        << " v_z=" << mpVelocities[j+2];
		//    std::cout <<  " max: " << MaxPositionChange << std::endl;
		//    assert(false);
		//}

		element->position.x += mpVelocities[j] * timeStep;
		element->position.y += mpVelocities[j + 1] * timeStep;

		//std::cout << "element: " << element->mGlobalIndex
		//    << " x=" << element->position.x
		//    << " y=" << element->position.y
		//    << " z=" << element->position.z
		//    << " dX=(" << mpVelocities[j] << ", " << mpVelocities[j+1]
		//        << ", " << mpVelocities[j+2] << ")"
		//    << std::endl;

		if (!is2D)
			element->position.z += mpVelocities[j + 2] * timeStep;
	}
	//FO << "ModelCellsSpherical.cpp / updatePosition()	ends" << endl;
	//FO.close();
	// end Compute new position based on forces
}

void ModelCellsSpherical::Simulate()
{
	step++;
	// TODO is this too slow?
	//if (biolink->getTimeInDays(time - lastTimeConsoleUpdate) > 1.0)
	//{
	//    //console->info("SimName: {} Time: {} number of cells: {} number of cells2: {}.",
	//    //              xmlName, biolink->getTimeInDays(time), cells.size(), cells2->size());
	//    lastTimeConsoleUpdate = time;
	//    compute();
	//    //printParameters();
	//}

	switch (mScenario)
	{
	case ScenarioReadFromData:
	case ScenarioSingleCell:
	case ScenarioEmbeddingMedium:
	case ScenarioSpherePacking:
	case ScenarioSphere:
		SimulateGeneric();
		break;
	case ScenarioGlucoseDiffusion:
		simulateWithNutrients();
		break;
	case ScenarioTestDiffusion:
		SimulateDiffusionTest();
		break;
	case ScenarioTwoCells:
		SimulateTwoCells();
		break;
	default:
		unreachable("Unknown simulation scenario!");
		break;
	}
}

void ModelCellsSpherical::AddDirectedMotion(CellSpherical *cell)
{
	// ToDo:  change into a parameter
	double someFactor = 190;
	Vector3f directedForceVector = -cell->position;
	// linear dependency on distance:
	double distance = directedForceVector.Normalize(); // |directedForcVector| = 1
	directedForceVector.y = 0;
	directedForceVector.z = 0;
	double directedMotionForce = someFactor * (time - 100) * distance;
	directedForceVector *= directedMotionForce;
	cell->directedForce += directedForceVector;
}

void ModelCellsSpherical::Simulate(int numberOfSimulationSteps)
{
	for (int i = 1; i<numberOfSimulationSteps; i++)
		Simulate();
}

void ModelCellsSpherical::RemoveCell(CellSpherical * cell)
{
	//console->info("Removing cell {}.", cell->mGlobalIndex);
	std::vector<CellSpherical *>::iterator found
		= std::find(cells.begin(), cells.end(), cell);

	cells.erase(found);
	cells2->remove(cell);
	mpArena->removeObject(cell->GLObject());

	delete cell;
}

// Grows all cells (using the given time delta) and applies cell divisions
void ModelCellsSpherical::GrowAndDivide()
{
	//remove cells which are labed to delete
	if (mNextKillTime && time >= mNextKillTime)
	{
		if (time >= mKillEnd)
		{
			mNextKillTime = 0.;
			while (mCellsToDelete.size())
			{
				CellSpherical * cell = mCellsToDelete.back();
				mCellsToDelete.pop_back();
				RemoveCell(cell);
			}
		}
		else
		{
			std::vector< CellSpherical * > killedCells;

			double probabilityToDie = mKillRate / mCellsToDelete.size();

			std::vector<CellSpherical * >::iterator doomedCellsIt;
			for (doomedCellsIt = mCellsToDelete.begin();
				doomedCellsIt != mCellsToDelete.end();
				++doomedCellsIt)
			{
				double dieThrow = mRandom.GetRandomUniform01();
				if (dieThrow < probabilityToDie)
					killedCells.push_back(*doomedCellsIt);
			}

			std::vector<CellSpherical * >::iterator found;
			for (doomedCellsIt = killedCells.begin();
				doomedCellsIt != killedCells.end();
				++doomedCellsIt)
			{
				found = std::find(mCellsToDelete.begin(), mCellsToDelete.end(), *doomedCellsIt);
				if (found != mCellsToDelete.end())
					mCellsToDelete.erase(found);

				RemoveCell(*doomedCellsIt);
			}

			if (!mCellsToDelete.size())
				mNextKillTime = 0.;
			else
				mNextKillTime += mKillInterval; // every mKillInterval
		}
	}

	// For all cells
	for (std::size_t i = 0; i< cells.size(); i++)
	{
		CellSpherical * cell = cells[i];

		// Cells may wake up, but can only become quiescent after completing the
		// cell cycle
		cell->WakeUp();

		if (cell->isQuiescent() || cell->isNecrotic())
			continue;

		// Grow the cell
		cell->Grow();

		// Divide the cell if possible
		if (cell->CanDivide())
		{
			CellSpherical * newCell = static_cast<CellSpherical *>(cell->Divide());
			AddCell(newCell);
		}
	}

	if (time >= nextProliferationUpdate)
		nextProliferationUpdate += 300.;
}

//! Adds and initializes a new cell
//! This method is used by the model to add cells at the given coordinates.
void ModelCellsSpherical::AddCell(double x, double y, double z)
{
	CellSpherical * newCell = new CellSpherical(x, y, z, defaultInitialCellRadius);
	InitCell(newCell);
	AddCell(newCell);
}

//! Adds and initializes a new cell
//! This method is used by the model to add cells at the given coordinates.
void ModelCellsSpherical::AddPolarCell(double x, double y, double z)
{
#pragma message("WARNING: use proper constructor here")
	CellSphericalPolar * newCell = new CellSphericalPolar(x, y, z);

	InitCell(newCell);

	double lx, ly, lz;
	mRandom.GetRandomUnitVector(&lx, &ly, &lz);

	newCell->mPolarDirection.Set(lx, ly, lz);
	newCell->mNewPolarDirection.Set(lx, ly, lz);

	// preliminary: hardcode opening angle (~33% adhesive surface area)
	// to be replaced by pre-calculated value from a surface area percentage
	// (which is provided as a model parameter).
	newCell->SetPoleRegionAngle(POLAR_REGION_ANGLE);

	// Add newly created cell to population
	AddCell(newCell);
}

void ModelCellsSpherical::InitCell(CellSpherical * newCell)
{
	newCell->AttachToModel(this);
	newCell->setState(Cell::State::StateQuiescent, false);
	newCell->SimObjId = no_cells++;
}

//! Adds an already created cell.  Useful when loading data.
void ModelCellsSpherical::AddCell(CellSpherical * newCell)
{
	UpdateParametersForCell(newCell);
	UpdateCellsStaining(newCell);
	newCell->AttachToModel(this);

	cells.push_back(newCell);
	cells2->add(newCell);

	mpArena->addObject(newCell->GLObject());

	//std::string outputFileName = "debugFile.txt";
	//std::ofstream FO;
	//FO.open(outputFileName.c_str(), std::ios::out | std::ios::app);
	//FO << "ModelCellsSpherical.cpp / AddCell(CellSpherical *)	add new cell/object transpartent: " << newCell->GLObject()->isTransparent() << std::endl;
	//FO.close();
}

//Add cells from mxf file
void ModelCellsSpherical::AddReadCells(std::string & filename)
{
	setlocale(LC_ALL, "C");

	//read mxf cell position
	std::fstream f;

	const int l = 512;
	char cstring[l];

	for (unsigned int i = 0; i < l; i++)
		cstring[i] = ' ';

	std::string s = filename.c_str();

	f.open(filename.c_str(), std::ios::in);

	//search postion
	bool search = 1;
	while (!f.eof() && search)
	{
		if (strncmp(cstring, "<CellList>", 10) == 0)
			search = 0;
		f.getline(cstring, sizeof(cstring));
	}

	//read position
	search = 1;
	while (!f.eof() && search)
	{
		if (strncmp(cstring, "</CellList>", 11) == 0)
			search = 0;
		else
		{
			char *pch;
			pch = strtok(cstring, " \t");
			pch = strtok(nullptr, " ");
			pch = strtok(nullptr, "\"");
			pch = strtok(nullptr, "\"");

			double tmp_position[3];

			int index = 0;

			while (pch != nullptr && index < 3)
			{

				tmp_position[index] = atof(pch);
				pch = strtok(nullptr, "\"");
				pch = strtok(nullptr, "\"");
				index++;
			}
			if (mUsePolarCells)
				AddPolarCell(tmp_position[0], tmp_position[1], tmp_position[2]);
			else
				AddCell(tmp_position[0], tmp_position[1], tmp_position[2]);

			cells.back()->setQuality(10, 10);
			cells.back()->setState(Cell::StateQuiescent);

			f.getline(cstring, sizeof(cstring));
		}
	}

	f.close();
}

void ModelCellsSpherical::AddLobuleShape(double radius, double height, double, double, double, double, double, bool)
{
	//Add lobule shape
	if (!is2D)
	{
		double d0[3][3] = { { 0,0,-height*0.5 },{ -radius*std::sqrt(3.)*0.5,-radius*0.5,-height*0.5 },{ 0,-radius,-height*0.5 } };
		AddBarrierTriangle(d0, 0);
		double u0[3][3] = { { 0,0,height*0.5 },{ 0,-radius,height*0.5 },{ -radius*std::sqrt(3.)*0.5,-radius*0.5,height*0.5 } };
		AddBarrierTriangle(u0, 0);
		double d1[3][3] = { { 0,0,-height*0.5 },{ -radius*std::sqrt(3.)*0.5,radius*0.5,-height*0.5 },{ -radius*std::sqrt(3.)*0.5,-radius*0.5,-height*0.5 } };
		AddBarrierTriangle(d1, 0);
		double u1[3][3] = { { 0,0,height*0.5 },{ -radius*std::sqrt(3.)*0.5,-radius*0.5,height*0.5 },{ -radius*std::sqrt(3.)*0.5,radius*0.5,height*0.5 } };
		AddBarrierTriangle(u1, 0);
		double d2[3][3] = { { 0,0,-height*0.5 },{ 0,radius,-height*0.5 },{ -radius*std::sqrt(3.)*0.5,radius*0.5,-height*0.5 } };
		AddBarrierTriangle(d2, 0);
		double u2[3][3] = { { 0,0,height*0.5 },{ -radius*std::sqrt(3.)*0.5,radius*0.5,height*0.5 },{ 0,radius,height*0.5 } };
		AddBarrierTriangle(u2, 0);
		double d3[3][3] = { { 0,0,-height*0.5 },{ radius*std::sqrt(3.)*0.5,radius*0.5,-height*0.5 },{ 0,radius,-height*0.5 } };
		AddBarrierTriangle(d3, 0);
		double u3[3][3] = { { 0,0,height*0.5 },{ 0,radius,height*0.5 },{ radius*std::sqrt(3.)*0.5,radius*0.5,height*0.5 } };
		AddBarrierTriangle(u3, 0);
		double d4[3][3] = { { 0,0,-height*0.5 },{ radius*std::sqrt(3.)*0.5,-radius*0.5,-height*0.5 },{ radius*std::sqrt(3.)*0.5,radius*0.5,-height*0.5 } };
		AddBarrierTriangle(d4, 0);
		double u4[3][3] = { { 0,0,height*0.5 },{ radius*std::sqrt(3.)*0.5,radius*0.5,height*0.5 },{ radius*std::sqrt(3.)*0.5,-radius*0.5,height*0.5 } };
		AddBarrierTriangle(u4, 0);
		double d5[3][3] = { { 0,0,-height*0.5 },{ 0,-radius,-height*0.5 } ,{ radius*std::sqrt(3.)*0.5,-radius*0.5,-height*0.5 } };
		AddBarrierTriangle(d5, 0);
		double u5[3][3] = { { 0,0,height*0.5 },{ radius*std::sqrt(3.)*0.5,-radius*0.5,height*0.5 },{ 0,-radius,height*0.5 } };
		AddBarrierTriangle(u5, 0);
	}

	double s0[3][3] = { { -radius*std::sqrt(3.)*0.5,-radius*0.5,-height*0.5 },{ 0,-radius,height*0.5 },{ 0,-radius,-height*0.5 } };
	AddBarrierTriangle(s0);
	double s1[3][3] = { { -radius*std::sqrt(3.)*0.5,-radius*0.5,-height*0.5 },{ -radius*std::sqrt(3.)*0.5,-radius*0.5,height*0.5 },{ 0,-radius,height*0.5 } };
	AddBarrierTriangle(s1);
	double s2[3][3] = { { -radius*std::sqrt(3.)*0.5,radius*0.5,-height*0.5 },{ -radius*std::sqrt(3.)*0.5,-radius*0.5,height*0.5 },{ -radius*std::sqrt(3.)*0.5,-radius*0.5,-height*0.5 } };
	AddBarrierTriangle(s2);
	double s3[3][3] = { { -radius*std::sqrt(3.)*0.5,radius*0.5,-height*0.5 },{ -radius*std::sqrt(3.)*0.5,radius*0.5,height*0.5 },{ -radius*std::sqrt(3.)*0.5,-radius*0.5,height*0.5 } };
	AddBarrierTriangle(s3);
	double s4[3][3] = { { 0,radius,-height*0.5 },{ -radius*std::sqrt(3.)*0.5,radius*0.5,height*0.5 },{ -radius*std::sqrt(3.)*0.5,radius*0.5,-height*0.5 } };
	AddBarrierTriangle(s4);
	double s5[3][3] = { { 0,radius,-height*0.5 },{ 0,radius,height*0.5 },{ -radius*std::sqrt(3.)*0.5,radius*0.5,height*0.5 } };
	AddBarrierTriangle(s5);
	double s6[3][3] = { { radius*std::sqrt(3.)*0.5,radius*0.5,-height*0.5 },{ 0,radius,height*0.5 },{ 0,radius,-height*0.5 } };
	AddBarrierTriangle(s6);
	double s7[3][3] = { { radius*std::sqrt(3.)*0.5,radius*0.5,-height*0.5 },{ radius*std::sqrt(3.)*0.5,radius*0.5,height*0.5 },{ 0,radius,height*0.5 } };
	AddBarrierTriangle(s7);
	double s8[3][3] = { { radius*std::sqrt(3.)*0.5,-radius*0.5,-height*0.5 },{ radius*std::sqrt(3.)*0.5,radius*0.5,height*0.5 },{ radius*std::sqrt(3.)*0.5,radius*0.5,-height*0.5 } };
	AddBarrierTriangle(s8);
	double s9[3][3] = { { radius*std::sqrt(3.)*0.5,-radius*0.5,-height*0.5 },{ radius*std::sqrt(3.)*0.5,-radius*0.5,height*0.5 },{ radius*std::sqrt(3.)*0.5,radius*0.5,height*0.5 } };
	AddBarrierTriangle(s9);
	double s10[3][3] = { { 0,-radius,-height*0.5 },{ radius*std::sqrt(3.)*0.5,-radius*0.5,height*0.5 },{ radius*std::sqrt(3.)*0.5,-radius*0.5,-height*0.5 } };
	AddBarrierTriangle(s10);
	double s11[3][3] = { { 0,-radius,-height*0.5 },{ 0,-radius,height*0.5 },{ radius*std::sqrt(3.)*0.5,-radius*0.5,height*0.5 } };
	AddBarrierTriangle(s11);
}

void ModelCellsSpherical::AddBarrierTriangle(double mPoints[][3], bool visible,
	double epsilon, double r, double g,
	double b, double a)
{
	ModelElementBarrierTriangle *mebt = new ModelElementBarrierTriangle(0, 0, 0);
	for (int i = 0; i < 3; i++)
		mebt->setPoint(mPoints[i][0], mPoints[i][1], mPoints[i][2], i);
	mebt->setNormalVector();
	mebt->setBoundingBox(epsilon);
	mebt->SetColor(r, g, b, a);
	mebt->freeSurfaceArea = 1;
	cells2->add(mebt);
	mBarriers.push_back(mebt);
	if (visible)
		this->mpArena->addObject(mebt->GLObject());
}

void ModelCellsSpherical::AddBloodVesselNetwork(GraphSphere * network)
{
	mpGraphBloodVesselNetwork = network;

	mpGraphBloodVesselNetwork->defaultPoissonRatio = defaultPoissonRatioSinusoids;
	mpGraphBloodVesselNetwork->defaultYoungModulus = defaultYoungModulusSinusoids;

	mpGraphBloodVesselNetwork->setBoundingBoxList(cells2);

	mpGraphBloodVesselNetwork->SampleGraph(0.15);

	for (unsigned int i = 0; i< mpGraphBloodVesselNetwork->mvNode.size(); i++)
		mpArena->addObject(mpGraphBloodVesselNetwork->mvNode[i]->GLObject());

	//calc dimension of lobule
	double shift = mpGraphBloodVesselNetwork->calcDimensions(mlobule_radius, mlobule_height);

	//shift height to center
	mpGraphBloodVesselNetwork->shift(-shift, 2);

	if (mReadCells)
	{
		for (unsigned int i = 0; i < cells.size(); i++)
		{
			cells[i]->position.z -= shift;
			cells[i]->boundingBox()->zmin -= shift;
			cells[i]->boundingBox()->zmax -= shift;
		}
	}

	//set maximal values in each dimensions
	mpGraphBloodVesselNetwork->setBoundingBox();

	//connect nodes
	mpGraphBloodVesselNetwork->ConnectNodes();
}

void ModelCellsSpherical::AddBloodVesselNetwork(std::string & filename, int filetype)
{
	GraphSphere * network = new GraphSphere();

	network->setBoundingBoxList(cells2);

	network->read(filename.c_str(), filetype);

	AddBloodVesselNetwork(network);
}

void ModelCellsSpherical::AddShape(double radius, double height, double r, double g, double b, double a, double epsilon, bool visible)
{
	if (is2D)
		mlobule_height = 1;

	if (mLobuleShape == ModelCellsSpherical::Hexagonal)
		AddLobuleShape(mlobule_radius, mlobule_height, 1, 1, 0, 0.5, 1);
	else if (mLobuleShape == ModelCellsSpherical::Quadric)
		AddQuadricShape(mlobule_radius, mlobule_height, 1, 1, 0, 0.5, 1);
}

void ModelCellsSpherical::AddQuadricShape(double radius, double height, double r, double g, double b, double a, double epsilon, bool visible)
{
	//Add lobule shape
	if (!is2D)
	{
		double d0[3][3] = { { -radius,-radius,-height*0.5 },{ radius,radius,-height*0.5 },{ -radius,radius,-height*0.5 } };
		AddBarrierTriangle(d0, 0);
		double d1[3][3] = { { -radius,-radius,-height*0.5 },{ radius,-radius,-height*0.5 },{ radius,radius,-height*0.5 } };
		AddBarrierTriangle(d1, 0);
		double u0[3][3] = { { -radius,-radius,height*0.5 },{ -radius,radius,height*0.5 },{ radius,radius,height*0.5 } };
		AddBarrierTriangle(u0, 0);
		double u1[3][3] = { { -radius,-radius,height*0.5 },{ radius,radius,height*0.5 },{ radius,-radius,height*0.5 } };
		AddBarrierTriangle(u1, 0);
	}

	double s0[3][3] = { { -radius,-radius,-height*0.5 },{ -radius,radius,height*0.5 },{ -radius,-radius,height*0.5 } };
	AddBarrierTriangle(s0);
	double s1[3][3] = { { -radius,radius,-height*0.5 },{ -radius,radius,height*0.5 },{ -radius,-radius,-height*0.5 } };
	AddBarrierTriangle(s1);
	double s2[3][3] = { { radius,-radius,-height*0.5 },{ radius,-radius,height*0.5 },{ radius,radius,height*0.5 } };
	AddBarrierTriangle(s2);
	double s3[3][3] = { { radius,radius,-height*0.5 },{ radius,-radius,-height*0.5 },{ radius,radius,height*0.5 } };
	AddBarrierTriangle(s3);
	double s4[3][3] = { { -radius,radius,-height*0.5 },{ radius,radius,height*0.5 },{ -radius,radius,height*0.5 } };
	AddBarrierTriangle(s4);
	double s5[3][3] = { { radius,radius,height*0.5 },{ -radius,radius,-height*0.5 },{ radius,radius,-height*0.5 } };
	AddBarrierTriangle(s5);
	double s6[3][3] = { { -radius,-radius,-height*0.5 },{ -radius,-radius,height*0.5 },{ radius,-radius,height*0.5 } };
	AddBarrierTriangle(s6);
	double s7[3][3] = { { radius,-radius,height*0.5 },{ radius,-radius,-height*0.5 },{ -radius,-radius,-height*0.5 } };
	AddBarrierTriangle(s7);
}

void ModelCellsSpherical::InitForces()
{
	// For all cells
	//std::string outputFileName = "debugFile.txt";
	//std::ofstream FO;
	//FO.open(outputFileName.c_str(), std::ios::out | std::ios::app);
	//FO << "ModelCellsSpherical.cpp / InitForces()	starts" << endl;
	for (unsigned int i = 0; i < cells2->EndIndex(); i++)
	{
		// dead cells might remain in cells2 as nullptr pointers!
		if (ModelElement *element = cells2->element(i))
			element->Reset();
	}
	//FO << "ModelCellsSpherical.cpp / InitForces()	after element reset." << endl;

	mpInteractionFrictionMatrix->Reset();
	//FO << "ModelCellsSpherical.cpp / InitForces()	after interactionFrictionMatrix reset." << endl;
	mProblemAllocationSize = mSolver->setProblemSize();
	//FO << "ModelCellsSpherical.cpp / InitForces()	after mSolver setProblemSize." << endl;
	mpVelocities = mSolver->getVelocity();
	//FO << "ModelCellsSpherical.cpp / InitForces()	after mSolver getVelocity." << endl;

	if (mpInteractionJKR)
	{
		if (cells2->size() > mJKRElementDoneSize)
		{
			std::size_t difference = cells2->size() - mJKRElementDoneSize;
			difference = ((std::size_t) std::ceil(1.25 * difference / 1024)) << 10;

			mJKRElementDoneSize += difference;
			// for JKR only
			mpJKRElementDone =
				(bool *)realloc(mpJKRElementDone, mJKRElementDoneSize * sizeof(bool));
		}

		memset(mpJKRElementDone, 0, mJKRElementDoneSize * sizeof(bool));
	}
	//FO << "ModelCellsSpherical.cpp / InitForces()	after interactionJKR." << endl;
	//FO.close();
}

void ModelCellsSpherical::UpdateInteractions()
{
	// For all cells
	//std::string outputFileName = "debugFile.txt";
	//std::ofstream FO;
	//FO.open(outputFileName.c_str(), std::ios::out | std::ios::app);
	//FO << "ModelCellsSpherical.cpp / UpdateInteractions()	starts" << endl;
	cells2->update();
	//FO << "ModelCellsSpherical.cpp / UpdateInteractions()	after cells2 update" << endl;

	BoundingBoxList::iterator cellIterator = cells2->begin();

	//FO << "ModelCellsSpherical.cpp / UpdateInteractions()	while starts" << endl;
	while (cellIterator != cells2->end())
	{
		//cout << "while" << endl;
		ModelElement *object_1 = *cellIterator;

		//cout << "interating things" << endl;
		//cout << "cell2:" << cells2 << endl;
		//cout << "cell2_size:" << cells2->size() << endl;
		//cout << "object:" << object_1 << endl;
		//FO << "ModelCellsSpherical.cpp / UpdateInteractions()	object " << object_1 << endl;
		CSListContainer< unsigned long > & interactingElements = object_1->mIntersectingList;

		// prepare interactingElements for JKR
		size_t counter = 0;
		size_t switchToNonEstablishedContacts = object_1->mEstablishedContacts.size();

		//cout << "process" << endl;
		//FO << "ModelCellsSpherical.cpp / UpdateInteractions()		process" << endl;
		if (mpInteractionJKR)
		{
			//mpInteractionJKR->time = biolink->getTimeInDays(time);
			// see if elements in the newly interactingElements are in the establishedContacts
			//  push_back into establishedContacts only if not.
			unsigned long * iterator;
			for (iterator = interactingElements.begin();
				iterator != interactingElements.end();
				++iterator)
			{
				//cout << "iterator:" << *iterator << endl;
				//FO << "ModelCellsSpherical.cpp / UpdateInteractions()		iterator: " << *iterator << endl;
				if (!mpJKRElementDone[*iterator])
				{
					if (object_1->mEstablishedContacts.find(*iterator) == object_1->mEstablishedContacts.end())
						object_1->mEstablishedContacts.push_back(*iterator);
				}
			}
			// the contents of establishedContacts takes over the function of interactingElements
			interactingElements.swap(object_1->mEstablishedContacts);

			// the first elements in interactingElements are from object_1->mEstablishedcontacts,
			// switch to false, when switchToNonEstablishedContacts is reached.
			mpInteractionJKR->mContactEstablished = true;
		}

		unsigned long * otherIterator = interactingElements.begin();

		while (otherIterator != interactingElements.end())
		{
			ModelElement * object_2 = cells2->element(*otherIterator);

			if (mpInteractionHertz)
			{
				(*mpInteractionHertz)(object_1, object_2);

				// if there's no contact area, there's no Hertz interaction,
				// next please!
				if (mpInteractionHertz->mContactArea == 0)
				{
					++otherIterator;
					continue;
				}
			}
			else  // JKR!
			{
				if (counter++ == switchToNonEstablishedContacts)
					mpInteractionJKR->mContactEstablished = false;

				if (mpJKRElementDone[*otherIterator])
				{
					++otherIterator;
					continue;
				}

				(*mpInteractionJKR)(object_1, object_2);

				if (mpInteractionJKR->mContactArea == 0)
				{
					++otherIterator;
					continue;
				}
			}

			// if ( useDetailedCellCellFriction )
			(*mpInteractionFrictionMatrix)(object_1, object_2);

			++otherIterator;
		}

		if (mpInteractionJKR)
			mpJKRElementDone[object_1->mGlobalIndex] = true;

		if (object_1->mType != ModelElement::TypeBarrierTriangle)
		{
			ModelElementSphere* cell1 = static_cast<ModelElementSphere*>(object_1);
			cell1->freeSurfaceArea = 4 * M_PI * cell1->mRadius * cell1->mRadius
				- cell1->surfaceContactArea;

			if (cell1->freeSurfaceArea <0)
				cell1->freeSurfaceArea = 0;
		}

		++cellIterator;
	}
	//FO << "ModelCellsSpherical.cpp / UpdateInteractions()	after while" << endl;
	//FO.close();
}

void ModelCellsSpherical::setFrictionMatrix()
{
	// bug file
	//std::string outputFileName = "debugFile.txt";
	//std::ofstream FO;
	//FO.open(outputFileName.c_str(), std::ios::out | std::ios::app);
	//FO << "ModelCellsSpherical.cpp / setFrictionMatrix()	step " << step << endl;
	//
	mpFrictionMatrices = mpInteractionFrictionMatrix->mpFrictionMatrices;
	//FO << "ModelCellsSpherical.cpp / setFrictionMatrix()	after mpFrictionMatrices." << endl;
	mSolver->setFrictionMatrix(mpFrictionMatrices);
	//FO << "ModelCellsSpherical.cpp / setFrictionMatrix()	after mSolver." << endl;
	//FO.close();
}

void ModelCellsSpherical::UpdateDishContactForce(CellSpherical * cell)
{
	assert(!is2D);
	// we assume here that the petry dish is in the xz-plane with z = 0
	Vector3f closest_point = { cell->position.x, cell->position.y, 0 };

	double distance = CSModelTools::GetDistance3D(closest_point, cell->position);

	double E_eff = 1.0
		/ ((1.0 - cell->poissonRatio * cell->poissonRatio)
			/ cell->youngModulus
			+ (1.0 - 0.25)
			/ 1e9);

	// shouldnt this be per area of the contact?
	double adhesion_constant = 3 * singleBondEnergy * adhesionDensity;

	double dcrit = std::pow(3.0 * cell->mRadius *
		std::pow(M_PI * adhesion_constant / (8.0 * E_eff), 2.0), 0.33333333);

	double delta = cell->mRadius - distance;

	if (delta < -dcrit)
	{
		std::cout << "Cell(" << cell->mGlobalIndex << "): NO CONTACT TO SUBSTRATE!"
			<< "\nDistance is:" << distance << " z:" << cell->position.z
			<< std::endl;
		return;
	}

	double mForce, mContactArea;

	//update JKR contact force and area.
	CSModelTools::GetJKRForce(E_eff, cell->mRadius, delta, adhesion_constant,
		mForce, mContactArea, dcrit);

	double force_abs = fabs(mForce);
	cell->lastForce += mForce;
	cell->accumulatedForceAbsolute += force_abs;
	cell->accumulatedPressure += force_abs / mContactArea;

	cell->directedForce.Add(0, 0, mForce);
	cell->surfaceContactArea += mContactArea;
}

void ModelCellsSpherical::DumbbellMetropolis()
{
	// bug file
	//std::string outputFileName = "debugFile.txt";
	//std::ofstream FO;
	//FO.open(outputFileName.c_str(), std::ios::out | std::ios::app);
	//FO << "ModelCellsSpherical.cpp / DumbbellMetropolis()	step " << step << ", number of cells: " << cells.size() << endl;
	const double referenceEnergy = 1.;

	// using the same condition as for the timestep: displacements must be
	// smaller than maxDisplacementRadius*initialRadius = 0.1*0.5 the
	// displacement due to a rotation alpha with maximum dumbbell extents
	// 2*divisionDistance is D = 2*divisionDistance * sin(alpha/2)
	// Setting D to maxDisplacementRadius*initialRadius gives us the maximum angle:
	// maxAngle = 2*arcsin (maxDisplacementRadius*initialRadius/(2*divisionDistance))
	//          = 2*arcsin (0.05/0.8) ~= 0.1250815
	const double maxAngle = 2 * std::asin(mMaximumDisplacementCellRadius*defaultInitialCellRadius / (2 * defaultDivisionDistance));
	//FO << "ModelCellsSpherical.cpp / DumbbellMetropolis()	maxAngle: " << maxAngle << endl;

	// memory for saving the rotations, in case they have to be reverted:
	std::vector< std::pair< CellSpherical*, std::pair<Vector3f, double> > > rotationRecord;
	//std::vector< std::pair< CellSpherical*, std::pair<Vector3f, Vector3f> > > rotationRecord;
	bool accepted = false;

	double energy = 0.;
	for (CellSpherical * cell : cells)
	{
		if (cell->mDaughterCell) // quickest way to get only one of the dumbbell partners.
		{
			//std::cout << "Compute dumbbell energies!"
			//    << " (" << cell->mpDivisionPartner->SimObjId << ", "
			//    << cell->SimObjId << ")."
			//    << std::endl;

			//FO << "ModelCellsSpherical.cpp / DumbbellMetropolis()	cell has mDaughter." << endl;
			//FO << "ModelCellsSpherical.cpp / DumbbellMetropolis()	cell start on oldInteractionEnergyCellDivisionPartner." << endl;
			double oldInteractionEnergyCellDivisionPartner = GetInteractionEnergy(cell->mpDivisionPartner);
			//FO << "ModelCellsSpherical.cpp / DumbbellMetropolis()	cell oldInteractionEnergyCellDivisionPartner: " << oldInteractionEnergyCellDivisionPartner << endl;
			//FO << "ModelCellsSpherical.cpp / DumbbellMetropolis()	cell start on oldInteractionEnergyCell." << endl;
			double oldInteractionEnergyCell = GetInteractionEnergy(cell);
			//FO << "ModelCellsSpherical.cpp / DumbbellMetropolis()	cell oldInteractionEnergyCell: " << oldInteractionEnergyCell << endl;
			double oldInteractionEnergy = oldInteractionEnergyCellDivisionPartner + oldInteractionEnergyCell;

			// a random rotation of a vector consists of a rotational
			// axis perpendicular to the vector and an angle.  Both axis
			// and angle have to be picked randomly, while the
			// distribution for the axis needs to be isotropic on the
			// one-circle in the plane perpendicular to the vector.
			Vector3f axis;
			if (is2D)
			{
				short sign = (mRandom.GetRandomUniform11() >= 0) ? 1 : -1;
				// if the dumbbell has to stay in the xy-plane, the
				// rotation axis is the z-axis:
				axis.Set(0., 0., sign*1.);
			}
			else
			{
				Vector3f dumbbellAxis = cell->position - cell->mpDivisionPartner->position;
				double x = (dumbbellAxis[0]) ? dumbbellAxis[0] : 1.;
				axis.Set(-dumbbellAxis[1] / x, 1., 0);
				axis.Normalize();

				axis = Rotate(axis, dumbbellAxis, 2 * M_PI * mRandom.GetRandomUniform01());
				axis = Normalize(axis);
			}

			// sampling from a distribution between 0 and maxAngle
			// uniformly distributed on a spherical cap with an opening
			// angle of maxAngle
			// (cf. CellSphericalPolar::probeRandomRotation()):
			// Using the density function 2*M_PI*r*sin(theta) (r=1,
			// 0<=theta<thetaMax) we have a distribution of F(theta) =
			// (1-cos(theta))/(1-cos(thetaMax)).  So, we can generate a
			// random number x uniformly distributed over [0,1) and
			// calculate our angle with the inverse function theta =
			// G(x) = acos( 1 - x * (1 - cos(thetaMax)) );
			double angle = std::acos(1 - mRandom.GetRandomUniform01() * (1 - std::cos(maxAngle)));
			//FO << "ModelCellsSpherical.cpp / DumbbellMetropolis()	cell angle: " << angle << endl;

			//rotationRecord.push_back(make_pair(cell, make_pair(cell->position, cell->mpDivisionPartner->position)));
			cell->RotateDumbbell(axis, angle);
			//FO << "ModelCellsSpherical.cpp / DumbbellMetropolis()	cell after rotateDumbbell." << endl;

			/*if (isnan(cell->position))
			{
				std::cout << "cell:" << cell->position << std::endl;
				std::cout << "cell:" << cell->mpDivisionPartner->position << std::endl;
			}*/

			std::pair<Vector3f, double> rot(axis, angle);
			rotationRecord.push_back(std::pair<CellSpherical*, std::pair<Vector3f, double> >(cell, rot));

			//FO << "ModelCellsSpherical.cpp / DumbbellMetropolis()	cell start on newEnergyCell." << endl;
			double energyCell = GetInteractionEnergy(cell);
			//FO << "ModelCellsSpherical.cpp / DumbbellMetropolis()	newEnergyCell: " << energyCell << endl;
			//FO << "ModelCellsSpherical.cpp / DumbbellMetropolis()	cell start on newEnergyCellDivisionPartner." << endl;
			double energyCellDivisionPartner = GetInteractionEnergy(cell->mpDivisionPartner);
			//FO << "ModelCellsSpherical.cpp / DumbbellMetropolis()	newEnergyCellDivisionPartner: " << energyCellDivisionPartner << endl;
			double newEnergy = energyCell + energyCellDivisionPartner;
			//FO << "ModelCellsSpherical.cpp / DumbbellMetropolis()	cell newEnergy: " << newEnergy << endl;
			//energy += oldInteractionEnergy
			//    - ( GetInteractionEnergy(cell) + GetInteractionEnergy(cell->mpDivisionPartner) );
			energy += (newEnergy - oldInteractionEnergy);

			//std::cout << "\tE:" << energy << " angle:" << angle << std::endl;
		}
		else
		{
			//FO << "ModelCellsSpherical.cpp / DumbbellMetropolis()	cell has no mDaughter." << endl;
		}
	}

	assert(!std::isnan(energy));
	assert(std::isfinite(energy));

	//FO << "ModelCellsSpherical.cpp / DumbbellMetropolis()	after assertEnergy." << endl;

	if (energy > 0)
	{
		if (mRandom.GetRandomUniform01() >= std::exp(-energy / referenceEnergy))
		{
			for (auto rot : rotationRecord) // reverse rotations;
			{
				//rot.first->position = rot.second.first;
				//rot.first->mpDivisionPartner->position = rot.second.second;
				Vector3f axis = -1. * rot.second.first;
				rot.first->RotateDumbbell(axis, rot.second.second);
			}
		}
		//else
		//    std::cout << "Rotation accepted" << std::endl;
	}
	//FO << "ModelCellsSpherical.cpp / DumbbellMetropolis()	after ifEnergy." << endl;
	//FO.close();
}

double ModelCellsSpherical::GetInteractionEnergy(CellSpherical * cell)
{
	// bug file
	//std::string outputFileName = "debugFile.txt";
	//std::ofstream FO;
	//FO.open(outputFileName.c_str(), std::ios::out | std::ios::app);
	//FO << "ModelCellsSpherical.cpp / GetInteractionEnergy()	step " << step << endl;
	//
	CSListContainer<ModelElement *> possibleContacts = cells2->GetContainedElements(*cell->getBoundingBox());

	//FO << "ModelCellsSpherical.cpp / GetInteractionEnergy()	after possibleContacts." << endl;

	double energy = 0.0;

	if (mpInteractionHertz != nullptr)
	{
		//FO << "ModelCellsSpherical.cpp / GetInteractionEnergy()	InteractionHertz" << endl;
		CSInteractionHertzEnergy potential;
		potential.mpSingleBondEnergy = mpInteractionHertz->mpSingleBondEnergy;
		potential.mpAdhesionDensity = mpInteractionHertz->mpAdhesionDensity;

		for (ModelElement** pContact = possibleContacts.begin(); pContact != possibleContacts.end(); ++pContact)
		{
			if (*pContact == static_cast<ModelElement*>(cell)) continue;
			potential(cell, *pContact);
			energy += potential.mEnergy;
		}
	}
	else if (mpInteractionJKR != nullptr)
	{
		//FO << "ModelCellsSpherical.cpp / GetInteractionEnergy()	InteractionJKR" << endl;
		CSInteractionJKREnergy potential;
		potential.mpSingleBondEnergy = mpInteractionJKR->mpSingleBondEnergy;
		potential.mpAdhesionDensity = mpInteractionJKR->mpAdhesionDensity;

		//FO << "ModelCellsSpherical.cpp / GetInteractionEnergy()	after potentialEnergy." << endl;

		for (ModelElement** pContact = possibleContacts.begin(); pContact != possibleContacts.end(); ++pContact)
		{
			if (*pContact == static_cast<ModelElement*>(cell)) continue;

			// potential has to know if the two cells are in contact
			// This should be mEstablishedContacts then contact is established
			// and something can happen
			auto result = cell->mEstablishedContacts.find((*pContact)->mGlobalIndex);
			//FO << "ModelCellsSpherical.cpp / GetInteractionEnergy()		after result." << endl;
			potential.mContactEstablished = (result != cell->mEstablishedContacts.end());
			//FO << "ModelCellsSpherical.cpp / GetInteractionEnergy()		after contactEstablished." << endl;

			potential(cell, *pContact);
			//FO << "ModelCellsSpherical.cpp / GetInteractionEnergy()		after potential." << endl;
			energy += potential.mEnergy;
			//FO << "ModelCellsSpherical.cpp / GetInteractionEnergy()		after potentialMEnergy: " << potential.mEnergy << endl;
		}
		//FO << "ModelCellsSpherical.cpp / GetInteractionEnergy()	after forContactLoop." << endl;
	}
	//FO << "ModelCellsSpherical.cpp / GetInteractionEnergy()	done and return energy: " << energy << endl;
	//FO.close();

	return energy;
}

void ModelCellsSpherical::UpdateCellFriction(CellSpherical * cell)
{
	// Add cell-cell friction
	//cell->frictionCoefficient += gammaCellsParallel * cell->surfaceContactArea;

	// Set surface area to the remaining surface area (not in contact with other cells)
	// Full surface of (spherical) cell is 4*pi*R^2
	if (is2D)
	{
		double substrate_contact = CSModelTools::GetContactAreaHertz(cell->mRadius, cell->position.z);
		cell->surfaceContactArea += substrate_contact;
		cell->frictionCoefficient += gammaECM * substrate_contact;
	}

	// set to zero for now
	cell->surfaceContactArea = 4 * M_PI * cell->mRadius * cell->mRadius - cell->surfaceContactArea;
	//cell->surfaceContactArea = 4 * M_PI * 0.5 * 0.5;
	// we will just always assume that we have no contacts
	//- cell->surfaceContactArea;

	// Add cell-ecm friction
	if (cell->surfaceContactArea > 0) // If there is still some free surface available
	{
		cell->frictionCoefficient += gammaMedium * cell->surfaceContactArea;
	}

	assert(cell->frictionCoefficient >= 0.0);
}

void ModelCellsSpherical::UpdateForcesLangevin(CellSpherical *cell)
{
	// Update sigma and langevin forces (after friction is completely known)
	double s;

	// For all cells
	// Earlier version: s^2 / timestep

	double absoluteForce;

	Vector3f unitVector(0, 0, 0);
	if (is2D)
		mRandom.GetRandomUnitVector(&unitVector.x, &unitVector.y);
	else
		mRandom.GetRandomUnitVector(&unitVector.x, &unitVector.y, &unitVector.z);

	//double friction = 0.1;
	s = cell->frictionCoefficient *
		sqrt(2. * dimension * defaultDiffusionConstantCells / timeStep);
	//s = 0.1 * sqrt(6. * defaultDiffusionConstantCells / timeStep);

	unitVector *= s;

	//s = cell->frictionCoefficient * sqrt(2 * defaultDiffusionConstantCells / timeStep);
	//cell->mLangevinForce.x = mRandom.GetRandomGauss(0.0, s);
	//cell->mLangevinForce.y = mRandom.GetRandomGauss(0.0, s);

	//if (is2D)
	//{
	//    cell->mLangevinForce.z = 0;
	//    absoluteForce = sqrt(cell->mLangevinForce.x*cell->mLangevinForce.x
	//                        +cell->mLangevinForce.y*cell->mLangevinForce.y );
	//}
	//else
	//{
	//    cell->mLangevinForce.z = mRandom.GetRandomGauss(0.0, s);
	//    absoluteForce = cell->mLangevinForce.Norm();
	//}

	cell->mLangevinForce = unitVector;
	absoluteForce = unitVector.Norm();

	if (cell->mpDivisionPartner)
	{
		cell->directedForce.Add(cell->mLangevinForce.x / 2,
			cell->mLangevinForce.y / 2,
			cell->mLangevinForce.z / 2);

		cell->lastForceAbsolute += absoluteForce / 2.;

		cell->mpDivisionPartner->directedForce.Add(cell->mLangevinForce.x / 2,
			cell->mLangevinForce.y / 2,
			cell->mLangevinForce.z / 2);

		cell->mpDivisionPartner->lastForceAbsolute += absoluteForce / 2.;
	}
	else
	{
		cell->directedForce.Add(cell->mLangevinForce.x,
			cell->mLangevinForce.y,
			cell->mLangevinForce.z);

		cell->lastForceAbsolute += absoluteForce;
	}

	//std::cout << "F=" << cell->mLangevinForce
	//    << " F2=" << unitVector
	//     << " ||F||=" << absoluteForce
	//     << " ||F2||=" << unitVector.Norm()
	//     << " s2=" << s2
	//     << " s=" << s
	//     << " D=" << defaultDiffusionConstantCells
	//     << " dt=" << timeStep
	//     << " gamma=" << cell->frictionCoefficient
	//     << "\n";

	//cell->lastForceAbsolute += absoluteForce;

	// Review?  how do the langevin forces contribute to pressure?  Omitted for now.
}

//end Forces

//this method is only used when one click on "apply parameters to all existing cells" in the GUI
void ModelCellsSpherical::UpdateParametersForAllCells()
{
	UpdateParametersFromBioLink();

	for (unsigned int i = 0; i<cells.size(); i++)
	{
		UpdateParametersForCell(cells[i]);
	}
}

void ModelCellsSpherical::UpdateParametersForCell(CellSpherical *cell)
{
	// t1m:  divisionRadius and initialRadius have to be initialized after a
	//       data read.  Will override any individual setting.

#pragma region Size and growth specifics

	cell->divisionRadius = defaultDivisionCellRadius;
	cell->initialRadius = defaultInitialCellRadius;
	cell->defaultCellCycleTime = defaultCellCycleTime;
	cell->defaultCellCycleTimeStandardDeviation = defaultCellCycleTimeStandardDeviation;
	cell->defaultDivisionDistance = defaultDivisionDistance;
	cell->UpdateDeltaRadius(timeStep);

#pragma endregion

#pragma region Biophysics

	cell->youngModulus = defaultYoungModulusCells;
	cell->poissonRatio = defaultPoissonRatioCells;
	cell->mQuiescentMin = mQuiescentMin;
	cell->mQuiescentMax = mQuiescentMax;
	cell->mUseDumbbell = mUseDumbbell;
	cell->mForceThreshold = biolink->getForceDimensionless(mForceThreshold);
	cell->SetStateModel(StateMechanisms(mStateChoice));
	cell->mOverlapThreshold = overlap_threshold;
	cell->mOverlapThresholdMin = overlap_threshold_min;
	cell->mDensityCritical = density_threshold;

#pragma endregion

#pragma region Dumbbell specifics

	if (mUseDumbbell)
	{
		cell->mDumbbellPeriod = mDefaultDumbbellPeriod;
		if (cell->mpDivisionPartner)
			cell->mElongationFactor = defaultDivisionDistance / cell->mDumbbellPeriod;
	}

#pragma endregion
}

//
// TODO a lot of these functions could be autogenerated
//
void ModelCellsSpherical::printParameters() const
{
	std::stringstream ss;

	ss << "Parameters";
	ss << "\nLength Scale " << biolink->length_scale;
	ss << "\nCell Diameter " << biolink->diameter_bioCells;
	ss << "\ndefaultCellCycleTime " << defaultCellCycleTime;
	ss << "\ndefaultCellCycleTimeSD " << defaultCellCycleTimeStandardDeviation;
	ss << "\nQuiescentMin " << mQuiescentMin;
	ss << "\nQuiescentMax " << mQuiescentMax;
	ss << "\nYoung " << defaultYoungModulusCells;
	ss << "\nPoisson " << defaultPoissonRatioCells;
	ss << "\nYoungSinusoid " << defaultYoungModulusSinusoids;
	ss << "\nPoissonSinusoid " << defaultPoissonRatioSinusoids;
	ss << "\nInitialRadius " << defaultInitialCellRadius;
	ss << "\nInitialDividideRadius " << defaultDivisionCellRadius;
	ss << "\nDefaultDivisionDistance " << defaultDivisionDistance;
	ss << "\ndefaultDiffusion " << defaultDiffusionConstantCells;
	ss << "\nsingleBondEnergy " << singleBondEnergy;
	ss << "\nadhesionDensity " << adhesionDensity;
	ss << "\ngammaCellPerpendicular " << gammaCellsPerpendicular;
	ss << "\ngammaCellParallel " << gammaCellsParallel;
	ss << "\ngammaECM " << gammaECM;
	ss << "\ngammaMedium " << gammaMedium;
	ss << "\nKillRate " << mKillRate;
	ss << "\nStateModelChoice " << mStateChoice;
	ss << "\nKillInterval " << mKillInterval;
	ss << "\nNextKillTime " << mNextKillTime;
	ss << "\nKillPeriod " << mKillPeriod;
	ss << "\nKillEnd " << mKillEnd;
	ss << "\nOverlapThreshold " << overlap_threshold;
	ss << "\nOverlapThresholdMin " << overlap_threshold_min;
	ss << std::boolalpha << "\nis2D " << is2D;
	ss << std::boolalpha << "\nuseCCFriction " << mUseCCFriction;
	ss << "\n";
	std::cout << ss.rdbuf();
}

void ModelCellsSpherical::UpdateParametersFromBioLink()
{
	std::cout << "Biolink diameter=" << biolink->diameter_bioCells << std::endl;
	std::cout << "Biolink length_scale=" << biolink->length_scale << std::endl;

	// Cycle time
	defaultCellCycleTime
		= biolink->scaleBiologyToInternal(biolink->cycletime_bio,
			BiologyLink::ScaleTime);

	// Cycle time deviation (Gaussian distribution)
	defaultCellCycleTimeStandardDeviation
		= biolink->scaleBiologyToInternal(biolink->cycletime_stddev_bio,
			BiologyLink::ScaleTime);

	//min and max pressure for proliferation
	mQuiescentMin
		= biolink->scaleBiologyToInternal(biolink->quiescentMinPressure_bio,
			BiologyLink::ScalePressure);
	mQuiescentMax
		= biolink->scaleBiologyToInternal(biolink->quiescentMaxPressure_bio,
			BiologyLink::ScalePressure);

	// remark: the defaultcelldiameter don't need to be update, it is always
	// 0.5 by definition (we set the characteristic length scale so that cell diameter is 1)

	// Young modulus
	defaultYoungModulusCells = biolink->getYoungModuleDimensionless(biolink->youngModulus_bioCells);

	// Poisson ratio
	defaultPoissonRatioCells = biolink->poisson_ratio_bioCells;

	// Young modulus
	defaultYoungModulusSinusoids = biolink->getYoungModuleDimensionless(biolink->youngModulus_bioSinusoids);

	// Poisson ratio
	defaultPoissonRatioSinusoids = biolink->poisson_ratio_bioSinusoids;

	// Cell radius
	//defaultInitialCellRadius =
	//    biolink->scaleBiologyToInternal(biolink->diameter_bioCells, BiologyLink::ScaleLength);
	//defaultInitialCellRadius /= 2.0;

	defaultInitialCellRadius = 0.7;
	std::cout << "defaultInitialCellRadius=" << defaultInitialCellRadius << std::endl;

	// double volume then the cell is ready for division
	defaultDivisionCellRadius = std::pow(2, 0.33333333) * defaultInitialCellRadius;
	defaultDivisionDistance = 0.9 * defaultInitialCellRadius;

	std::cout << "defaultDivisionCellRadius=" << defaultDivisionCellRadius << std::endl;

	// Diffusion constant
	defaultDiffusionConstantCells
		= biolink->scaleBiologyToInternal(biolink->diffusion_constant_cells_bio,
			BiologyLink::ScaleDiffusivity);

	// Get single bond energy
	singleBondEnergy
		= biolink->scaleBiologyToInternal(biolink->single_bond_energy_bio,
			BiologyLink::ScaleEnergy);

	// Set cell-cell adhesion density
	adhesionDensity
		= biolink->scaleBiologyToInternal(biolink->adhesion_to_cells_bio,
			BiologyLink::ScaleAdhesionDensity);

	// Cell-cell friction coefficient
	gammaCellsPerpendicular
		= biolink->scaleBiologyToInternal(biolink->cell_gamma_perpendicular_bio,
			BiologyLink::ScaleFrictionCoefficient);

	gammaCellsParallel
		= biolink->scaleBiologyToInternal(biolink->cell_gamma_parallel_bio,
			BiologyLink::ScaleFrictionCoefficient);

	// Cell-ECM friction coefficient
	gammaECM
		= biolink->scaleBiologyToInternal(biolink->ecm_gamma_bio,
			BiologyLink::ScaleFrictionCoefficient);

	gammaMedium
		= biolink->scaleBiologyToInternal(biolink->medium_gamma_bio,
			BiologyLink::ScaleFrictionCoefficient);

	//set values depend on other variables
}

void ModelCellsSpherical::RegisterParameters()
{
	if (mpParameters)
		return;

	mpParameters = new CSParameterContext(name);

	CSParameterChoice * choiceContactModel = new CSParameterChoice(mContactModelNames, 0);

	mpParameters->setParameter("Contact Model", CSParameter::Choice, choiceContactModel, "");

	// length scale
	mpParameters->addParameter("Length Scale", CSParameter::Double, &biolink->length_scale, "m");

	// Default cell diameter
	mpParameters->addParameter("Cell Diameter", CSParameter::Double, &biolink->diameter_bioCells, "m");

	// cycle time
	mpParameters->addParameter("Cell Cycle Time", CSParameter::Double, &biolink->cycletime_bio, "s");

	// std deviation of cycle times
	mpParameters->addParameter("Cycle Time SD", CSParameter::Double, &biolink->cycletime_stddev_bio, "s");

	// If to use a dumbbell of two spheres for a period of time before the cell cleaves into two
	mpParameters->addParameter("Use Dumbbell Division", CSParameter::Bool, &mUseDumbbell, "");

	// dumbbell time
	mpParameters->addParameter("Dumbbell time", CSParameter::Double, &this->mDefaultDumbbellPeriod, "s");

	// Young Modulus Cells
	mpParameters->addParameter("Young Modulus Cells", CSParameter::Double, &biolink->youngModulus_bioCells, "Pa");

	// Poisson Ratio Cells
	mpParameters->addParameter("Poisson Ratio Cells", CSParameter::Double, &biolink->poisson_ratio_bioCells, "");

	// Young Modulus Sinusoids
	mpParameters->addParameter("Young Modulus Sinusoids", CSParameter::Double, &biolink->youngModulus_bioSinusoids, "Pa");

	// Poisson Ratio Sinusoids
	mpParameters->addParameter("Poisson Ratio Sinusoids", CSParameter::Double, &biolink->poisson_ratio_bioSinusoids, "");

	// diffusion constant for cells
	mpParameters->addParameter("Diffusion Constant", CSParameter::Double, &biolink->diffusion_constant_cells_bio, "m^2/s");

	// cell-cell adhesion
	mpParameters->addParameter("Cell-Cell Adhesion", CSParameter::Double, &biolink->adhesion_to_cells_bio, "m^-2");

	// If cells should be allocated as CellSphericalPolar or simply CellSpherical objects
	mpParameters->addParameter("Use Cell Polarity", CSParameter::Bool, &mUsePolarCells, "");

	// cell-cell gamma for shear friction
	mpParameters->addParameter("Cell-Cell gamma for shear friction", CSParameter::Double, &biolink->cell_gamma_perpendicular_bio, "Ns/m^3");

	// cell-cell gamma for collision friction
	mpParameters->addParameter("Cell-Cell gamma for normal friction", CSParameter::Double, &biolink->cell_gamma_parallel_bio, "Ns/m^3");

	mpParameters->addParameter("Use cell-cell friction", CSParameter::Bool, &mUseCCFriction, "");

	// Cell-vs-Extracellular Matrix gamma
	mpParameters->addParameter("Cell-ECM gamma", CSParameter::Double, &biolink->ecm_gamma_bio, "Ns/m^3");
	mpParameters->addParameter("Cell-Medium gamma", CSParameter::Double, &biolink->medium_gamma_bio, "Ns/m^3");

	// for ScenarioEmbeddingmedium:
	mpParameters->addParameter("Cell Population Radius (Embedding Medium)", CSParameter::Double, &mPopulationRadius, "Cell Diameters");

	mpParameters->addParameter("Cell Distance in Initial Population (Embedding Medium)", CSParameter::Double, &mPopulationInitialDistance, "Cell Diameters");

	CSParameterChoice * choiceLobuleShape = new CSParameterChoice(mLobuleShapeNames, 0);
	mpParameters->setParameter("Lobule Shape", CSParameter::Choice, choiceLobuleShape, "");
	mpParameters->addParameter("Lobule shape : radius", CSParameter::Double, &this->mlobule_radius, "");
	mpParameters->addParameter("Lobule shape : height", CSParameter::Double, &this->mlobule_height, "");

	mpParameters->addParameter("Blood Vessel Network", CSParameter::Bool, &this->mBloodVesselNetwork, "");
	mpParameters->addParameter("Path to Blood Vessel Network", CSParameter::FileName, &this->mBloodVesselNetworkPath, "");

	mpParameters->addParameter("Read Cell Position", CSParameter::Bool, &this->mReadCells, "");

	mpParameters->addParameter("Pressure Quiescent->Proliferating", CSParameter::Double, &biolink->quiescentMinPressure_bio, "Pa");
	mpParameters->addParameter("Pressure Proliferating->Quiescent", CSParameter::Double, &biolink->quiescentMaxPressure_bio, "Pa");
	mpParameters->setParameter("Overlap threshold Proliferating->Quiescent", CSParameter::Double, &overlap_threshold, "');");
	mpParameters->setParameter("Overlap threshold Quiescent->Proliferating", CSParameter::Double, &overlap_threshold_min, "');");
	mpParameters->setParameter("Density threshold Proliferating->Quiescent", CSParameter::Double, &density_threshold, "');");
	mpParameters->setParameter("Overlap duration Proliferating->Quiescent", CSParameter::Double, &overlap_duration_quiescence, "s");
	mpParameters->setParameter("Overlap duration Quiescent->Necrotic", CSParameter::Double, &overlap_duration_necrotic, "s");

	mpParameters->setParameter("stateMechanism", CSParameter::Int, &mStateChoice, "");
	mpParameters->setParameter("Force Threshold", CSParameter::Double, &mForceThreshold, "nN");
	mpParameters->setParameter("Density Threshold", CSParameter::Double, &mDensityThreshold, "");

	mpParameters->addParameter("dynamic timestep", CSParameter::Bool, &this->mUseDynamicalTimeSteps, "");
	mpParameters->addParameter("time step", CSParameter::Double, &this->timeStep, "s");
	mpParameters->addParameter("time step scaling factor", CSParameter::Int, &this->mTimeScalingFactor, "");
	mpParameters->addParameter("max Displacement", CSParameter::Double, &this->mMaximumDisplacementCellRadius, "");
	mpParameters->addParameter("dynamic timestep min", CSParameter::Double, &this->mDynamicTimeStepMin, "s");
	mpParameters->addParameter("dynamic timestep max", CSParameter::Double, &this->mDynamicTimeStepMax, "s");

	mpParameters->addParameter("output path", CSParameter::DirName, &this->mOutputPath, "");
	mpParameters->addParameter("output suffix", CSParameter::String, &this->mOutputPrefix, "");
	mpParameters->addParameter("observe division", CSParameter::Bool, &this->observeDivision, "");
}

void ModelCellsSpherical::InitParameters(CSParameterContext *parms)
{
	if (!parms)
	{
		if (!mpParameters) {
			RegisterParameters();
		}
		DefaultParameters(mpParameters);
		return;
	}

	std::vector<CSParameter *> parameters = parms->getParameters();

	std::vector<CSParameter *>::const_iterator parmsIt;

	for (parmsIt = parameters.begin(); parmsIt != parameters.end(); ++parmsIt)
	{
		if ((*parmsIt)->name() == "Contact Model")
		{
			CSParameter * foundParm = mpParameters->findParameter("Contact Model");
			((CSParameterChoice *)foundParm->dataPointer())->setCurrentIndex(((CSParameterChoice *)(*parmsIt)->dataPointer())->currentIndex());
		}
		else if ((*parmsIt)->name() == "Single Bond Energy")
			biolink->single_bond_energy_bio = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Cell Diameter")
			biolink->diameter_bioCells = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Length Scale")
			biolink->length_scale = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Cell Cycle Time")
			biolink->cycletime_bio = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Cycle Time SD")
			biolink->cycletime_stddev_bio = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Use Dumbbell Division")
			this->mUseDumbbell = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Dumbbell time")
			this->mDefaultDumbbellPeriod = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Young Modulus Cells")
			biolink->youngModulus_bioCells = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Poisson Ratio Cells")
			biolink->poisson_ratio_bioCells = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Young Modulus Sinusoids")
			biolink->youngModulus_bioSinusoids = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Poisson Ratio Sinusoids")
			biolink->poisson_ratio_bioSinusoids = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Diffusion Constant")
			biolink->diffusion_constant_cells_bio = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Cell-Cell Adhesion")
			biolink->adhesion_to_cells_bio = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Use Cell Polarity")
			this->mUsePolarCells = *(bool*)(*parmsIt)->dataPointer();
		else if ((*parmsIt)->name() == "Cell-Cell gamma for shear friction")
			biolink->cell_gamma_perpendicular_bio = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Cell-Cell gamma for normal friction")
			biolink->cell_gamma_parallel_bio = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Use cell-cell friction")
			this->mUseCCFriction = *(bool*)(*parmsIt)->dataPointer();
		else if ((*parmsIt)->name() == "Cell-Medium gamma")
			biolink->medium_gamma_bio = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Cell-ECM gamma")
			biolink->ecm_gamma_bio = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Cell Population Radius (Embedding Medium)")
			mPopulationRadius = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Cell Distance in Initial Population (Embedding Medium)")
			mPopulationInitialDistance = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Lobule shape : radius")
			this->mlobule_radius = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Lobule shape : height")
			this->mlobule_height = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Blood Vessel Network")
			this->mBloodVesselNetwork = *(bool *)(*parmsIt)->dataPointer();
		else if ((*parmsIt)->name() == "Path to Blood Vessel Network")
			this->mBloodVesselNetworkPath = (*parmsIt)->dataString();
		else if ((*parmsIt)->name() == "Read Cell Position")
			this->mReadCells = *(bool *)(*parmsIt)->dataPointer();
		else if ((*parmsIt)->name() == "Pressure Quiescent->Proliferating")
			biolink->quiescentMinPressure_bio = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Pressure Proliferating->Quiescent")
			biolink->quiescentMaxPressure_bio = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Overlap threshold Proliferating->Quiescent")
			overlap_threshold = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Density threshold Proliferating->Quiescent")
			density_threshold = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Overlap threshold Quiescent->Proliferating")
			overlap_threshold_min = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Overlap duration Proliferating->Quiescent")
			overlap_duration_quiescence = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Overlap duration Quiescent->Necrotic")
			overlap_duration_necrotic = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "time step")
			this->timeStep = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "dynamic timestep")
			this->mUseDynamicalTimeSteps = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "time step scaling factor")
			this->mTimeScalingFactor = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "max Displacement")
			this->mMaximumDisplacementCellRadius = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "dynamic timestep min")
			this->mDynamicTimeStepMin = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "dynamic timestep max")
			this->mDynamicTimeStepMax = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "output suffix")
			this->mOutputPrefix = (*parmsIt)->dataString();
		else if ((*parmsIt)->name() == "observe division")
			this->observeDivision = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "output path")
			this->mOutputPath = (*parmsIt)->dataString();
		else if ((*parmsIt)->name() == "stateMechanism")
			mStateChoice = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Force Threshold")
			mForceThreshold = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Density Threshold")
			mDensityThreshold = (*parmsIt)->value();
		else if ((*parmsIt)->name() == "Lobule Shape")
		{
			CSParameter * foundParm = mpParameters->findParameter("Lobule Shape");
			((CSParameterChoice *)foundParm->dataPointer())->setCurrentIndex(((CSParameterChoice *)(*parmsIt)->dataPointer())->currentIndex());
		}
		else {
			std::cerr << "ModelCellsSpherical::InitParameters:  Unknown parameter given:" << "\t"
				<< (*parmsIt)->name() << std::endl;
		}
	}
}

void
ModelCellsSpherical::DefaultParameters(CSParameterContext *emptyContext)
{
	if (!emptyContext)
		return;

	// what interaction model to use:
	CSParameterChoice * choiceContactModel = new CSParameterChoice(mContactModelNames, 1);
	emptyContext->setParameter("Contact Model", CSParameter::Choice, choiceContactModel, "");

	// length scale
	emptyContext->setParameter("Cell Diameter", CSParameter::Double, new double(1.4e-5), "m");

	// cycle time
	emptyContext->setParameter("Cell Cycle Time", CSParameter::Double, new double(79200.), "s");

	// std deviation of cycle times
	emptyContext->setParameter("Cycle Time SD", CSParameter::Double, new double(7200.), "s");

	// use an extending dumbbell for a period of time before the cell cleaves into two
	emptyContext->setParameter("Use Dumbbell Division", CSParameter::Bool, new bool(true), "");

	//dumbbell time period
	emptyContext->addParameter("Dumbbell time", CSParameter::Double, new double(1800.), "s");

	// Young Modulus
	emptyContext->setParameter("Young Modulus Cells", CSParameter::Double, new double(450.), "Pa");

	// Poisson Ratio
	emptyContext->setParameter("Poisson Ratio Cells", CSParameter::Double, new double(.4), "");

	// Young Modulus
	emptyContext->setParameter("Young Modulus Sinusoids", CSParameter::Double, new double(600.), "Pa");

	// Poisson Ratio
	emptyContext->setParameter("Poisson Ratio Sinusoids", CSParameter::Double, new double(.4), "");

	// diffusion constant for cells
	emptyContext->setParameter("Diffusion Constant", CSParameter::Double, new double(1.e-16), "m^2/s");

	// cell-cell adhesion
	emptyContext->setParameter("Cell-Cell Adhesion", CSParameter::Double, new double(1.e15), "m^-2");

	// use simple CellSphericals or CellSphericalPolar
	emptyContext->setParameter("Use Cell Polarity", CSParameter::Bool, new bool(false), "");

	// cell-cell gamma
	emptyContext->setParameter("Cell-Cell gamma for shear friction", CSParameter::Double, new double(1.e8), "Ns/m^3");
	emptyContext->setParameter("Cell-Cell gamma for normal friction", CSParameter::Double, new double(1.e8), "Ns/m^3");
	emptyContext->setParameter("Use cell-cell friction", CSParameter::Bool, new bool(true), "");
	// Cell-vs-Extracellular Matrix gamma

	emptyContext->setParameter("Cell-Medium gamma", CSParameter::Double, new double(1.e9), "Ns/m^3");

	CSParameterChoice * choiceLobuleShape = new CSParameterChoice(mLobuleShapeNames, 0);
	// for ScenarioEmbeddingmedium:
	emptyContext->setParameter("Cell Population Radius (Embedding Medium)", CSParameter::Double, new double(20), "Cell Diameters");

	emptyContext->setParameter("Cell Distance in Initial Population (Embedding Medium)", CSParameter::Double, new double(.95), "Cell Diameters");


	emptyContext->setParameter("Lobule Shape", CSParameter::Choice, choiceLobuleShape, "");

	emptyContext->setParameter("Lobule shape : radius", CSParameter::Double, new double(10), "");
	emptyContext->setParameter("Lobule shape : height", CSParameter::Double, new double(10), "");

	emptyContext->setParameter("Blood Vessel Network", CSParameter::Bool, new bool(false), "");
	emptyContext->setParameter("Path to Blood Vessel Network", CSParameter::FileName, new std::string("../../input/default2_md000011.mxf"), "");

	emptyContext->setParameter("Read Cell Position", CSParameter::Bool, new bool(false), "");

	emptyContext->setParameter("Pressure Quiescent->Proliferating", CSParameter::Double, new double(80), "Pa");
	emptyContext->setParameter("Pressure Proliferating->Quiescent", CSParameter::Double, new double(100), "Pa");
	emptyContext->setParameter("Overlap threshold Proliferating->Quiescent", CSParameter::Double, new double(.81), "');");
	emptyContext->setParameter("Overlap duration Proliferating->Quiescent", CSParameter::Double, new double(100), "s");
	emptyContext->setParameter("Overlap duration Quiescent->Necrotic", CSParameter::Double, new double(604800), "s");
	emptyContext->setParameter("Oxygen level for Proliferating->Quiescent", CSParameter::Double, new double(0.111), "");
	emptyContext->setParameter("Oxygen level for Quiescent->Necrotic", CSParameter::Double, new double(0.0333), "");
	emptyContext->setParameter("Oxygen uptake rate", CSParameter::Double, new double(.16), "");
	emptyContext->setParameter("Initial oxygen concentration", CSParameter::Double, new double(16.5), "");
	emptyContext->setParameter("time step", CSParameter::Double, new double(1.), "s");
	emptyContext->setParameter("dynamic timestep", CSParameter::Bool, new bool(true), "");
	emptyContext->setParameter("time step scaling factor", CSParameter::Int, new unsigned int(2), "");
	emptyContext->setParameter("max Displacement", CSParameter::Double, new double(0.1), "");
	emptyContext->setParameter("dynamic timestep min", CSParameter::Double, new double(1.), "s");
	emptyContext->setParameter("dynamic timestep max", CSParameter::Double, new double(256.), "s");
	emptyContext->setParameter("Single bond energy", CSParameter::Double, new double(1e-17), "J");

	emptyContext->addParameter("output path", CSParameter::DirName, new std::string("./output/"), "");
	emptyContext->addParameter("output suffix", CSParameter::String, new std::string("sphericalCells_sim1000"), "");
	emptyContext->addParameter("observe division", CSParameter::Bool, new bool(false), "");

}

CSParameterContext * ModelCellsSpherical::GetParameters(std::string contextName)
{
	if (!mpParameters)
	{
		RegisterParameters();
		return mpParameters;
	}

	if (!contextName.size())
		return mpParameters;

	if (contextName == name)
		return mpParameters;

	return mpParameters->findContext(contextName);
}

void ModelCellsSpherical::xmlParameters(QXmlStreamReader * xml,
	std::stringstream & warnings,
	std::stringstream & errors)
{
	Q_ASSERT(xml->name() == "parameters" && xml->isStartElement());
	CSParameterContextTemporary * parms = nullptr;

	// load defaults
	CSParameterContextTemporary defaultParms("");
	DefaultParameters(&defaultParms);
	InitParameters(&defaultParms);

	parms = (CSParameterContextTemporary *)CSParameterContext::createFromXML(xml);
	if (parms)
	{
		InitParameters((CSParameterContext *)parms);
		delete parms;
	}
}

void ModelCellsSpherical::xmlSimulation(QXmlStreamReader * xml,
	std::stringstream & warnings,
	std::stringstream & errors)
{
	Q_ASSERT(xml->name() == "simulation" && xml->isStartElement());
	CSParameterContextTemporary * parms = nullptr;

	parms = (CSParameterContextTemporary *)CSParameterContext::createFromXML(xml);
	if (parms)
	{
		CSParameter * untilDays = parms->findParameter("Simulation time");
		if (untilDays) {
			simulateUntilDays = untilDays->value();
		}
		CSParameter * observe = parms->findParameter("Observation interval");
		enableObservation = false;
		if (observe)
		{
			observeEveryDays = observe->value();
			enableObservation = observe->value() > 0;
		}

		CSParameter * scenario = parms->findParameter("Scenario");
		if (scenario)
		{
			mScenario = (ModelCellsSpherical::Scenario)(int)scenario->value();
		}
		CSParameter * ncells = parms->findParameter("Number of cells");
		if (ncells)
			nof_cells = ncells->value();
		CSParameter * dist = parms->findParameter("Two cells - distance");
		if (dist)
			scen_two_cells_dist = dist->value();
		CSParameter * radius = parms->findParameter("Initial radius");
		if (radius)
			initial_radius = radius->value();
		CSParameter * mindist = parms->findParameter("Minimum cell distance");
		if (mindist)
			minimum_cell_distance = mindist->value();
		delete parms;
	}
}

ModelCellsSpherical * ModelCellsSpherical::createFromXML(QXmlStreamReader * xml,
	std::stringstream & /* errors */,
	std::stringstream & warnings)
{
	Q_ASSERT(xml->name() == "Model" && xml->isStartElement());
	ModelCellsSpherical * model = new ModelCellsSpherical();

	// to be on the safe side query 'dimensions' again (even though done in Model::createFromXML())
	unsigned int dimensions = xml->attributes().value("dimensions").toString().toUInt();
	model->is2D = (dimensions == 2) ? true : false;
	return model;
}

// a shortcut for writing Parameter entries directly w/o a CSParameterContext
#define WRITE_XML_DOUBLE_PARM( name, value, unit )              \
{                                                           \
    xml->writeStartElement( "Parameter" );                      \
    xml->writeAttribute( "name", name );                        \
    xml->writeAttribute( "type", "Double" );                    \
    xml->writeAttribute( "value", QString("%1").arg(value) );   \
    xml->writeAttribute( "unit", unit );                        \
    xml->writeEndElement(); }{}

void ModelCellsSpherical::writeXML(QXmlStreamWriter *xml) const
{
#undef QT_NO_CAST_FROM_ASCII

	xml->writeStartElement("Model");
	xml->writeAttribute("type", xmlType.c_str());
	xml->writeAttribute("name", name.c_str());
	xml->writeAttribute("dimensions", is2D ? "2" : "3");

	xml->writeStartElement("simulation");
	WRITE_XML_DOUBLE_PARM("Simulation time", simulateUntilDays, "d");
	if (enableObservation)
		WRITE_XML_DOUBLE_PARM("Observation interval", observeEveryDays, "d");
	xml->writeEndElement();

	// parameters:
	xml->writeStartElement("parameters");

	if (mpParameters)
		mpParameters->writeXML(xml, true);

	xml->writeEndElement();

	// cells:  left here for debugging purposes only - OBSOLETE
	xml->writeStartElement("cells");

	std::vector<CellSpherical *>::const_iterator cell = cells.begin();

	while (cell != cells.end())
	{
		(*cell)->writeXML(xml);
		++cell;
	}

	xml->writeEndElement(); // cells

	xml->writeEndElement();
}

void ModelCellsSpherical::readModelData(H5::H5File * inputFile,
	std::stringstream & errors,
	std::stringstream & warnings)
{
	hsize_t numElements;  // used for reading data extents and allocating space
	cells.clear();
	std::string groupString = "/" + name;

	try
	{
		// read in the state of the Marsaglia random number generator:
		std::string dataSetRNGString = groupString + "/MarsagliaRNG";

		H5::DataSet dataSetRNG =
			inputFile->openDataSet(dataSetRNGString);

		H5::CompType rngDataType;
		Random::HDF5DataFormat(rngDataType);

		dataSetRNG.read(&mRandom, rngDataType);
	}
	catch (...)
	{
		std::cout << "No Random Number Generator data in Data file.  Using default seed\n";
		mRandom.Init();
	}



	// reading in barriers:
	std::string dataSetBarriersString =
		groupString + "/ModelElementBarrierTriangles";

	try
	{
		H5::DataSet dataSetBarriers =
			inputFile->openDataSet(dataSetBarriersString);

		H5::CompType dataTypeBarriers =
			dataSetBarriers.getCompType();

		H5::CompType readTypeBarrier =
			ModelElementBarrierTriangle::ParseHDF5DataFormat(dataTypeBarriers,
				errors,
				warnings);

		H5::DataSpace outputSpaceBarriers = dataSetBarriers.getSpace();

		outputSpaceBarriers.getSimpleExtentDims(&numElements, nullptr);

		ModelElementBarrierTriangle * readBarriers =
			new ModelElementBarrierTriangle[numElements];

		dataSetBarriers.read(readBarriers, readTypeBarrier);

		for (unsigned int i = 0; i<numElements; ++i)
		{
			ModelElementBarrierTriangle * barrier = readBarriers + i;
			barrier->setNormalVector();
			barrier->setBoundingBox(.1);
			mBarriers.push_back(barrier);
			cells2->add(barrier);
			mpArena->addObject(barrier->GLObject());
		}
	}
	catch (...)
	{
		std::cout << "No data set found named \"" << dataSetBarriersString << "\".\n";
	}

	// reading in vessel graph / vessel spheres:
	std::string dataSetVesselSpheresString =
		groupString + "/ModelElementVesselSpheres";

	try
	{
		H5::DataSet dataSetVesselSpheres =
			inputFile->openDataSet(dataSetVesselSpheresString);

		H5::CompType dataTypeVesselSphere =
			dataSetVesselSpheres.getCompType();

		H5::CompType readTypeVesselSphere =
			ModelElementVesselSphere::ParseHDF5DataFormat(dataTypeVesselSphere,
				errors,
				warnings);

		H5::DataSpace outputSpaceVesselSpheres = dataSetVesselSpheres.getSpace();

		outputSpaceVesselSpheres.getSimpleExtentDims(&numElements, nullptr);

		ModelElementVesselSphere * readVesselSpheres =
			new ModelElementVesselSphere[numElements];

		dataSetVesselSpheres.read(readVesselSpheres, readTypeVesselSphere);

		if (mpGraphBloodVesselNetwork)
			delete mpGraphBloodVesselNetwork;

		GraphSphere * network = new GraphSphere();

		for (unsigned int i = 0; i<numElements; ++i)
		{
			ModelElementVesselSphere * vesselSphere = readVesselSpheres + i;
			network->mvNode.push_back(vesselSphere);
		}

		AddBloodVesselNetwork(network);
	}
	catch (...)
	{
		std::cout << "No data set found named \"" << dataSetVesselSpheresString << "\".\n";
	}


	try
	{
		// read in cells' data:
		std::string dataSetCellsSphericalString =
			groupString + "/CellsSpherical";

		H5::DataSet dataSetCellsSpherical =
			inputFile->openDataSet(dataSetCellsSphericalString);

		// ToDo:  set defaults for non-essential data fields

		// build up the data type in a flexible way

		// The data layout of the saved data
		H5::CompType dataTypeCellSpherical =
			dataSetCellsSpherical.getCompType();

		// the 'read type', i.e. the data layout for our read-in buffer:
		H5::CompType readType
			= CellSpherical::ParseHDF5DataFormat(dataTypeCellSpherical,
				errors,
				warnings);

		// how many data elements of type dataTypeCellSpherical are in the dataSet
		H5::DataSpace outputSpace = dataSetCellsSpherical.getSpace();
		outputSpace.getSimpleExtentDims(&numElements, nullptr);

		// allocate the read-out buffer
		CellSpherical * readCells = new CellSpherical[numElements];


		// the actual reading of the data into the readout array
		dataSetCellsSpherical.read(readCells, readType);

		// add the cells to this Model:
		for (unsigned int i = 0; i<numElements; ++i)
		{
			CellSpherical *cell = readCells + i;
			AddCell(cell);

			// The following is what UpdateParametersForAllCells would do
			// minus SetCycletimeGaussClamped which would use the Random number
			// generator (and change its internal state) and overwrite the cycle
			// time which has been read in.
			cell->divisionRadius = defaultDivisionCellRadius;
			cell->initialRadius = defaultInitialCellRadius;
			cell->UpdateDeltaRadius(timeStep);
			cell->youngModulus = defaultYoungModulusCells;
			cell->poissonRatio = defaultPoissonRatioCells;
			cell->mQuiescentMin = mQuiescentMin;
			cell->mQuiescentMax = mQuiescentMax;


			// why?  The cellcycleTime will be read in!
			// cell->SetCycleTimeGaussClamped(defaultCellCycleTime, defaultCellCycleTimeStandardDeviation);//added by eugenio
			//cell->setState(Cell::StateQuiescent, false);
			//cell->setState(Cell::StateDividing, true);
			cell->maxoverlap = 100;

			if ((long)cell->mDivisionPartnerIndex != -1)
			{
				cell->mpDivisionPartner = readCells + cell->mDivisionPartnerIndex;
				cell->mDaughterCell = (cell->mDivisionPartnerIndex < i);
				cell->setState(Cell::StateDividing);
			}
			else
				cell->mpDivisionPartner = nullptr;
		}
	}
	catch (...)
	{
	}

	try
	{
		// read in cells' data:
		std::string dataSetCellsSphericalString =
			groupString + "/CellsSphericalPolar";

		H5::DataSet dataSetCellsSpherical =
			inputFile->openDataSet(dataSetCellsSphericalString);

		// ToDo:  set defaults for non-essential data fields

		// build up the data type in a flexible way

		// The data layout of the saved data
		H5::CompType dataTypeCellSpherical =
			dataSetCellsSpherical.getCompType();

		// the 'read type', i.e. the data layout for our read-in buffer:
		H5::CompType readType
			= CellSphericalPolar::ParseHDF5DataFormat(dataTypeCellSpherical,
				errors,
				warnings);

		// how many data elements of type dataTypeCellSpherical are in the dataSet
		H5::DataSpace outputSpace = dataSetCellsSpherical.getSpace();
		outputSpace.getSimpleExtentDims(&numElements, nullptr);

		CellSphericalPolar * readCells = new CellSphericalPolar[numElements];

		// the actual reading of the data into the readout array
		dataSetCellsSpherical.read(readCells, readType);

		// add the cells to this Model:
		for (unsigned int i = 0; i<numElements; ++i)
		{
			CellSphericalPolar *cell = readCells + i;
			AddCell(cell);

			if ((long)cell->mDivisionPartnerIndex != -1)
			{
				cell->mpDivisionPartner = readCells + cell->mDivisionPartnerIndex;
				cell->mDaughterCell = (cell->mDivisionPartnerIndex < i);
				cell->setState(Cell::StateDividing);
			}
			else
				cell->mpDivisionPartner = nullptr;
		}

	}
	catch (...)
	{
	}

	UpdateParametersForAllCells();
}

void ModelCellsSpherical::writeHDF5(H5::H5File * outputFile) const
{
	// Todo: error:
	if (!outputFile)
		return;

	// dimensions - have to be handed over as pointers:
	hsize_t dims[] = { 1 };
	hsize_t cellsSphericalMaxDims[] = { cells.size() };
	// the memory space of a single data element
	H5::DataSpace memspace(1, dims, nullptr);

	std::string groupString = "/" + name;
	H5::Group modelGroup(outputFile->createGroup(groupString));

	// preparing the Marsaglia random generator to be saved:
	H5::CompType rngDataType(sizeof(Random));
	Random::HDF5DataFormat(rngDataType);

	H5::DataSpace rngDataSpace(1, dims, nullptr);

	H5::DataSet rngDataSet(outputFile->createDataSet(groupString + "/MarsagliaRNG",
		rngDataType,
		rngDataSpace));

	rngDataSet.write(&mRandom, rngDataType, memspace, rngDataSpace);


	// creating the compound type and defining its composition
	// (by calling CellSpherical::HDF5Dataformat())

	H5::DSetCreatPropList cparms;

	hsize_t offset = 0;
	hsize_t count = 1;

	// hack for now until std::map< ModelElement::Type,
	// std::vector<ModelElement*> > is finished and merged:
	if (cells[0]->mType == ModelElement::TypeCellSpherical)
	{
		H5::CompType cellSphericalType(sizeof(CellSpherical));
		CellSpherical::HDF5DataFormat(cellSphericalType);

		H5::DataSpace cellSphericalDataSpace(1, cellsSphericalMaxDims, nullptr);

		H5::DSetCreatPropList cparms;
		cparms.setChunk(1, dims);
		H5::DataSet cellSphericalDataSet(outputFile->createDataSet(groupString + "/CellsSpherical",
			cellSphericalType,
			cellSphericalDataSpace,
			cparms));

		hsize_t offset = 0;
		hsize_t count = 1;

		std::vector<CellSpherical *>::const_iterator cIt;

		for (cIt = cells.begin(); cIt != cells.end(); ++cIt)
		{
			// preliminary:  determine dumbbell divisionPartnerIndex
			if ((*cIt)->mpDivisionPartner)
			{
				std::vector<CellSpherical *>::const_iterator found = find(cells.begin(), cells.end(), (*cIt)->mpDivisionPartner);
				if (found != cells.end())
					(*cIt)->mDivisionPartnerIndex = found - cells.begin();
				else
					(*cIt)->mDivisionPartnerIndex = -1;
			}
			else
				(*cIt)->mDivisionPartnerIndex = -1;

			cellSphericalDataSpace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
			cellSphericalDataSet.write(*cIt, cellSphericalType, memspace, cellSphericalDataSpace);
			outputFile->flush(H5F_SCOPE_LOCAL);
			offset++;
		}
	}
	else // hopefully only CellSphericalPolar were allocated :)
	{
		H5::CompType cellSphericalType(sizeof(CellSphericalPolar));
		CellSphericalPolar::HDF5DataFormat(cellSphericalType);

		H5::DataSpace cellSphericalDataSpace(1, cellsSphericalMaxDims, nullptr);

		cparms.setChunk(1, dims);
		H5::DataSet cellSphericalDataSet(outputFile->createDataSet(groupString + "/CellsSphericalPolar",
			cellSphericalType,
			cellSphericalDataSpace,
			cparms));

		std::vector<CellSpherical *>::const_iterator cIt;

		for (cIt = cells.begin(); cIt != cells.end(); ++cIt)
		{
			// preliminary:  determine dumbbell divisionPartnerIndex
			if ((*cIt)->mpDivisionPartner)
			{
				std::vector<CellSpherical *>::const_iterator found = find(cells.begin(), cells.end(), (*cIt)->mpDivisionPartner);
				if (found != cells.end())
					(*cIt)->mDivisionPartnerIndex = found - cells.begin();
				else
					(*cIt)->mDivisionPartnerIndex = -1;
			}
			else
				(*cIt)->mDivisionPartnerIndex = -1;

			cellSphericalDataSpace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
			cellSphericalDataSet.write(*cIt, cellSphericalType, memspace, cellSphericalDataSpace);
			outputFile->flush(H5F_SCOPE_LOCAL);
			offset++;
		}
	}


	H5::CompType vesselSphereType(sizeof(ModelElementVesselSphere));
	ModelElementVesselSphere::HDF5DataFormat(vesselSphereType);

	hsize_t vesselSpheresMaxDims[] = { (mpGraphBloodVesselNetwork)
		? mpGraphBloodVesselNetwork->mvNode.size()
		: 0
	};
	H5::DataSpace vesselSpheresDataSpace(1, vesselSpheresMaxDims, nullptr);

	H5::DataSet vesselSpheresDataSet(outputFile->createDataSet(groupString + "/ModelElementVesselSpheres",
		vesselSphereType,
		vesselSpheresDataSpace,
		cparms));

	if (mpGraphBloodVesselNetwork)
	{
		offset = 0;
		std::vector<ModelElementVesselSphere *>::const_iterator vesselIt;
		for (vesselIt = mpGraphBloodVesselNetwork->mvNode.begin();
			vesselIt != mpGraphBloodVesselNetwork->mvNode.end();
			++vesselIt)
		{
			vesselSpheresDataSpace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
			vesselSpheresDataSet.write(*vesselIt,
				vesselSphereType,
				memspace,
				vesselSpheresDataSpace);
			outputFile->flush(H5F_SCOPE_LOCAL);
			++offset;
		}
	}

	H5::CompType barrierType(sizeof(ModelElementBarrierTriangle));
	ModelElementBarrierTriangle::HDF5DataFormat(barrierType);

	hsize_t barrierMaxDims[] = { mBarriers.size() };
	H5::DataSpace barrierDataSpace(1, barrierMaxDims, nullptr);

	H5::DataSet barrierDataSet(outputFile->createDataSet(groupString + "/ModelElementBarrierTriangles",
		barrierType,
		barrierDataSpace,
		cparms));

	offset = 0;
	ModelElementBarrierTriangle ** barrIt;
	for (barrIt = mBarriers.begin(); barrIt != mBarriers.end(); ++barrIt)
	{
		barrierDataSpace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
		barrierDataSet.write(*barrIt, barrierType, memspace, barrierDataSpace);
		outputFile->flush(H5F_SCOPE_LOCAL);
		++offset;
	}

}

#pragma region Visualization

void ModelCellsSpherical::UpdateCellsStaining(int mode)
{
#pragma region 0: All white
	if (mode == 0)
	{
		for (unsigned int i = 0; i<cells.size(); i++)
		{
			if (cells.at(i)->mSubType == SubType::ETERNAL_QUIESCENT) 
				cells.at(i)->SetColor(227. / 255., 179. / 255., 94. / 255.);
			else
			{
				cells.at(i)->SetColor(1, 1, 1);
			}
		}
	}
#pragma endregion

#pragma region 1: Cell cycle state

	else if (mode == 1)
	{
		for (unsigned int i = 0; i<cells.size(); i++)
		{
			if (cells.at(i)->mSubType == SubType::ETERNAL_QUIESCENT) cells.at(i)->SetColor(227. / 255., 179. / 255., 94. / 255.);
			else
			{
				if (cells.at(i)->getState(Cell::StateQuiescent))
				{
					cells.at(i)->SetColor(0.5, 0.5, 0.5);
				}
				else if (cells.at(i)->getState(Cell::StateNecrotic))
				{
					cells.at(i)->SetColor(0, 0, 0);
				}
				else if (cells.at(i)->getState(Cell::StateDividing))
				{
					cells.at(i)->SetColor(0, 1, 0);
				}
				else
				{
					cells.at(i)->SetColor(1, 1, 1);
				}
			}
		}
	}
#pragma endregion

#pragma region 2: Last absolute force

	else if (mode == 2)
	{
		for (unsigned int i = 0; i<cells.size(); i++)
		{
			if (cells.at(i)->mSubType == SubType::ETERNAL_QUIESCENT) cells.at(i)->SetColor(227. / 255., 179. / 255., 94. / 255.);
			else
			{
				(*core->tools).color->CreateSuperTrafficHue(cells.at(i)->lastPressure, 0., 20000.);
				cells.at(i)->SetColor((*core->tools).color->r, (*core->tools).color->g, (*core->tools).color->b);
			}
		}
	}

#pragma endregion

#pragma region 3: Cell volume

	else if (mode == 3)
	{
		for (unsigned int i = 0; i<cells.size(); i++)
		{
			if (cells.at(i)->mSubType == SubType::ETERNAL_QUIESCENT) cells.at(i)->SetColor(227. / 255., 179. / 255., 94. / 255.);
			else
			{
				// Prelim: Later: Actual volume calculation
				(*core->tools).color->CreateTrafficHue(cells.at(i)->mRadius, defaultInitialCellRadius, defaultDivisionCellRadius);
				cells.at(i)->SetColor((*core->tools).color->r, (*core->tools).color->g, (*core->tools).color->b);
			}
		}
	}


#pragma endregion

	else if (mode == 4) // relative position within a lobule
	{
		for (unsigned int i = 0; i<cells.size(); ++i)
		{
			(*core->tools).color->CreateTrafficHue(cells[i]->mLayerIndicator, 0, 1);
			cells[i]->SetColor((*core->tools).color->r, (*core->tools).color->g, (*core->tools).color->b);
		}
	}
	else if (mode == 5)
	{
		for (unsigned int i = 0; i<cells.size(); ++i)
		{
			if (cells[i]->mLesionEdge)
				cells[i]->SetColor(1, 0, 0);
			else
				cells[i]->SetColor(1, 1, 1);
		}
	}
	// Remember colormode
	lastColormodeCells = mode;
}

void ModelCellsSpherical::UpdateCellsStaining()
{
	UpdateCellsStaining(lastColormodeCells);
}

void ModelCellsSpherical::UpdateCellsStaining(CellSpherical *cell)
{
#pragma region 0: All white

	switch (lastColormodeCells)
	{
	case 0:
		if (cell->mSubType == SubType::ETERNAL_QUIESCENT) cell->SetColor(227. / 255., 179. / 255., 94. / 255.);
		else
		{
			cell->SetColor(1, 1, 1);
		}
		break;
#pragma endregion

#pragma region 1: Cell cycle state

	case 1:
		if (cell->mSubType == SubType::ETERNAL_QUIESCENT) cell->SetColor(227. / 255., 179. / 255., 94. / 255.);
		else
		{
			if (cell->getState(Cell::StateQuiescent))
			{
				cell->SetColor(0.5, 0.5, 0.5);
			}
			else
			{
				cell->SetColor(1, 1, 1);
			}
		}
		break;

#pragma endregion

#pragma region 2: Last absolute force

	case 2:
		if (cell->mSubType == SubType::ETERNAL_QUIESCENT) cell->SetColor(227. / 255., 179. / 255., 94. / 255.);
		else
		{
			(*core->tools).color->CreateSuperTrafficHue(cell->lastPressure, 0., 20000.);
			cell->SetColor((*core->tools).color->r, (*core->tools).color->g, (*core->tools).color->b);
		}
		break;

#pragma endregion

#pragma region 3: Cell volume

	case 3:
		if (cell->mSubType == SubType::ETERNAL_QUIESCENT) cell->SetColor(227. / 255., 179. / 255., 94. / 255.);
		else
		{
			// Prelim: Later: Actual volume calculation
			(*core->tools).color->CreateTrafficHue(cell->mRadius, defaultInitialCellRadius, defaultDivisionCellRadius);
			cell->SetColor((*core->tools).color->r, (*core->tools).color->g, (*core->tools).color->b);
		}
		break;
	}
	cell->SetAlpha(0.5);
#pragma endregion
}

void ModelCellsSpherical::setVisible(bool visible)
{

	for (unsigned int i = 0; i < this->cells.size(); i++)
	{
		if (this->cells[i]->mVisible == true)
		{
			mpArena->removeObject(cells[i]->GLObject());
		}
		else
		{
			mpArena->addObject(this->cells[i]->GLObject());
		}

		this->cells[i]->mVisible = visible;
	}

	mpArena->draw();

}

void ModelCellsSpherical::writeXML2()
{
}

void ModelCellsSpherical::spherePacking(double radiusSphere)
{
	int counter = 0;
	double radius = 0.5;
	int max = (radiusSphere + 3) * 2;
	double max_h = radiusSphere;

	int I = max;
	int J = max;
	int K = max;

	Vector3f pos;

	for (int i = 0; i < I; i++)
	{
		for (int j = 0; j < J; j++)
		{
			for (int k = 0; k < K; k++)
			{
				pos.x = (2.*i + ((j + k) % 2))         * radius - max_h;
				pos.y = sqrt(3.)*(j + 1. / 3.*(k % 2)) * radius - max_h;
				pos.z = 2.*sqrt(6.) / 3.*k         * radius - max_h;

				if (pos.Norm() < radiusSphere) {
					AddCell(pos.x, pos.y, pos.z);
					std::cerr << "add cell : " << this->cells.size() << std::endl;
					counter++;
				}
				else
					std::cerr << "does not add cell (cell outside)" << std::endl;
			}
		}
	}
}
