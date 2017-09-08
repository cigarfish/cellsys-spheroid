#ifndef MODEL_CELLS_SPHERICAL_H
#define MODEL_CELLS_SPHERICAL_H

#include "../../model/Model/CSModel.h"
#include "../../tools/Tools.h"
#include "../../gui/CSGLArena.h"
#include "../../gui/GLTools/CSGLBar.h"
#include "../Cell/CellSpherical.h"
#include "../Cell/CellSphericalPolar.h"
#include "../Elements/ModelElementBarrierTriangle.h"
#include "../BiologyLink.h"
#include "../../../tools/random/Random.h"
#include "../../../tools/dataIO/vtp/CSVTPWriter.h"
#include "../../../tools/math/LinearSolver.h"
#include "../../Observation/Observation.h"

#include "../../BasicDatatypes/GraphSphere.h"

#include "../../../Core.h"
#include "../../../tools/math/CGP.h"

#include <vector>
#include <string>
#include <sstream>
#include <functional>

class CSParameterContext;
class QXmlStreamWriter;
class QXmlStreamReader;
class CSInteractionHertz;
class CSInteractionJKR;
class CSInteractionFrictionMatrix;
class ModelElementBarrierTriangle;
class BoundingBoxList;
class CSVTPWriter;

//! Model class for 2D monolayer models
class ModelCellsSpherical : public CSModel
{
	// preliminary:  declare Cell classes as friends to access parameters
	// to be replaced by cell-local parameters.
	friend class CellSpherical;
	friend class CellSphericalPolar;
	friend class monolayer;

#pragma region Parameters

public:
	//! The model type.
	static const std::string xmlType;

	enum ContactModel { ContactModelHertz, ContactModelJKR };
	enum LobuleShape { None, Hexagonal, Quadric };

	ContactModel mContactModel;
	LobuleShape mLobuleShape;

	static const std::string mContactModelNames[];
	static const std::string mLobuleShapeNames[];

	//! The individual model instantiation name (in case of multiple instances):
	std::string xmlName;

	//! Change the name
	// ...and the parameter context's name!
	void SetName(std::string newName)
	{
		if (newName != name)
		{
			core->models[newName] = this;
			if (name.size())
				core->models.erase(core->models.find(name));
			xmlName = newName;
			CSModel::SetName(newName);
		}
	};

	enum Scenario {
		ScenarioReadFromData = 0,
		ScenarioSingleCell,
		ScenarioGlucoseDiffusion,
		ScenarioEmbeddingMedium,
		ScenarioSpherePacking,
		ScenarioSphere,
		ScenarioTwoCells,
		ScenarioTestDiffusion,
		ScenarioTestManyCells
	};

	// Initialisation functions for the different scenarios
	void resetReadFromData();
	void resetSingleCell();
	void resetEmbeddingMedium();
	void resetSpherePacking();
	void resetSphere();
	void resetTwoCells();
	void resetTestManyCells();

	bool mManyCellsRandom = false;

	Scenario mScenario;
	void SetScenario(Scenario scen) { mScenario = scen; };

	// print internal model parameters debug
	// XXX MOVE THIS TO A TYPE OF REPORTER?
	void printParameters() const;

	// BioLink -> (default) parameters in Model class
	void UpdateParametersFromBioLink();

	// BioLink -> Dimensionless (default) parameters in Model class
	void UpdateDimensionlessParametersFromBioLink();

	// Push all dimensionless (default) parameters in Model class -> cell
	void UpdateParametersForCell(CellSpherical *);

	// Push all dimensionless (default) parameters in Model class -> All existing cells
	void UpdateParametersForAllCells();

	bool observeDivision;
	bool enableObservation;

protected:
	//border for lobule shape
	double mlobule_radius;
	double mlobule_height;

	bool mBloodVesselNetwork;
	std::string mBloodVesselNetworkPath;
	GraphSphere * mpGraphBloodVesselNetwork;

	bool mReadCells;

	// THIS SHOULD BE LIKE A KILL PLUGIN OR STH LIKE THAT!!
	//! Container used when a set of cells are to be killed.
	std::vector<CellSpherical*> mCellsToDelete;

	double mKillRate;
	//! Time interval for when to remove cells
	double mKillInterval;
	//! Next absolute time for when cells are removed, incremented by mKillInterval
	double mNextKillTime;
	//! The time period over which chosen cells will die. This period will be
	//! divided into mKillInterval large intervals.
	double mKillPeriod;
	//! Absolute time when the killing period will end.  After this time, all
	//! remaining cells will be removed.
	double mKillEnd;

	CSListContainer<ModelElementBarrierTriangle *> mBarriers;

	double nextProliferationUpdate;

	// XXX MOVE THESE SOMEWHERE ELSE
	// Parameters for overlap-based contact inhibited growth and necrosis
	double overlap_threshold = 0.75;
	double overlap_threshold_min = 0.75;
	double overlap_duration_quiescence = 0;
	double overlap_duration_necrotic = 604800;
	double density_threshold = 1.1;

	// Parameters for oxygen-based quiescence en necrosis;
	double oxygen_threshold_quiescence = 0.111;
	double oxygen_threshold_necrotic = 0.0333;

	// Parameters for sphere/disk initialization
	double initial_radius, minimum_cell_distance;
	int nof_cells;

	// Solver for the equation of motion
	LinearSolver * mSolver;

public:
	//! Default radius of the cells at t=0 and right after cell division
	double defaultInitialCellRadius = 0.55;

	//! Mean cell cycle time
	double defaultCellCycleTime;

	//! Standard deviation of cell cycle time (Gaussian distribution)
	double defaultCellCycleTimeStandardDeviation;

	//! Default radius beyond which the cells divide into two daughter cells
	double defaultDivisionCellRadius = 0.693;

	//! Default distance of each of the two daughter cells center to the center of the mother cell
	double defaultDivisionDistance = 0.45;

	//! Grows all cells and applies cell divisions (assuming one timestep)
	virtual void GrowAndDivide();

	virtual void AddPolarCell(double x, double y, double z);
	virtual void AddCell(double x, double y, double z);

	//! Helper method to set several defaults for AddCell(.,.,.) and
	//! AddPolarCell(.,.,.).
	void InitCell(CellSpherical *);

	//! Adds an existing cell
	virtual void AddCell(CellSpherical * newCell);

	//! Removes a cell from all containers.
	void RemoveCell(CellSpherical *);

	//Add Cells from MXF file
	void AddReadCells(std::string & filename);

	//Add shape of Lobule
	void AddBarrierTriangle(double mPoints[][3], bool visible = 1, double epsilon = 0.1, double r = 1, double g = 1, double b = 0, double a = 0.5);

	void AddBloodVesselNetwork(GraphSphere *);
	void AddBloodVesselNetwork(std::string & filename, int filetype = 2);

	void AddShape(double radius = 10, double height = 10, double r = 1, double g = 1, double b = 0, double a = 0.5, double epsilon = 0.1, bool visible = 1);
	void AddLobuleShape(double radius = 10, double height = 10, double r = 1, double g = 1, double b = 0, double a = 0.5, double epsilon = 0.1, bool visible = 1);
	void AddQuadricShape(double radius = 10, double height = 10, double r = 1, double g = 1, double b = 0, double a = 0.5, double epsilon = 0.1, bool visible = 1);


	double defaultYoungModulusSinusoids;      //! Default young modulus of sinusoids (dimensionless)
	double defaultPoissonRatioSinusoids;      //! Default poisson radtio of sinusoids

	double defaultYoungModulusCells;          //! Default young modulus of cells (dimensionless)
	double defaultPoissonRatioCells;          //! Default poisson radtio of cells
	double defaultDiffusionConstantCells;     //! Default diffusion constant of cells

											  // I'D LOVE TO MOVE THESE TO A XML DEFINITION FOR THE INTERACTIONS
											  //! Single bond energy (dimensionless)
	double singleBondEnergy;

	//! Default cell-cell adhesion density
	double adhesionDensity;

	//! Friction coefficient with other cells
	double gammaCells;

	double gammaCellsParallel;
	double gammaCellsPerpendicular;

	//! Fricition coefficient with ECM
	double gammaECM;
	double gammaMedium;

	// THESE SHOULD BE MOVED TO THE CELLS XXX
	//pressure for quiescent
	double mQuiescentMin;
	double mQuiescentMax;

	// XXX GET RID OF THESE
	// parameters for ScenarioEmbeddingmedium
	// Radius of edge of initial population in ScenarioEmbeddingmedium
	double mPopulationRadius;
	// distance between cells in initial configuration in Scenarioembeddingmedium
	double mPopulationInitialDistance;

	//! Resets directional forces
	void InitForces();
	void AddDirectedMotion(CellSpherical*);

	// Possible CSInteractions used in this model
	//  HertzForce interaction -> CSInteractionHertz
	//  ToDo:  JKR interaction
	//  FrictionMatrix from Hertz or JKR interactions -> CSInteractionFrictionMatrixHertz
	//  ...
	// These CSInteraction objects are not allocated in the constructor,
	//  but should be initialized in SetupSimulation depending on the options
	//  set by the user in the GUI or data bundle (xml).
	CSInteractionHertz          * mpInteractionHertz;
	CSInteractionJKR            * mpInteractionJKR;
	CSInteractionFrictionMatrix * mpInteractionFrictionMatrix;

	//! Loop over all cell-cell interactions.
	//! Calling CSInteractionHertz(.,.) and
	//! CSInteractionFrictionMatrixHertz(.,.)
	//! on every interacting Cell-Cell pair.
	virtual void UpdateInteractions();

	// XXX get rid of these?
	int mStateChoice;
	double mDensityThreshold;
	double mForceThreshold;

	void UpdateVelocitiesNeglectingCCFriction();

	//! Method to apply the Metropolis algorithm to the rotation of cells in
	//! dumbbell division phase.  The relevant energy is computed from the
	//! interactions.  Their values are computed by the interaction class
	//! CSInteractionHertzEnergy in GetInteractionEnergy().
	void DumbbellMetropolis();

	//! Calculate the interaction energy (Hertz or JKR, corresponding on the
	//! chosen interaction) by summing up the energies of all interactions of
	//! the cell given as argument.
	double GetInteractionEnergy(CellSpherical *);

	//! Loop over all cells.
	//! Calling UpdateCellFriction(.) and
	//! UpdateForcesLangevin(.) on every cell in
	//! std::vector<CellSpherical*> cells
	virtual void UpdateSingleCellEffects();

	//! Prepares cell friction based on accumulated contact areas
	virtual void UpdateCellFriction(CellSpherical *);

	//! Adds random Langevin-based forces
	virtual void UpdateForcesLangevin(CellSpherical *);

	//! Method to rescale the contribution of the Langevin force when adapting
	//! time step size.  Since the Langevin force scales with sqrt(1/timeStep).
	//! \param timeStepFactor  The factor by which the timeStep changed.
	virtual void rescaleLangevinContrib(double timeStepFactor);

public:
	Observation * observe;

	double scen_two_cells_dist;

	// XXX is this still needed?
	//! Flag for deciding to use CellSphericalPolar instead of CellSpherical.
	bool mUsePolarCells;

	//! Flag for deciding to use dumbbell or instantaneous division.
	bool mUseDumbbell;
	double mDefaultDumbbellPeriod;

	//  Flag for deciding to use cell-cell friction or not (Eugenio)
	bool mUseCCFriction;

	//! Flag for choosing dynamical step size adaptation per maximum allowed
	//! displacement for each time step.
	bool mUseDynamicalTimeSteps;
	double mDynamicTimeStepMin;//min and max Time step for dynamic timesteps
	double mDynamicTimeStepMax;

	//! Scaling factor for dynamical timeSteps.
	unsigned int mTimeScalingFactor;

	//! A variable to determine the maximum velocity in solveSystem() for time
	//! step adaptation.
	double mVelocityMaxSquared;

	//! A variable to hold the square of the maximum allowed displacement:
	//  (0.1 * defaultInitialCellRadius)^2.  Constant during simulation.
	double mMaximumDisplacementSquared;
	double mMaximumDisplacementCellRadius; //max allowed displacement

protected:

#pragma endregion

	// One simulation step. t -> t+1
	virtual void Simulate();

	void updatePosition();
	void setFrictionMatrix();
	void simulateCommon();
	void simulateWithNutrients();
	void SimulateGeneric();
	void SimulateTwoCells();
	void SimulateDiffusionTest();

	//! Array for keeping track which cells have been handled in UpdateInteractions.
	//  This is only necessary for the JKR contact model, in order not to account for contacts twice.
	bool * mpJKRElementDone;
	unsigned long int mJKRElementDoneSize;

	//! Vector holding the velocities (for solving the equation of motion)
	//  size = 3*cells.size()
	// TODO make this const!
	double * mpVelocities;
	unsigned long int mProblemAllocationSize;

	// vectors for solving the overdamped equation of motion with friction
	// i.e. if mpInteractionFrictionMatrix != NULL;
	// declared as members, so that reallocation has to happen rarely

	CSParameterContext * mpParameters;
	CSParameterContext * mpParametersVisualization;

	virtual void RegisterParameters();
	void UpdateDishContactForce(CellSpherical * cell);

public:

	void setSolverTolerance(const double value);

	//! Main cell population
	std::vector<CellSpherical *> cells;

	// boundary elements
	std::vector<ModelElementBarrierTriangle *> barrier;

	// Default constructor
	ModelCellsSpherical();
	virtual ~ModelCellsSpherical();

#pragma region Model interface (Reset, Simulation control)

	// Preliminary init
	virtual void ResetDerived();

	// virtual method from class model (to deprecate InitSimulateInThread).
	// Sets parameters for simulation control
	virtual void SetupSimulation();

	// Method to do multiple simulation steps
	void Simulate(int numberOfSimulationSteps);

#pragma endregion

#pragma region Cell population staining

	//! Stores the last color mode that was called from the gui
	int lastColormodeCells;

	//! Update cells colors using that last color mode
	virtual void UpdateCellsStaining();

	//! Update cell color using a new color mode
	virtual void UpdateCellsStaining(int mode);

	//! Update one cell's color using the last color mode;
	virtual void UpdateCellsStaining(CellSpherical *);

	void setVisible(bool visible);

	virtual CSParameterContext * GetParameters(std::string contextName = "");
	virtual void InitParameters(CSParameterContext * fromContext = 0);

	static void DefaultParameters(CSParameterContext *);

	virtual void writeXML(QXmlStreamWriter *) const;

	static ModelCellsSpherical * createFromXML(QXmlStreamReader *xmlReader,
		std::stringstream & errors,
		std::stringstream & warnings);

	virtual void writeHDF5(H5::H5File * outputFile) const;

	void writeSpheroidData(bool first);

	// method for reading the hdf5 data into an already allocated *Model
	// used after the createFromXML has created a certain Model and initialised
	// the parameters and, most importantly, its name.
	virtual void readModelData(H5::H5File * inputFile,
		std::stringstream & errors,
		std::stringstream & warnings);

	std::vector<double *> mFrictionMatrices;

	// The pointer to the array for all friction matrices created by the
	// cell-cell parallel/perpendicular friction calculation in
	// mpInteractionFrictionMatrix used by the conjugated gradient solver to
	// solve an over-damped equation of motion:
	// This should be 'connected' to the corresponding pointer in
	// mpInteractionFrictionMatrix:
	//   mpFrictionMatrices = mpInteractionFrictionMatrix->mpFrictionMatrices
	// after allocation of mpInteractionFrictionMatrix;
	double * mpFrictionMatrices;

	void writeXML2();

	void spherePacking(double radiusSphere);

	unsigned long mCellsAtLastUpdate;

	chronos::Chronos * timer;

	size_t step;

private:

	virtual void xmlParameters(QXmlStreamReader * xml,
		std::stringstream & warnings,
		std::stringstream & errors);

	virtual void xmlSimulation(QXmlStreamReader * xml,
		std::stringstream & warnings,
		std::stringstream & errors);
};


inline
void
ModelCellsSpherical::UpdateSingleCellEffects()
{
	// MOVE TO C++11
	std::vector<CellSpherical *>::iterator cell;

	for (cell = cells.begin(); cell != cells.end(); ++cell)
	{
		(*cell)->lastForceAbsolute = (*cell)->accumulatedForceAbsolute;
		(*cell)->lastPressure = (*cell)->accumulatedPressure;

		//UpdateDishContactForce(*cell);
		UpdateCellFriction(*cell);
		if (defaultDiffusionConstantCells != 0)
			UpdateForcesLangevin(*cell);

		if ((*cell)->mType == ModelElement::TypeCellSphericalPolar)
			static_cast<CellSphericalPolar *>(*cell)->OrientationMetropolis();

		// Applying forces on cells in dumbbell division phase onto the center of
		// mass of the two division partners.
		// This will be done twice unfortunately, the second time
		// c1Force==c2Force.  It will, though, remove the necessity to distribute
		// every single force contribution onto the division partners.
		if ((*cell)->mpDivisionPartner)
		{
			Vector3f c1Force = (*cell)->directedForce;
			Vector3f c2Force = (*cell)->mpDivisionPartner->directedForce;

			(*cell)->directedForce = 0.5 * (c1Force + c2Force);
			(*cell)->mpDivisionPartner->directedForce = (*cell)->directedForce;

			Vector3f c1Langevin = (*cell)->mLangevinForce;
			Vector3f c2Langevin = (*cell)->mpDivisionPartner->mLangevinForce;

			(*cell)->mLangevinForce = 0.5 * (c1Langevin + c2Langevin);
			(*cell)->mpDivisionPartner->mLangevinForce = (*cell)->mLangevinForce;
		}
	}

	// XXX MOVE THIS SOMEWHERE
	if (mScenario == ScenarioTwoCells)
	{
		if (time > 100 && time < 1000)
		{
			AddDirectedMotion(cells[1]);
			AddDirectedMotion(cells[0]);
		}
	}
}

inline void
ModelCellsSpherical::rescaleLangevinContrib(const double timeStepFactor)
{
	std::vector<CellSpherical *>::iterator cell;

	double timeScale = sqrt(1. / timeStepFactor);
	double scaleFactor = (1 - sqrt(1 / timeStepFactor));
	//std::cout << " timeStepFactor=" << timeStepFactor << " timeScale=" << timeScale
	//    << " scale=" << scaleFactor << std::endl;

	for (cell = cells.begin(); cell != cells.end(); ++cell)
	{
		(*cell)->directedForce.x -= (*cell)->mLangevinForce.x * scaleFactor;
		(*cell)->directedForce.y -= (*cell)->mLangevinForce.y * scaleFactor;
		(*cell)->directedForce.z -= (*cell)->mLangevinForce.z * scaleFactor;

		// update Langevin forces otherwise if we scale again, we get a wrong
		// result
		//std::cout << "Cell: " << (*cell)->mGlobalIndex << " F_L=("
		//    << (*cell)->mLangevinForce.x << " ,"
		//    << (*cell)->mLangevinForce.y << " ,"
		//    << (*cell)->mLangevinForce.z <<")." << std::endl;
		(*cell)->mLangevinForce.x *= timeScale;
		(*cell)->mLangevinForce.y *= timeScale;
		(*cell)->mLangevinForce.z *= timeScale;
	}
}

#endif
