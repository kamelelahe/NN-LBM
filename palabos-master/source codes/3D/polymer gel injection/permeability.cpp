/* This code is developed with the help of Palabos library
  Developed and modified by Elahe Kamel
  kamelElahe@gmail.com

  The most recent release of Palabos can be downloaded at
  <https://palabos.unige.ch/>

 * Contact for Palabos:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 */

 /* This code is used to simulated polymer gel fluid flow in 3D porous media
 polymer gel is a complex fluid with both shear-thining and time-thickening behavior
 if the simulation input parameters for n (power-law index) and TD_factor (time dependency factor) set to 1,
 we can also se this code for simulationg newtonian fluid flow in porous media
 */

#include "palabos3D.h"
#include "palabos3D.hh"

#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>


using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor

//This function calculates pressure gradient 
class PressureGradient {
public:
	PressureGradient(T deltaP_, plint nx_) : deltaP(deltaP_), nx(nx_)
	{ }
	void operator() (plint iX, plint iY, plint iZ, T& density, Array<T, 3>& velocity) const
	{
		velocity.resetToZero();
		density = (T)1 - deltaP * DESCRIPTOR<T>::invCs2 / (T)(nx - 1) * (T)iX;
	}
private:
	T deltaP;
	plint nx;
};

//this part is used to read 2D geometry file with .dat extension.  
void readGeometry(std::string fNameIn, std::string fNameOut, MultiScalarField3D<int>& geometry)
{
	const plint nx = geometry.getNx();
	const plint ny = geometry.getNy();
	const plint nz = geometry.getNz();

	Box3D sliceBox(0, 0, 0, ny - 1, 0, nz - 1);
	std::unique_ptr<MultiScalarField3D<int> > slice = generateMultiScalarField<int>(geometry, sliceBox);
	plb_ifstream geometryFile(fNameIn.c_str());
	for (plint iX = 0; iX < nx - 1; ++iX) {
		if (!geometryFile.is_open()) {
			pcout << "Error: could not open geometry file " << fNameIn << std::endl;
			exit(EXIT_FAILURE);
		}
		geometryFile >> *slice;
		copy(*slice, slice->getBoundingBox(), geometry, Box3D(iX, iX, 0, ny - 1, 0, nz - 1));
	}

	{
		VtkImageOutput3D<T> vtkOut("porousMedium", 1.0);
		vtkOut.writeData<float>(*copyConvert<int, T>(geometry, geometry.getBoundingBox()), "tag", 1.0);
	}

	{
		std::unique_ptr<MultiScalarField3D<T> > floatTags = copyConvert<int, T>(geometry, geometry.getBoundingBox());
		std::vector<T> isoLevels;
		isoLevels.push_back(0.5);
		typedef TriangleSet<T>::Triangle Triangle;
		std::vector<Triangle> triangles;
		Box3D domain = floatTags->getBoundingBox().enlarge(-1);
		domain.x0++;
		domain.x1--;
		isoSurfaceMarchingCube(triangles, *floatTags, isoLevels, domain);
		TriangleSet<T> set(triangles);
		std::string outDir = fNameOut + "/";
		set.writeBinarySTL(outDir + "porousMedium.stl");
	}
}

//set up initial state of porous media
void porousMediaSetup(MultiBlockLattice3D<T, DESCRIPTOR>& lattice,
	OnLatticeBoundaryCondition3D<T, DESCRIPTOR>* boundaryCondition,
	MultiScalarField3D<int>& geometry, T deltaP)
{
	//get dimentions of porous media
	const plint nx = lattice.getNx();
	const plint ny = lattice.getNy();
	const plint nz = lattice.getNz();

	//definition of inlet and outlet condition and location
	//inlet
	pcout << "Definition of inlet/outlet." << std::endl;
	Box3D inlet(0, 0, 1, ny - 2, 1, nz - 2);
	boundaryCondition->addPressureBoundary0N(inlet, lattice);
	setBoundaryDensity(lattice, inlet, (T)1.);
	//outlet
	Box3D outlet(nx - 1, nx - 1, 1, ny - 2, 1, nz - 2);
	boundaryCondition->addPressureBoundary0P(outlet, lattice);
	setBoundaryDensity(lattice, outlet, (T)1. - deltaP * DESCRIPTOR<T>::invCs2);

	//definition of inner boundaries
	pcout << "Definition of the geometry." << std::endl;
	// Where "geometry" evaluates to 1, use bounce-back.
	defineDynamics(lattice, geometry, new BounceBack<T, DESCRIPTOR>(), 1);
	// Where "geometry" evaluates to 2, use no-dynamics (which does nothing).
	defineDynamics(lattice, geometry, new NoDynamics<T, DESCRIPTOR>(), 2);

	//initialize porous media at equilibrium state
	pcout << "Initilization of rho and u." << std::endl;
	initializeAtEquilibrium(lattice, lattice.getBoundingBox(), PressureGradient(deltaP, nx));

	lattice.initialize();
	delete boundaryCondition;
}

//this function is exactly as previous one. but instead of inilitizing porous media at equilibrium state, we keep the previous state of it
void NNporousMediaSetup(MultiBlockLattice3D<T, DESCRIPTOR>& lattice,
	OnLatticeBoundaryCondition3D<T, DESCRIPTOR>* boundaryCondition,
	MultiScalarField3D<int>& geometry, T deltaP)
{
	//get dimentions of porous media
	const plint nx = lattice.getNx();
	const plint ny = lattice.getNy();
	const plint nz = lattice.getNz();

	//definition of inlet and outlet condition and location
	//inlet
	pcout << "Definition of inlet/outlet." << std::endl;
	Box3D inlet(0, 0, 1, ny - 2, 1, nz - 2);
	boundaryCondition->addPressureBoundary0N(inlet, lattice);
	setBoundaryDensity(lattice, inlet, (T)1.);

	//outlet
	Box3D outlet(nx - 1, nx - 1, 1, ny - 2, 1, nz - 2);
	boundaryCondition->addPressureBoundary0P(outlet, lattice);
	setBoundaryDensity(lattice, outlet, (T)1. - deltaP * DESCRIPTOR<T>::invCs2);

	//definition of inner boundaries
	pcout << "Definition of the geometry." << std::endl;
	// Where "geometry" evaluates to 1, use bounce-back.
	defineDynamics(lattice, geometry, new BounceBack<T, DESCRIPTOR>(), 1);
	// Where "geometry" evaluates to 2, use no-dynamics (which does nothing).
	defineDynamics(lattice, geometry, new NoDynamics<T, DESCRIPTOR>(), 2);

	lattice.initialize();
	delete boundaryCondition;
}


//this function generates output image for omega field in gif format
void writeGifs(MultiBlockLattice3D<T, DESCRIPTOR>& lattice, plint iter)
{
	const plint nx = lattice.getNx();
	const plint ny = lattice.getNy();
	const plint nz = lattice.getNz();

	const plint imSize = 600;
	ImageWriter<T> imageWriter("leeloo");

	imageWriter.writeScaledGif(createFileName("omega_half", iter, 6), *computeOmega(lattice, Box3D(nx / 2, nx / 2, 0, ny - 1, 0, nz - 1)),
		imSize, imSize);
	imageWriter.writeScaledGif(createFileName("ux_half", iter, 6),
		*computeVelocityNorm(lattice, Box3D(nx / 2, nx / 2, 0, ny - 1, 0, nz - 1)),
		imSize, imSize);
	imageWriter.writeScaledGif(createFileName("omega_ydir", iter, 6), *computeOmega(lattice, Box3D(0, nx - 1, ny / 2, ny / 2, 0, nz - 1)),
		imSize, imSize);
	imageWriter.writeScaledGif(createFileName("ux_ydir", iter, 6),
		*computeVelocityNorm(lattice, Box3D(0, nx - 1, ny / 2, ny / 2, 0, nz - 1)),
		imSize, imSize);

}

//write vtk files for simulation output
void writeVTK(MultiBlockLattice3D<T, DESCRIPTOR>& lattice, plint iter)
{
	VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), 1.);
	vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", 1.);
	vtkOut.writeData<3, float>(*computeVelocity(lattice), "velocity", 1.);
	vtkOut.writeData< float>(*computeOmega(lattice), "omega", 1.);
	vtkOut.writeData<float>(*computeDensity(lattice), "density", 1.);
	vtkOut.writeData<6, float>(*computePiNeq(lattice), "piNeq", 1.);
	vtkOut.writeData<6, float>(*computeShearStress(lattice), "ShearStress ", 1.);
}

//this function calculates the permeability of porous media
T computePermeability(MultiBlockLattice3D<T, DESCRIPTOR>& lattice, T nu, T deltaP, Box3D domain)
{
	pcout << "Computing the permeability." << std::endl;

	// Compute only the x-direction of the velocity (direction of the flow).
	plint xComponent = 0;
	plint nx = lattice.getNx();

	T meanU = computeAverage(*computeVelocityComponent(lattice, domain, xComponent));

	pcout << "Average velocity     = " << meanU << std::endl;
	pcout << "Lattice viscosity nu = " << nu << std::endl;
	pcout << "Grad P               = " << deltaP / (T)(nx - 1) << std::endl;
	pcout << "Permeability         = " << nu * meanU / (deltaP / (T)(nx - 1)) << std::endl;

	return meanU;
}

int main(int argc, char** argv)
{
	plbInit(&argc, &argv);


	//for this simulation, 12 input values re required as below
	if (argc != 12) {
		pcout << "Error missing some input parameter\n";
		pcout << "The structure is :\n";
		pcout << "1. Input file name.\n";
		pcout << "2. Output directory name.\n";
		pcout << "3. number of cells in X direction.\n";
		pcout << "4. number of cells in Y direction.\n";
		pcout << "5. number of cells in Z direction.\n";
		pcout << "6. Delta P .\n";
		pcout << "7. n (power-law exponent) .\n";
		pcout << "8. Nu 0 .\n";
		pcout << "9. carreau number .\n";
		pcout << "10. Nu inf .\n";
		pcout << "11. time dependancy factor .\n";
		pcout << "Example: " << argv[0] << " benheimer.dat tmp/ 300 300 300 0.0005 1 10 1 1\n";
		exit(EXIT_FAILURE);
	}

	//assign user input variables to  local variables
	std::string fNameIn = argv[1];
	std::string fNameOut = argv[2];
	const plint nx = atoi(argv[3]);
	const plint N = atoi(argv[3]);
	const plint ny = atoi(argv[4]);
	const plint nz = atoi(argv[5]);
	const T deltaP = atof(argv[6]);
	const T n = atof(argv[7]);
	const T Nu0 = atof(argv[8]);
	const T Cu = atof(argv[9]);
	const T NuInf = atof(argv[10]);
	const T timeDepend = atof(argv[11]);

	global::directories().setOutputDir(fNameOut + "/");

	//calculate dynamic variables
	const T tau = Nu0 * DESCRIPTOR<T>::invCs2 + (T)0.5;
	const T omega = (T)1 / tau;
	const T nu = Nu0;
	MultiTensorField3D<T, 2> velField(nx, ny, nz);

	pcout << "Creation of the lattice." << std::endl;
	MultiBlockLattice3D<T, DESCRIPTOR> lattice(nx, ny, nz, new BGKdynamics<T, DESCRIPTOR>(omega));

	// Switch off periodicity.
	lattice.periodicity().toggleAll(false);

	pcout << "Reading the geometry file." << std::endl;
	MultiScalarField3D<int> geometry(nx, ny, nz);
	readGeometry(fNameIn, fNameOut, geometry);

	pcout << "nu = " << nu << std::endl;
	pcout << "tau = " << tau << std::endl;
	pcout << "omega = " << omega << std::endl;
	pcout << "deltaP = " << deltaP << std::endl;
	pcout << "nx = " << lattice.getNx() << std::endl;
	pcout << "ny = " << lattice.getNy() << std::endl;
	pcout << "nz = " << lattice.getNz() << std::endl;

	//set up porous media before starting first part of simulation(newtonian fluid)
	porousMediaSetup(lattice, createLocalBoundaryCondition3D<T, DESCRIPTOR>(), geometry, deltaP);

	//starting fisrt part of simulation. newtonian fluid flow.
	//we did this part because in order to run non-Newtonian fluid flow in porous media, maximum velocity is required.
	//so this simulation provides us with rough stimation of average velocity
	pcout << "Simulation begins" << std::endl;

	T new_avg_f, old_avg_f, relE_f1;
	const plint maxT = 10000;
	T meanVel, meanP, meanRho;
	plint xComponent = 0;
	plint iT = 0;
	T conv = 10e-12;
	bool continueSimulation = true;
	lattice.toggleInternalStatistics(false);

	for (; iT < maxT; ++iT) {

		lattice.toggleInternalStatistics(true);
		//collision and streaming
		lattice.collideAndStream();
		new_avg_f = getStoredAverageEnergy(lattice);
		lattice.toggleInternalStatistics(false);

		relE_f1 = std::fabs(old_avg_f - new_avg_f) * 100 / old_avg_f;

		/*if (iT % 5000 == 0 ) {
			writeVTK(lattice,iT);

		}
		if (iT % 5000 == 0) {
			writeGifs(lattice,iT);

		}*/

		//if convergence crieteria met, stop simulation
		if (relE_f1 < conv) {
			break;
		}

		old_avg_f = new_avg_f; // store new properties
		if (iT % 5000 == 0 && iT > 0) {
			pcout << "Iteration " << iT << std::endl;
			pcout << "-----------------" << std::endl;
			pcout << "Relative difference of Energy: " << std::setprecision(3) << relE_f1 << " %" << std::endl;
			pcout << "The preliminary permeability is: " << std::endl;
			computePermeability(lattice, nu, deltaP, lattice.getBoundingBox());
			pcout << "**********************************************" << std::endl;
		}
	}

	pcout << "starting the simulation of non-Newtonian fluid flow" << std::endl;
	//compute maximumm velocity after  simulation
	const T uMax = computeMax(*computeVelocityComponent(lattice, lattice.getBoundingBox(), xComponent));;
	pcout << "max vel is " << uMax << std::endl;
	//calculate reynolds number from maximum velocity
	T Re = (uMax * N) / Nu0;

	CarreauFlowParam<T> parameters(
		(T)uMax,    // uMax
		(T)Re,    // Re
		(T)Cu,     // Cu
		(T)NuInf,      // NuInf (Only the case nuInf = 0 has been implemented for the velocity inlet bc)
		(T)n,      // n
		N,          // N
		(T)timeDepend,            //time dependency factor
		1.,         // lx
		1.,          // ly
		1.           // lz

	);
	//write a log file for carreau flow properties
	writeLogFile(parameters, "Carreau Poseuille Flow");

	//set carreau parameters
	global::CarreauParameters().setNu0(parameters.getLatticeNu0());
	global::CarreauParameters().setNuInf(parameters.getLatticeNuInf());
	global::CarreauParameters().setLambda(parameters.getLatticeLambda());
	global::CarreauParameters().setExponent(parameters.getExponent());
	global::CarreauParameters().setTDfactor(parameters.getTDfactor());

	//Set carreau dynamic over the domain
	setCompositeDynamics(
		lattice,
		lattice.getBoundingBox(),
		new CarreauDynamics<T, DESCRIPTOR, 1>(new BGKdynamics<T, DESCRIPTOR>(omega)));

	//set up porous media for non-Newtonaian geometry
	NNporousMediaSetup(lattice, createLocalBoundaryCondition3D<T, DESCRIPTOR>(), geometry, deltaP);

	const plint maxTT = 30001;
	//starting simualtion of polymer gel fluid flow in porous media
	for (; iT < maxTT; ++iT) {

		lattice.toggleInternalStatistics(true);
		lattice.collideAndStream();
		new_avg_f = getStoredAverageEnergy(lattice);
		lattice.toggleInternalStatistics(false);

		relE_f1 = std::fabs(old_avg_f - new_avg_f) * 100 / old_avg_f;
		if (iT % 30000 == 0 && iT > 0) {
			writeVTK(lattice, iT);
			writeGifs(lattice, iT);
		}

		if (relE_f1 < conv) {
			break;
		}

		old_avg_f = new_avg_f; // store new properties
		if (iT % 5000 == 0 && iT > 0) {
			pcout << "Iteration " << iT << std::endl;
			pcout << "-----------------" << std::endl;
			pcout << "Relative difference of Energy: " << std::setprecision(3) << relE_f1 << " %" << std::endl;
			pcout << "Average Omega is " << computeAverage(*computeOmega(lattice, lattice.getBoundingBox())) << std::endl;

			pcout << "**********************************************" << std::endl;
		}
	}

	pcout << "End of simulation at iteration " << iT << std::endl;

	//compute important parameters
	meanVel = computeAverage(*computeVelocityComponent(lattice, lattice.getBoundingBox(), xComponent));
	meanRho = computeAverage(*computeDensity(lattice, lattice.getBoundingBox()));
	T meanOmega = computeAverage(*computeOmega(lattice, lattice.getBoundingBox()));
	T meanTau = (T)1 / meanOmega;
	T meanVisco = ((T)2 * meanTau - (T)1) / (T)6;
	meanP = meanRho * DESCRIPTOR<T>::invCs2;
	std::string fullName = global::directories().getLogOutDir() + "PMLog.dat";

	//write down results in a file
	plb_ofstream ofile(fullName.c_str(), continueSimulation ? std::ofstream::app : std::ofstream::out);

	ofile << "porous Media" << "\n\n";
	ofile << "nx = " << lattice.getNx() << "\n";
	ofile << "ny = " << lattice.getNy() << "\n";
	ofile << "deltaP = " << deltaP << "\n";
	ofile << "Average nu = " << meanVisco << "\n";
	ofile << "Average omega = " << meanOmega << "\n";
	ofile << "Tau = " << meanTau << "\n";
	ofile << "Average velocity = " << meanVel << "\n";
	ofile << "mean Desnity= " << meanRho << "\n";
	ofile << "mean pressure= " << meanP << "\n";
	ofile << "Grad P = " << deltaP / (T)(nx - 1) << "\n";
	ofile << "Permeability = " << meanVisco * meanVel / (deltaP / (T)(nx - 1)) << "\n";
	
	pcout << "Permeability:" << std::endl << std::endl;
	computePermeability(lattice, nu, deltaP, lattice.getBoundingBox());
	pcout << std::endl;

	//pcout << "Writing VTK file ..." << std::endl << std::endl;
	//writeVTK(lattice,iT);
	pcout << "Finished!" << std::endl << std::endl;

	return 0;
}
