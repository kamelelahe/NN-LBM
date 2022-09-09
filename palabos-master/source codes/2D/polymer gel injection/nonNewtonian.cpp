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

/* This code is used to simulated polymer gel fluid flow in 2D porous media
polymer gel is a complex fluid with both shear-thining and time-thickening behavior
if the simulation input parameters for n (power-law index) and TD_factor (time dependency factor) set to 1,
we can also se this code for simulationg newtonian fluid flow in porous media*/

#include "palabos2D.h"
#include "palabos2D.hh"

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
#define DESCRIPTOR descriptors::D2Q9Descriptor

//This function calculates pressure gradient 
class PressureGradient {
public:
	PressureGradient(T deltaP_, plint nx_) : deltaP(deltaP_), nx(nx_)
	{ }
	void operator() (plint iX, plint iY, T& density, Array<T, 2>& velocity) const
	{
		velocity.resetToZero();
		density = (T)1 - deltaP * DESCRIPTOR<T>::invCs2 / (T)(nx - 1) * (T)iX;

	}
private:
	T deltaP;
	plint nx;
};

//this part is used to read 2D geometry file with .dat extension.  
void readGeometry(std::string fNameIn, std::string fNameOut, MultiScalarField2D<int>& geometry)
{
	const plint nx = geometry.getNx();
	const plint ny = geometry.getNy();

	Box2D sliceBox(0, 0, 0, ny - 1);
	std::unique_ptr<MultiScalarField2D<int> > slice = generateMultiScalarField<int>(geometry, sliceBox);
	plb_ifstream geometryFile(fNameIn.c_str());
	for (plint iX = 0; iX < nx - 1; ++iX) {
		if (!geometryFile.is_open()) {
			pcout << "Error: could not open geometry file " << fNameIn << std::endl;
			exit(EXIT_FAILURE);
		}
		geometryFile >> *slice;
		copy(*slice, slice->getBoundingBox(), geometry, Box2D(iX, iX, 0, ny - 1));
	}

}

//set up initial state of porous media
void porousMediaSetup(MultiBlockLattice2D<T, DESCRIPTOR>& lattice,
	OnLatticeBoundaryCondition2D<T, DESCRIPTOR>* boundaryCondition,
	MultiScalarField2D<int>& geometry, T deltaP)
{
	//get dimentions of porous media
	const plint nx = lattice.getNx();
	const plint ny = lattice.getNy();

	//definition of inlet and outlet condition and location
	//inlet
	pcout << "Definition of inlet/outlet." << std::endl;
	Box2D inlet(0, 0, 1, ny - 2);
	boundaryCondition->addPressureBoundary0N(inlet, lattice);
	setBoundaryDensity(lattice, inlet, (T)1.);
	//outlet
	Box2D outlet(nx - 1, nx - 1, 1, ny - 2);
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
void NNporousMediaSetup(MultiBlockLattice2D<T, DESCRIPTOR>& lattice,
	OnLatticeBoundaryCondition2D<T, DESCRIPTOR>* boundaryCondition,
	MultiScalarField2D<int>& geometry, T deltaP)
{
	//get dimentions of porous media
	const plint nx = lattice.getNx();
	const plint ny = lattice.getNy();

	//definition of inlet and outlet condition and location
	//inlet
	pcout << "Definition of inlet/outlet." << std::endl;
	Box2D inlet(0, 0, 1, ny - 2);
	boundaryCondition->addPressureBoundary0N(inlet, lattice);
	setBoundaryDensity(lattice, inlet, (T)1.);
	//outlet
	Box2D outlet(nx - 1, nx - 1, 1, ny - 2);
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
void writeGifs(MultiBlockLattice2D<T, DESCRIPTOR>& lattice, plint iter)
{
	const plint nx = lattice.getNx();
	const plint ny = lattice.getNy();

	const plint imSize = 600;
	ImageWriter<T> imageWriter("leeloo");

	// Write velocity-norm at x=0.
	imageWriter.writeScaledGif(createFileName("omega", iter, 6), *computeOmega(lattice, Box2D(0, nx - 1, 0, ny - 1)),
		imSize, imSize);
	imageWriter.writeScaledGif(createFileName("ux_inlet", iter, 6),
		*computeVelocityNorm(lattice, Box2D(0, nx - 1, 0, ny - 1)),
		imSize, imSize);
}

//This function writes simlation result of NormVelocity, NormShearRate and Omega in .dat file
std::string createFileNameDat(plint number, plint width, std::string name) {
	std::stringstream fNameStream;
	fNameStream << std::setfill('0') << std::setw(width) << number << name;
	return fNameStream.str();
}
void writeDat(std::string fNameOut, MultiBlockLattice2D<T, DESCRIPTOR>& lattice, plint iter)
{

	//generate unique name for Norm Velocity .dat file
	std::string Myname_v = createFileNameDat(iter, 6, "NormVelocity.dat");
	const char* fname_v = Myname_v.c_str();
	char result_v[100];
	const char* nname = fNameOut.c_str();
	strcpy(result_v, nname);
	strcat(result_v, fname_v);
	//write norm velocity in simulation step on generated file
	plb_ofstream fout(result_v);
	fout << *computeVelocity(lattice);
	fout.close();

	//generate unique name for Norm ShearRate .dat file
	std::string Myname_S = createFileNameDat(iter, 6, "NormShearRate.dat");
	const char* fname_S = Myname_S.c_str();
	char result_S[100];
	strcpy(result_S, nname);
	strcat(result_S, fname_S);
	//write Norm ShearRate in simulation step on generated file
	fout.open(result_S);
	fout << *computeNorm(*computeShearStress(lattice));
	fout.close();

	//generate unique name for Norm Omega .dat file
	std::string Myname_O = createFileNameDat(iter, 6, "Omega.dat");
	const char* fname_O = Myname_O.c_str();
	char result_O[100];
	strcpy(result_O, nname);
	strcat(result_O, fname_O);
	//write Omega in simulation step on generated file
	fout.open(result_O);
	fout << *computeOmega(lattice);
	fout.close();

}

//This function generates vtk file for velocity norm and velocity
void writeVTK(MultiBlockLattice2D<T, DESCRIPTOR>& lattice, plint iter)
{
	VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), 1.);
	vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", 1.);
	vtkOut.writeData<2, float>(*computeVelocity(lattice), "velocity", 1.);
}

//this function calculates the permeability od porous media
T computePermeability(MultiBlockLattice2D<T, DESCRIPTOR>& lattice, T nu, T deltaP, Box2D domain)
{
	pcout << "Computing the permeability." << std::endl;

	// Compute only the x-direction of the velocity (direction of the flow).
	plint xComponent = 0;
	plint nx = lattice.getNx();

	//calculate average value of velocity
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

	//for this simulation, 11 input values re required as below
	if (argc != 11) {
		pcout << "Error missing some input parameter\n";
		pcout << "The structure is :\n";
		pcout << "1. Input file name.\n";
		pcout << "2. Output directory name.\n";
		pcout << "3. number of cells in X direction.\n";
		pcout << "4. number of cells in Y direction.\n";
		pcout << "5. Delta P .\n";
		pcout << "6. n (power-law exponent) .\n";
		pcout << "7. Nu 0 .\n";
		pcout << "8. carreau number .\n";
		pcout << "9. Nu inf .\n";
		pcout << "10. time dependancy factor .\n";
		pcout << "Example: " << argv[0] << " benheimer.dat tmp/ 300 300 0.0005 1 10 1 1\n";
		exit(EXIT_FAILURE);
	}
	std::string fNameIn = argv[1];
	std::string fNameOut = argv[2];

	//assign user input variables to  local variables
	const plint nx = atoi(argv[3]);
	const plint N = atoi(argv[3]);
	const plint ny = atoi(argv[4]);
	const T deltaP = atof(argv[5]);
	const T n = atof(argv[6]);
	const T Nu0 = atof(argv[7]);
	const T Cu = atof(argv[8]);
	const T NuInf = atof(argv[9]);
	const T timeDepend = atof(argv[10]);

	global::directories().setOutputDir(fNameOut + "/");

	//calculate dynamic variables
	const T tau = Nu0 * DESCRIPTOR<T>::invCs2 + (T)0.5;
	const T omega = (T)1 / tau;
	const T nu = Nu0;
	MultiTensorField2D<T, 2> velField(nx, ny);

	pcout << "Creation of the lattice." << std::endl;
	MultiBlockLattice2D<T, DESCRIPTOR> lattice(nx, ny, new BGKdynamics<T, DESCRIPTOR>(omega));

	// Switch off periodicity.
	lattice.periodicity().toggleAll(false);

	pcout << "Reading the geometry file." << std::endl;
	MultiScalarField2D<int> geometry(nx, ny);
	readGeometry(fNameIn, fNameOut, geometry);

	pcout << "nu = " << nu << std::endl;
	pcout << "tau = " << tau << std::endl;
	pcout << "omega = " << omega << std::endl;
	pcout << "deltaP = " << deltaP << std::endl;
	pcout << "nx = " << lattice.getNx() << std::endl;
	pcout << "ny = " << lattice.getNy() << std::endl;

	//set up porous media before starting first part of simulation(newtonian fluid)
	porousMediaSetup(lattice, createLocalBoundaryCondition2D<T, DESCRIPTOR>(), geometry, deltaP);

	//starting fisrt part of simulation. newtonian fluid flow.
	//we did this part because in order to run non-Newtonian fluid flow in porous media, maximum velocity is required.
	//so this simulation provides us with rough stimation of average velocity
	pcout << "Simulation begins" << std::endl;
	T new_avg_f, old_avg_f, relE_f1;
	const plint maxT = 10001;
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
		//execute writeDate function whch writes some simulation result id dat file every 1000 iteration
		if (iT % 10000 == 0 && iT > 0) {
			writeDat(fNameOut, lattice, iT);

		}
		//if convergence crieteria met, stop simulation
		if (relE_f1 < conv) {
			break;
		}

		old_avg_f = new_avg_f;

		//console out result every 10000 iterations
		if (iT % 10000 == 0 && iT > 0) {
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

	//store carreau flow properties 
	CarreauFlowParam<T> parameters(
		(T)uMax,    // uMax
		(T)Re,    // Re
		(T)Cu,     // Cu
		(T)NuInf,      // NuInf (Only the case nuInf = 0 has been implemented for the velocity inlet bc)
		(T)n,      // n
		N,          // N
		(T)timeDepend,         //time dependency factor
		1.,         // lx
		1.          // ly

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
	NNporousMediaSetup(lattice, createLocalBoundaryCondition2D<T, DESCRIPTOR>(), geometry, deltaP);

	const plint maxTT = 30001;
	//starting simualtion of polymer gel fluid flow in porous media
	for (; iT < maxTT; ++iT) {

		lattice.toggleInternalStatistics(true);
		lattice.collideAndStream();
		new_avg_f = getStoredAverageEnergy(lattice);
		lattice.toggleInternalStatistics(false);

		relE_f1 = std::fabs(old_avg_f - new_avg_f) * 100 / old_avg_f;
		if (iT % 10000 == 0 && iT > 0) {
			writeDat(fNameOut, lattice, iT);

		}
		if (relE_f1 < conv) {
			break;
		}

		old_avg_f = new_avg_f; // store new properties
		if (iT % 10000 == 0 && iT > 0) {
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
	Box2D Myoutlet(nx - 1, nx - 1, 1, ny - 2);
	T q1 = computeSum(*computeVelocityNorm(lattice, Myoutlet));
	T UMeanOutlet = computeAverage(*computeVelocityNorm(lattice, Myoutlet));
	//T q2=computeSum(*computeVelocity(lattice,Myoutlet));
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
	ofile << "Average velocity Total = " << meanVel << "\n";
	ofile << "Average Desnity= " << meanRho << "\n";
	ofile << "Average pressure= " << meanP << "\n";
	ofile << "Grad P = " << deltaP / (T)(nx - 1) << "\n";
	ofile << "q Norm = " << q1 << "\n";
	ofile << "average velocity outlet = " << UMeanOutlet << "\n";;
	ofile << "Permeability = " << meanVisco * meanVel / (deltaP / (T)(nx - 1)) << "\n";
	ofile << "max Norm shear-rate" << computeMax(*computeNorm(*computeShearStress(lattice))) << "\n";
	ofile << "min Norm shear-rate" << computeMin(*computeNorm(*computeShearStress(lattice))) << "\n";
	ofile << "max Omega" << computeMax(*computeOmega(lattice)) << "\n";
	ofile << "min Omega" << computeMin(*computeOmega(lattice)) << "\n";
	ofile << "max velocityNorm" << computeMax(*computeVelocityNorm(lattice)) << "\n";
	ofile << "min velocityNorm" << computeMin(*computeVelocityNorm(lattice)) << "\n";

	pcout << "Permeability:" << std::endl << std::endl;
	computePermeability(lattice, nu, deltaP, lattice.getBoundingBox());
	pcout << std::endl;
	writeGifs(lattice, iT);

	pcout << "Finished!" << std::endl << std::endl;

	return 0;
}
