#pragma once

#include <vector>
#include <map>
#include <string>
#include "Matrix.h"
#include "rng.h"
#include "ConfigFile.h"
using std::map;
using std::vector;
using std::string;

class Simulation {

public:
	Simulation(const ConfigFile& config);

	~Simulation();

	void perform();

	void initialization();

	void colonization();

	void colonizationForest();

	void invasiveColonization();

	//---------------------------------------------------------------------------
	static float z01()   // float
	{
		return(float(rand())/(RAND_MAX+1));     //32767.0
	}

	//---------------------------------------------------------------------------
	float z02()
	{
		return (float)Zufall.rand_halfclosed01();
	}
	//----------------------------------------------------------------------------

protected:
	//! Initialize DX, DY as well as Disp
	void readDispersalMatrix(const std::string& dispFileName);

	//! Initialize NX, NY
	void readMap(const std::string& mapFileName);

	void writeParameterSet(vector<double>& parameterSet,vector<double>& minValue,vector<double>& maxValue,vector<double>& allometabolConst);

	void addSpeciesIfNew(vector<double>& parameterSet, map<vector<double>, Matrix<float> >& spInfo,map<vector<double>, 
		Matrix<float> >& spDisp, map<vector<double>, Matrix<float> >& spHabt, map<vector<double>, Matrix<float> >& spDispMain,
								 map<vector<double>, Matrix<float> >& spHabtMain);

	void addRandomSpeciesIfNew(vector<double>& parameterSet, map<vector<double>, Matrix<float> >& spInfo, map<vector<double>, 
		Matrix<float> >& spDisp, map<vector<double>, Matrix<float> >& spHabt);

	void createDispersalKernel(vector<double>& parameterSet, map<vector<double>, Matrix<float> >& spDisp);

	void habitatSuitability(vector<double>& parameterSet, map<vector<double>, Matrix<float> >& spHabt);

	double dispersalNegExpFunction(double alpha, double r);

	double dispersal2DtFunction(double alpha, double P, double r);

	double gamma (double x);
		
	void dispersalfromMainland(map<vector<double>, Matrix<float> >& spInfo, map<vector<double>, Matrix<float> >& spDisp, 
		Matrix<float>& K, map<vector<double>, Matrix<float> >& spHabt);

	bool matrixIsEmpty(Matrix<float>& matrix);

	void removeNoncolonizingSpecies(map<vector<double>, Matrix<float> >& spInfo, map<vector<double>, Matrix<float> >& spDisp,
								 map<vector<double>, Matrix<float> >& spHabt);

	void removeExtinctSpecies(map<vector<double>, Matrix<float> >& spAdult, 
	map<vector<double>, Matrix<float> >& spRecruit,	map<vector<double>, Matrix<float> >& spSeed,
	map<vector<double>, Matrix<float> >& spDisp, map<vector<double>, Matrix<float> >& spHabt);

	void addSpeciesIfcolonizing(map<vector<double>, Matrix<float> >& spInfoCol, map<vector<double>, Matrix<float> >& spAdult, 
		map<vector<double>, Matrix<float> >& spRecruit, map<vector<double>, Matrix<float> >& spSeed, 
		Matrix<float>& K,map<vector<double>, Matrix<float> >& spDispCol, map<vector<double>, Matrix<float> >& spDisp,
		map<vector<double>, Matrix<float> >& spHabtCol,map<vector<double>, Matrix<float> >& spHabt);

	void LookforCladoSpecies(map<vector<double>, Matrix<float> >& spAdult, 
	map<vector<double>, Matrix<float> >& spRecruit,	map<vector<double>, Matrix<float> >& spSeed,
	map<vector<double>, Matrix<float> >& spDisp, map<vector<double>, Matrix<float> >& spHabt);
	
	void LookforAnaSpecies(map<vector<double>, Matrix<float> >& spAdult, 
	map<vector<double>, Matrix<float> >& spRecruit,	map<vector<double>, Matrix<float> >& spSeed,
	map<vector<double>, Matrix<float> >& spDisp, map<vector<double>, Matrix<float> >& spHabt, map<vector<double>, Matrix<float> >& immGeneflow);
	
	//reproduction, dispersal, speciation and recruitment:
	void localdynamics(map<vector<double>, Matrix<float> >& spAdult, map<vector<double>, Matrix<float> >& spRecruit,
					map<vector<double>, Matrix<float> >& spSeed, Matrix<float>& K, 
		map<vector<double>, Matrix<float> >& spDisp,map<vector<double>, Matrix<float> >& spHabt, map<vector<double>, Matrix<float> >& spHabtMain);

	//reproduction, dispersal, speciation and recruitment with random parameter combination (Allometabolic = 0):
	void localdynamicsNoMTE(map<vector<double>, Matrix<float> >& spAdult, map<vector<double>, Matrix<float> >& spRecruit,
					map<vector<double>, Matrix<float> >& spSeed, Matrix<float>& K,
	               map<vector<double>, Matrix<float> >& spDisp, map<vector<double>, Matrix<float> >& spHabt, map<vector<double>, Matrix<float> >& spHabtMain);

	//The island dynamics AND the output:
	void islanddynamics(map<vector<double>, Matrix<float> >& spAdult, map<vector<double>, Matrix<float> >& spRecruit,
				map<vector<double>, Matrix<float> >& spSeed, Matrix<float>& K, 
		map<vector<double>, Matrix<float> >& spDisp, map<vector<double>, Matrix<float> >& spHabt, map<vector<double>, Matrix<float> >& spHabtMain, Matrix<float>& TempKelvinBareIsland);

	/** Parameters of the whole simulation **/

	unsigned int DY;   //# rows in dispersal matrix
	unsigned int DX;   //# columns in dispersal matrix

	unsigned int NY;   //# rows in landscape matrix
	unsigned int NX;   //# columns in landscape matrix

	//defining Pi:
	double PI;

	//dispersal radius:
	unsigned int kernelRadius;
	unsigned int replicates;
	unsigned int simulationTimesteps;
	unsigned int islandTimestep;
	unsigned int migrationEvents;
	unsigned int maxIslandRadius;
	unsigned int speciesOutputtimestep;
	std::string parameterRangesTree;
	std::string parameterRangesShrub;
	std::string parameterRangesHerb;
	string simulatedPFT;
	string dispKernel;
	string habSuitability;	
	string readMainlandSpeciespool;
	string treeTemperatureFeedback; 
	string colonize;
	string singleSpecies;
	string metabolicRateToCorrect; 
	string cladospeciation;
	string anaspeciation;
	string treeline;
	string altitudinalTemp;
	string boundaries;
	string vulnerableSpecialists;
	string fixedIslandaltitude;
	string dynamicGeology;
	unsigned int mainlandSpeciesNumber;
	unsigned int invasiveSpeciesNumber; 
	unsigned int randomInvasivecenters;
	unsigned int randomDisturbance, vulcanicDisturbance;
	unsigned int invasivesBegin, invasivesperTimestep, speciesAbundDistriBegin; 
	unsigned int climateChange, temperatureChangeInterval;
	float temperatureChange; 
	float pftTree, pftShrub, pftHerb;
	float disturbanceProb;
	float allometabolic;
	float temperature;
	float tempForFixedMetabolism;
	float Boltzmann;
	float E;
	float maxK;
    unsigned int verticalpreference,  lowlandDegradation,  degradationImpact, mountainRadius;
	int maxOptaltitude, minOptaltitude;
	float ldd;
	float nicheEvolution;

	Matrix<double> Disp; // dispersal matrix
	Matrix<int> State; // landscape matrix
	Matrix<int> gisdata2; // map containing the initial states

	// Random generator
	KW_RNG::RNG Zufall;

	// Gamma Function generator
	KW_RNG::RNG Gamma;

	/** Parameters of the current run **/
	// coordinates list for all island patches
	// is initialized in constructor
	std::vector<std::pair<unsigned int,unsigned int> > patchesPosition;

	// coordinates list for all Mainland patches
	// is initialized in constructor
	std::vector<std::pair<unsigned int,unsigned int> > patchesPositionMainLand;

	//interactions/model timestep:
	unsigned int timeStep;
	unsigned int replicate;


	//
	unsigned int growthStep;
	unsigned int erosionStep;
	unsigned int erosionStep2;
	unsigned int centralx;
	unsigned int centraly;

	//if we use logs for parameter values
	std::vector<bool> useLogForParam;

	//the min, max values are stored
	std::vector<double> minValueTree, maxValueTree, allometabolConstTree, minValueShrub, maxValueShrub, allometabolConstShrub, minValueHerb, maxValueHerb,allometabolConstHerb;

	//Carrying capacity matrix:
	Matrix<float> Ktree;

	//Dynamic teperature matrix:
	Matrix<float> TempKelvinBareIsland; //dynamic with island dynamics
	Matrix<float> TempKelvin; //dynamic with island dynamics and tree cover

	//For Habitat suitability:
	Matrix<int> checkedSuitability; 

	//For lowland degradation:
	Matrix<float> lowlandDegradationMatrix;

	//Species map for colonization from mainland: abundance matrices with parameters as key:
	map<vector<double>, Matrix<float> > speciesInfoCol;

	//map of sp-specific dispersal kernels for colonization from mainland: abundance matrices with parameters as key:
	map<vector<double>, Matrix<float> > speciesDispCol;

	//map of sp-specific suitable habitat for colonization from mainland: abundance matrices with parameters as key:
	map<vector<double>, Matrix<float> > speciesHabtCol;

	//Species map for the islands: abundance matrices with parameters as key:
	map<vector<double>, Matrix<float> > speciesAdult;
	//Species map for the islands: abundance matrices with parameters as key:
	map<vector<double>, Matrix<float> > speciesRecruit;
	//Species map for the islands: abundance matrices with parameters as key:
	map<vector<double>, Matrix<float> > speciesSeed;

	//map of sp-specific dispersal kernels for the islands: abundance matrices with parameters as key:
	map<vector<double>, Matrix<float> > speciesDisp;
	
	//map of sp-specific suitable habitat for colonization from mainland: abundance matrices with parameters as key:
	map<vector<double>, Matrix<float> > speciesHabt;

	//map of sp-specific gene flow, for protracted anagenesis (Rosindell & Phillimore 2011) :
	map<vector<double>, Matrix<float> > immGeneflow;

	//vector of all species from the mainland, defined as a vector of species specific parameters: 
	//std::vector<std::pair<unsigned int,std::vector<double>>> mainlandSpeciesPool;
	std::vector<std::vector<double>> mainlandSpeciesPool;

	//vector of all invasive species (same number of species as mainland), defined as a vector of species specific parameters: 
	std::vector<std::vector<double>> invasiveSpeciesPool;

	//map of sp-specific dispersal kernels for the islands: abundance matrices with parameters as key:
	map<vector<double>, Matrix<float> > speciesDispMain;
	
	//map of sp-specific suitable habitat for colonization from mainland: abundance matrices with parameters as key:
	map<vector<double>, Matrix<float> > speciesHabtMain;

	//save the prior map in order to be used for each replicate with the same Mainland species pool:
	map<vector<double>, Matrix<float> > speciesInitialHabtMain;


	//maps for the mutants
	map<vector<double>, Matrix<float> > speciesAdultMutant;
	map<vector<double>, Matrix<float> > speciesRecruitMutant;
	map<vector<double>, Matrix<float> > speciesSeedMutant;
	map<vector<double>, Matrix<float> > speciesDispMutant;
    map<vector<double>, Matrix<float> > speciesHabtMutant;

	//Output:
	std::string simulationName;
	unsigned int Col, DemExt, VulcanicExt,RandomExt, ErosionExt, SpecClado, SpecAna;
	//to exclude pre-species from this:
	unsigned int prespeciesExtinct;

	//Disturance output
	unsigned int disturbRadius, disturbOccurrence, epicenterX,epicenterY;

	//tracking single species in single species simulation:
	unsigned int currentSpecies;

};