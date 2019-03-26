#include "Simulation.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <map>
#include <sstream>
#include <stdlib.h>
#include <cstdlib>
#include <time.h>
#include <algorithm>
#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <armadillo>
#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif

using std::ifstream;
using std::ofstream;
using std::map;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::pair;
using namespace arma;

Simulation::Simulation(const ConfigFile& config)
{
	string dispFileName = config.read<string>("dispersalFile");
	string mapFileName = config.read<string>("mapFile");
	string ParameterRangesTree = config.read<string>("parameterRangesTree");
	string ParameterRangesShrub = config.read<string>("parameterRangesShrub");
	string ParameterRangesHerb = config.read<string>("parameterRangesHerb");
	string simuSpecifications=config.read<string>("SimuSpecifications"); 
	string SimulatedPFT=config.read<string>("SimulatedPFT"); 
	unsigned int SimulationTimesteps  = config.read<unsigned int>("SimulationTimesteps");
	unsigned int Replicates  = config.read<unsigned int>("Replicates");
	unsigned int IslandTimestep  = config.read<unsigned int>("IslandTimestep");
	unsigned int MigrationEvents  = config.read<unsigned int>("MigrationEvents");
	unsigned int MaxIslandRadius  = config.read<unsigned int>("MaxIslandRadius");
	unsigned int MainlandSpeciesNumber  = config.read<unsigned int>("MainlandSpeciesNumber"); 
	unsigned int InvasiveSpeciesNumber  = config.read<unsigned int>("InvasiveSpeciesNumber");
	unsigned int SpeciesOutputtimestep  = config.read<unsigned int>("SpeciesOutputtimestep");  
	cout<< "SimulationTimesteps" << SimulationTimesteps << endl;
	string DispersalKernel = config.read<string>("DispersalKernel");  
	string HabSuitability = config.read<string>("HabitatSuitability");  
	string ReadMainlandSpeciespool = config.read<string>("ReadMainlandSpeciespool"); 
	string TreeTemperatureFeedback = config.read<string>("TreeTemperatureFeedback"); 
	string Colonize = config.read<string>("Colonize"); 
	string SingleSpecies = config.read<string>("SingleSpecies"); 
	string MetabolicRateToCorrect = config.read<string>("MetabolicRateToCorrect");   
	string Treeline = config.read<string>("Treeline"); 
	string CladoSpeciation = config.read<string>("CladoSpeciation");
	string AnaSpeciation = config.read<string>("AnaSpeciation"); 
	string AltitudinalTemp = config.read<string>("AltitudinalTemp");  
	unsigned int KernelRadius = config.read<unsigned int>("KernelRadius");
	float PFTtree = config.read<float>("PFTtree");
	float PFTshrub = config.read<float>("PFTshrub");
	float PFTherb = config.read<float>("PFTherb"); 
	float DisturbanceProb = config.read<float>("DisturbanceProb");
	float Allometabolic = config.read<float>("Allometabolic");
	float Temperature = config.read<float>("Temperature"); 
	float TempForFixedMetabolism = config.read<float>("TempForFixedMetabolism");   
	string Boundaries = config.read<string>("Boundaries");
	string VulnerableSpecialists = config.read<string>("VulnerableSpecialists"); 
	string FixedIslandaltitude = config.read<string>("FixedIslandaltitude");
	string DynamicGeology = config.read<string>("DynamicGeology");
	float boltzmann = config.read<float>("Boltzmann");
	float e = config.read<float>("E");
	float MaxK = config.read<float>("MaxK");
	float LDD = config.read<float>("LDD");  
	unsigned int RandomDisturbance = config.read<unsigned int>("RandomDisturbance");
	unsigned int VulcanicDisturbance = config.read<unsigned int>("VulcanicDisturbance");  
	unsigned int InvasivesperTimestep = config.read<unsigned int>("InvasivesperTimestep");
	unsigned int InvasivesBegin = config.read<unsigned int>("InvasivesBegin");
	unsigned int RandomInvasivecenters = config.read<unsigned int>("RandomInvasivecenters");
	unsigned int ClimateChange = config.read<unsigned int>("ClimateChange");
	unsigned int TemperatureChangeInterval = config.read<unsigned int>("TemperatureChangeInterval"); 
	unsigned int Verticalpreference = config.read<unsigned int>("Verticalpreference"); 
	unsigned int MaxOptaltitude = config.read<unsigned int>("MaxOptaltitude"); 
	int MinOptaltitude = config.read<unsigned int>("MinOptaltitude");
	unsigned int LowlandDegradation = config.read<unsigned int>("LowlandDegradation"); 
	unsigned int DegradationImpact = config.read<unsigned int>("DegradationImpact"); 
	unsigned int SpeciesAbundDistriBegin = config.read<unsigned int>("SpeciesAbundDistriBegin");  
	unsigned int MountainRadius = config.read<unsigned int>("MountainRadius");  

	float TemperatureChange = config.read<float>("TemperatureChange"); 
	float NicheEvolution = config.read<float>("NicheEvolution");

	kernelRadius = KernelRadius;
	replicates=Replicates;
	simulationTimesteps= SimulationTimesteps;
	islandTimestep= IslandTimestep;
	migrationEvents= MigrationEvents;
	maxIslandRadius = MaxIslandRadius;
	speciesOutputtimestep=SpeciesOutputtimestep;
	parameterRangesTree = ParameterRangesTree;
	parameterRangesShrub = ParameterRangesShrub;
	parameterRangesHerb = ParameterRangesHerb;
	mainlandSpeciesNumber = MainlandSpeciesNumber;
	invasiveSpeciesNumber= InvasiveSpeciesNumber;
	dispKernel =DispersalKernel;
	simulatedPFT = SimulatedPFT;
	habSuitability=HabSuitability;
	pftTree= PFTtree;
	pftShrub = PFTshrub;
	pftHerb = PFTherb;
	randomDisturbance=RandomDisturbance;
	disturbanceProb=DisturbanceProb;
	vulcanicDisturbance=VulcanicDisturbance;
	allometabolic = Allometabolic;
	temperature = Temperature;
	tempForFixedMetabolism = TempForFixedMetabolism;
	boundaries= Boundaries;
	vulnerableSpecialists=VulnerableSpecialists;
	fixedIslandaltitude =FixedIslandaltitude;
	dynamicGeology = DynamicGeology;
	Boltzmann = boltzmann;
	E = e;
	maxK = MaxK;
	readMainlandSpeciespool=ReadMainlandSpeciespool;
	treeTemperatureFeedback =TreeTemperatureFeedback;
	colonize = Colonize;
	singleSpecies =SingleSpecies ;
	metabolicRateToCorrect=MetabolicRateToCorrect;
	invasivesperTimestep=InvasivesperTimestep;
	invasivesBegin = InvasivesBegin;
	randomInvasivecenters = RandomInvasivecenters;
	cladospeciation = CladoSpeciation;
	anaspeciation = AnaSpeciation;
	treeline = Treeline;
	temperatureChange =TemperatureChange;
	climateChange  = ClimateChange;
	temperatureChangeInterval =TemperatureChangeInterval;
	verticalpreference =Verticalpreference;
	maxOptaltitude = MaxOptaltitude;
	minOptaltitude = MinOptaltitude;
	altitudinalTemp = AltitudinalTemp;
	lowlandDegradation =LowlandDegradation ;
	degradationImpact =DegradationImpact ;
	speciesAbundDistriBegin = SpeciesAbundDistriBegin;
	ldd= LDD;
	nicheEvolution = NicheEvolution;
	mountainRadius =MountainRadius;
	//std::string spdisp;		
    //readDispersalMatrix(dispFileName);      // Dispersal matrix file
	readMap(mapFileName);  // landscape matrix file
	simulationName = simuSpecifications;
}
Simulation::~Simulation()
{
	/*if(model != 0) {
		delete model;
	}*/
}

void Simulation::readDispersalMatrix(const std::string& dispFileName) {
	ifstream dispFileStream(dispFileName.c_str());
	string nix;
	dispFileStream >> nix >> DX >> nix >> DY;
	Disp = Matrix<double> (DX,DY);
	for (unsigned int r=0;r<DY;++r)
	{
		for (unsigned int c=0;c<DX;++c) 
		{
			dispFileStream>>Disp(c, r);
		}
	}
	dispFileStream.close();	
	cout << "Disp:" << endl << Disp << endl;
}

void Simulation::readMap(const std::string& mapFileName) {
	ifstream mapFileStream(mapFileName.c_str());

	string nix;
	mapFileStream >> nix >> NX >> nix >> NY
		>> nix >> nix >> nix >> nix >> nix >> nix >> nix >> nix;

	gisdata2=Matrix<int> (NX,NY);
    
	patchesPosition.clear();

	patchesPositionMainLand.clear();

	for (unsigned int r=0;r<NY;++r)
	{
		for (unsigned int c=0;c<NX;++c)
		{
            mapFileStream >> gisdata2(c,r);

			if(gisdata2(c,r)>=1) 
			{
				//if the cell content is 1, add the cell to the patches
				patchesPosition.push_back(std::pair<int,int>(c,r));
			}
			if(gisdata2(c,r)==0) 
			{
				//if the cell content is 1, add the cell to the patches
				patchesPositionMainLand.push_back(std::pair<int,int>(c,r));
			}
		}
	}
	
	mapFileStream.close();
	cout << "islandsize:" << patchesPosition.size() << endl;
	cout << "mainlandsize:" << patchesPositionMainLand.size() << endl;
}

void Simulation::perform()
{
    replicate=1;
	while (replicate<=replicates)
	{ 
		//if (interactionID== interactionsPerFile*(outputfilenumber-1))
		
		if(singleSpecies=="true"){
			
			currentSpecies=0;
			for(unsigned int sp=0; sp<mainlandSpeciesNumber; sp++){
				std::stringstream DisturbanceOutputFileName;
				DisturbanceOutputFileName<< "Disturbance." + simulationName + ".Species." << sp+1 << "." << replicate << ".txt";

				ofstream exportd7(DisturbanceOutputFileName.str().c_str());
				exportd7<<"Replicate\tTimestep\tIslandSize\tOccurrence\tSize\tEpicenterX\tEpicenterY\n";   
							
				std::stringstream AlphaSpeciesOutputFileName;
				AlphaSpeciesOutputFileName <<  "AlphaSpecies." + simulationName + ".Species." << sp+1 << "." << replicate <<".txt";

				ofstream exportd9(AlphaSpeciesOutputFileName.str().c_str());
				exportd9<<"Replicate\tTimestep\tCellX\tCellY\tFenology\tSpeciesType\tLineage\tMotherSpID\tMass\tAnnual\tAbundA\tAbundR\tAbundS\n";  
		
				std::stringstream SpeciesOutputFileName;
				SpeciesOutputFileName << "Species." + simulationName + ".Species." << sp+1 << "." << replicate << ".txt";
		
				ofstream exportd3(SpeciesOutputFileName.str().c_str());
				exportd3<<"Replicate\tTimestep\tFenology\tSpeciesType\tAnnual\tAltitude\tSuitability\tGeneralism\tMass\tC\tAlpha\tTailFatness\tMseed\tCellNumber\tMeanAbundA\tMeanAbundR\tMeanAbundS\tverticalgeneralist\tIslandSide\tSpID\tMotherspID\tLineage\n";  
				//exportd3<<"Replicate\tTimestep\tFenology\tEndemic\tInvasive\tAnnual\tAltitude\tSuitability\tGeneralism\tMass\tC\tAlpha\tTailFatness\tMseed\tCellNumber\tMeanAbundA\tMeanAbundR\tMeanAbundS\tverticalgeneralist\tIslandSide\tSpID\tMotherspID\tLineage\n";  
			
				std::stringstream LineagesOutputFileName;
				LineagesOutputFileName <<  "Lineages." + simulationName + ".Species." << sp+1 << "." << replicate << ".txt";
		  
				ofstream exportd4(LineagesOutputFileName.str().c_str());
				exportd4<<"Replicate\tTimestep\tLineage\tFenology\tSpeciesType\tAnnual\tGeneralism\tMass\tC\tAlpha\tTailFatness\tMseed\tCellNumber\tIslandSide\tAltitude\tSuitability\tMutantID\tMotherspID\tverticalgeneralist\n";  
				
				std::stringstream NewSpeciesOutputFileName;
				NewSpeciesOutputFileName << "NewSpecies." + simulationName + ".Species." << sp+1 << "."  << replicate << ".txt";
		
				ofstream exportd11(NewSpeciesOutputFileName.str().c_str());
				exportd11<<"Replicate\tTimestep\tFenology\tAnnual\tAltitude\tSuitability\tGeneralism\tMass\tC\tAlpha\tTailFatness\tMseed\tverticalgeneralist\tIslandSide\tSpID\tMotherspID\tLineage\n";  
			

				initialization();
		  		int esportd3Step=0; //can be a configuration variable
				clock_t start = clock();

				while(timeStep<simulationTimesteps)
				{   				
					//Simulation:
					if (colonize=="island") colonization();
					if ((colonize=="forest")&&(timeStep==0)) colonizationForest();
					if((invasivesBegin>0)&&(timeStep>=invasivesBegin)) invasiveColonization();
					if(allometabolic==1) localdynamics(speciesAdult,speciesRecruit,speciesSeed, Ktree,speciesDisp, speciesHabt,speciesHabtMain);
					else localdynamicsNoMTE(speciesAdult,speciesRecruit,speciesSeed, Ktree,speciesDisp, speciesHabt,speciesHabtMain);
					int exportd3Tobewritten= 0;
					if(timeStep==speciesOutputtimestep*esportd3Step) 
					{
						esportd3Step++;
						exportd3Tobewritten= 1;
					}

					timeStep++;
					islanddynamics(speciesAdult,speciesRecruit,speciesSeed, Ktree,speciesDisp, speciesHabt,speciesHabtMain,TempKelvinBareIsland);

					if(exportd3Tobewritten==1)
					{//Output abundance and alpha richness:
						for(unsigned int i=0; i<patchesPosition.size(); i++) 
						{
							//cout << "speciesAdult.size()  " << speciesAdult.size()<< endl;
							for(map<vector<double>, Matrix<float> >::iterator itr = speciesAdult.begin(); itr != speciesAdult.end(); ++itr) {
								const vector<double>& parameters = itr->first; // key (== parameterSet)
								Matrix<float>& abundances = itr->second; // value (== abundances)
								if (timeStep>=speciesAbundDistriBegin){ //lowlandDegradation-speciesOutputtimestep*2
									exportd9<<replicate<<"\t"<<timeStep<<"\t"
									<<patchesPosition[i].first<<"\t"<<patchesPosition[i].second<<"\t"<<parameters[0]<<"\t"<<parameters[25]<<"\t"<<parameters[17]<<"\t"<<parameters[18]<<
									"\t"<<parameters[16]<<"\t"<<parameters[15]<<
									"\t"<<abundances(patchesPosition[i])<<"\t"<<speciesRecruit[parameters](patchesPosition[i])<<"\t"<<speciesSeed[parameters](patchesPosition[i])<<"\n";	
								}

							}

						}

					}
					for(map<vector<double>, Matrix<float> >::iterator itr = speciesAdult.begin(); itr != speciesAdult.end(); ++itr) 
					{	
						const vector<double>& parameters = itr->first; // key (== parameterSet)
						Matrix<float>& abundances = itr->second; // value (== abundances) 
						unsigned int cellnumber =0;
				
						float sumabundA =0;		
						float sumabundR =0;
						float sumabundS =0;
						if((parameters[25]==2)||(parameters[25]==1))
						{
						   exportd11<<replicate<<"\t"<<timeStep<<"\t"
							<<parameters[0]<<"\t"<<parameters[15]<<"\t"<<parameters[19]<<"\t"<<parameters[9]<<
							"\t"<<parameters[23]<<"\t"<<parameters[16]<<"\t"<<parameters[3]<<
							"\t"<<parameters[6]<<"\t"<<parameters[7]<<"\t"<<parameters[14]<<"\t"
							<<parameters[24]<<"\t"<<parameters[10]<<
							"\t"<<parameters[17]*1000+parameters[16]<<"\t"<<parameters[18]<<"\t"<<parameters[17]<<"\n";	
						}
						if(exportd3Tobewritten==1)
						{
							for(unsigned int i=0; i<patchesPosition.size(); i++)
							{
								bool counted = false;
								if(abundances(patchesPosition[i])>0) {
									if(!counted) {
										cellnumber++;	
										counted=true;
									}
									sumabundA += abundances(patchesPosition[i]);
								}

								if(speciesRecruit[parameters](patchesPosition[i])>0) {
									if(!counted) {
										cellnumber++;	
										counted=true;
									}	
									sumabundR += speciesRecruit[parameters](patchesPosition[i]);
								}
								if(speciesSeed[parameters](patchesPosition[i])>0) {
									if(!counted) {
										cellnumber++;	
										counted=true;
									}	
									sumabundS += speciesSeed[parameters](patchesPosition[i]);
								}
							}
						}	
						
						// radiated endemics:
						if((parameters[25]==2)||(parameters[25]==6)) 
						{
							if(exportd3Tobewritten==1)
							{
								exportd4<<replicate<<"\t"<<timeStep<<"\t"<<parameters[17]<<
								"\t"<<parameters[0]<<"\t"<<parameters[25]<<"\t"<<parameters[15]<<"\t"<<parameters[23]<<"\t"<<parameters[16]<<
								"\t"<<parameters[3]<<
								"\t"<<parameters[6]<<"\t"<<parameters[7]<<"\t"<<parameters[14]<<"\t"
								<<cellnumber<<"\t"<<parameters[10]<<"\t"<<parameters[19]<<"\t"<<parameters[9]<<
								"\t"<<parameters[17]*1000+parameters[16]*0.001<<"\t"<<parameters[18]<<"\t"<<parameters[24]<<"\n";	
							}
					
						}		
						//differentiated endemics:
						if(parameters[25]==3) 
						{
							if(exportd3Tobewritten==1)
							{
								exportd4<<replicate<<"\t"<<timeStep<<"\t"<<parameters[17]<<
								"\t"<<parameters[0]<<"\t"<<parameters[25]<<"\t"<<parameters[15]<<"\t"<<parameters[23]<<"\t"<<parameters[16]<<
								"\t"<<parameters[3]<<
								"\t"<<parameters[6]<<"\t"<<parameters[7]<<"\t"<<parameters[14]<<"\t"
								<<cellnumber<<"\t"<<parameters[10]<<"\t"<<parameters[19]<<"\t"<<parameters[9]<<
								"\t"<<parameters[17]*1000+parameters[16]*0.001<<"\t"<<parameters[18]<<"\t"<<parameters[24]<<"\n";	
							}
					
						}
						if(exportd3Tobewritten==1)
						{
							if (sumabundA!=0) sumabundA=sumabundA/cellnumber;
							if (sumabundR!=0) sumabundR=sumabundR/cellnumber;
							if (sumabundS!=0) sumabundS=sumabundS/cellnumber;

							exportd3<<replicate<<"\t"<<timeStep<<"\t"
							<<parameters[0]<<"\t"<<parameters[25]<<"\t"<<parameters[15]<<"\t"<<parameters[19]<<"\t"<<parameters[9]<<
							"\t"<<parameters[23]<<"\t"<<parameters[16]<<"\t"<<parameters[3]<<
							"\t"<<parameters[6]<<"\t"<<parameters[7]<<"\t"<<parameters[14]<<"\t"
							<<cellnumber<<"\t"<<sumabundA<<"\t"<<sumabundR<<"\t"<<sumabundS<<"\t"<<parameters[24]<<"\t"<<parameters[10]<<
							"\t"<<parameters[17]*1000+parameters[16]*0.001<<"\t"<<parameters[18]<<"\t"<<parameters[17]<<"\n";	
						}
					}
					
					if((vulcanicDisturbance==1)||(randomDisturbance==1)){
						exportd7<<replicate<<"\t"<<timeStep<<"\t"<<patchesPosition.size()<<"\t"<<disturbOccurrence<<"\t"<<				
							disturbRadius<<"\t"<<epicenterX<<"\t"<<epicenterY<<"\n";
						disturbRadius=0;
						disturbOccurrence=0;
						epicenterX=0;
						epicenterY=0;
					}
            
				}
		
				clock_t duration = clock() - start;
				cout << "Replicate " << replicate << " : " << duration << " ms" << endl;
				currentSpecies++;	
				exportd3.close();
				exportd4.close();
				exportd7.close();
				exportd9.close();
			}
			replicate++;
		}
		if(singleSpecies=="false") {
			std::stringstream GammaRichnessTOutputFileName;
			GammaRichnessTOutputFileName<< "GammaRichness.Total." + simulationName + "." << replicate << ".txt";

			ofstream exportd(GammaRichnessTOutputFileName.str().c_str());
			exportd<<"Replicate\tTimestep\tIslandSize\tGammaTreeRichness\tGammaShrubRich\tGammaHerbRich\tGammaTreeRichA\tGammaShrubRichA\tGammaHerbRichA\tCol\tDemExt\tRandomExt\tVulcanicExt\tErosionExt\tSpecClado\tSpecAna\tEndemicCladoTrees\tEndemicCladoShrubs\tEndemicCladoHerbs\tEndemicCladoTreesA\tEndemicCladoShrubsA\tEndemicCladoHerbsA\tEndemicAnaTrees\tEndemicAnaShrubs\tEndemicAnaHerbs\tEndemicAnaTreesA\tEndemicAnaShrubsA\tEndemicAnaHerbsA\n";   
		
			std::stringstream GammaRichnessInvasivesTOutputFileName;
			GammaRichnessInvasivesTOutputFileName<< "GammaRichness.Total.Invasives." + simulationName + "." << replicate << ".txt";

			ofstream exportd10(GammaRichnessInvasivesTOutputFileName.str().c_str());
			exportd10<<"Replicate\tTimestep\tIslandSize\tCol\tDemExt\tRandomExt\tVulcanicExt\tErosionExt\tSpecClado\tSpecAna\tInvasiveTrees\tInvasiveShrubs\tInvasiveHerbs\tInvasiveTreesA\tInvasiveShrubsA\tInvasiveHerbsA\n";   
		
			std::stringstream GammaRichnessRSOutputFileName;
			GammaRichnessRSOutputFileName<< "GammaRichness.NonAdults." + simulationName + "." << replicate << ".txt";

			ofstream exportd6(GammaRichnessRSOutputFileName.str().c_str());
			exportd6<<"Replicate\tTimestep\tIslandSize\tGammaTreeRichR\tGammaShrubRichR\tGammaHerbRichR\tGammaTreeRichS\tGammaShrubRichS\tGammaHerbRichS\tEndemicCladoTreesR\tEndemicCladoShrubsR\tEndemicCladoHerbsR\tEndemicCladoTreesS\tEndemicCladoShrubsS\tEndemicCladoHerbsS\tEndemicAnaTreesR\tEndemicAnaShrubsR\tEndemicAnaHerbsR\tEndemicAnaTreesS\tEndemicAnaShrubsS\tEndemicAnaHerbsS\n";   
	    
			std::stringstream DisturbanceOutputFileName;
			DisturbanceOutputFileName<< "Disturbance." + simulationName + "." << replicate << ".txt";

			ofstream exportd7(DisturbanceOutputFileName.str().c_str());
			exportd7<<"Replicate\tTimestep\tIslandSize\tOccurrence\tSize\tEpicenterX\tEpicenterY\n";   
							
			std::stringstream AlphaRichnessOutputFileName;
			AlphaRichnessOutputFileName <<  "AlphaRichness." + simulationName + "." << replicate <<".txt";

			ofstream exportd2(AlphaRichnessOutputFileName.str().c_str());
			exportd2<<"Replicate\tTimestep\tCellID\tAlphaTreeRichness\tAlphaShrubRichness\tAlphaHerbRichness\tAlphaETreeRichness\tAlphaEShrubRichness\tAlphaEHerbRichness\tK\tTotAbundTreeA\tTotAbundShrubA\tTotAbundHerbA\tTotAbundTreeR\tTotAbundShrubR\tTotAbundHerbR\tTotAbundTreeS\tTotAbundShrubS\tTotAbundHerbS\n";  
	
			std::stringstream AlphaRichnessInvasiveOutputFileName;
			AlphaRichnessInvasiveOutputFileName <<  "AlphaRichness.Invasive." + simulationName + "." << replicate <<".txt";

			ofstream exportd8(AlphaRichnessInvasiveOutputFileName.str().c_str());
			exportd8<<"Replicate\tTimestep\tCellID\tAlphaTreeRichness\tAlphaShrubRichness\tAlphaHerbRichness\tTotAbundTreeA\tTotAbundShrubA\tTotAbundHerbA\tTotAbundTreeR\tTotAbundShrubR\tTotAbundHerbR\tTotAbundTreeS\tTotAbundShrubS\tTotAbundHerbS\n";  
		
			std::stringstream AlphaSpeciesOutputFileName;
			AlphaSpeciesOutputFileName <<  "AlphaSpecies." + simulationName + "." << replicate <<".txt";

			ofstream exportd9(AlphaSpeciesOutputFileName.str().c_str());
			exportd9<<"Replicate\tTimestep\tCellX\tCellY\tFenology\tSpeciesType\tLineage\tMotherSpID\tMass\tAnnual\tAbundA\tAbundR\tAbundS\n";  
		
			std::stringstream SpeciesOutputFileName;
			SpeciesOutputFileName << "Species." + simulationName + "." << replicate << ".txt";
		
			ofstream exportd3(SpeciesOutputFileName.str().c_str());
			exportd3<<"Replicate\tTimestep\tFenology\tSpeciesType\tAnnual\tAltitude\tSuitability\tGeneralism\tMass\tC\tAlpha\tTailFatness\tMseed\tCellNumber\tMeanAbundA\tMeanAbundR\tMeanAbundS\tverticalgeneralist\tIslandSide\tSpID\tMotherspID\tLineage\n";  
			//exportd3<<"Replicate\tTimestep\tFenology\tEndemic\tInvasive\tAnnual\tAltitude\tSuitability\tGeneralism\tMass\tC\tAlpha\tTailFatness\tMseed\tCellNumber\tMeanAbundA\tMeanAbundR\tMeanAbundS\tverticalgeneralist\tIslandSide\tSpID\tMotherspID\tLineage\n";  
			
			std::stringstream NewSpeciesOutputFileName;
			NewSpeciesOutputFileName << "NewSpecies." + simulationName + "." << replicate << ".txt";
		
			ofstream exportd11(NewSpeciesOutputFileName.str().c_str());
			exportd11<<"Replicate\tTimestep\tFenology\tAnnual\tAltitude\tSuitability\tGeneralism\tMass\tC\tAlpha\tTailFatness\tMseed\tverticalgeneralist\tIslandSide\tSpID\tMotherspID\tLineage\n";  
			

			std::stringstream LineagesOutputFileName;
			LineagesOutputFileName <<  "Lineages." + simulationName + "." << replicate << ".txt";
		  
			ofstream exportd4(LineagesOutputFileName.str().c_str());
			exportd4<<"Replicate\tTimestep\tLineage\tFenology\tSpeciesType\tAnnual\tGeneralism\tMass\tC\tAlpha\tTailFatness\tMseed\tCellNumber\tIslandSide\tAltitude\tSuitability\tMutantID\tMotherspID\tverticalgeneralist\n";  
			//exportd4<<"Replicate\tTimestep\tLineage\tFenology\tInvasive\tAnnual\tGeneralism\tMass\tC\tAlpha\tTailFatness\tMseed\tCellNumber\tIslandSide\tAltitude\tSuitability\tMutantID\tMotherspID\tverticalgeneralist\n";  
		
			std::stringstream MainlandSpeciesPoolOutputFileName;
			MainlandSpeciesPoolOutputFileName <<  "MainlandSpeciesPool." + simulationName + ".txt";

			ofstream exportd5(MainlandSpeciesPoolOutputFileName.str().c_str());
			exportd5<<"Fenology\tRmax\tM\tGamma\tGenTime\tMutRate\talpha\tP\tKpftfactor\trandomTrend\tislandSide\tGrowth\tMrecruit\tGerm\tMseed\tAnnual\tMass\tLineage\tMotherspID\toptimalAltitude\tMassSeedling\tMassSeed\tKpftfactorSeedling\tgeneralist\tverticalgeneralist\tspeciesType\tyearOfSpec\n";  

		  	initialization();
			int esportd3Step=0; //can be a configuration variable
			clock_t start = clock();

			for (unsigned int i=0; i< mainlandSpeciesNumber; i++)
					{
					   exportd5<<mainlandSpeciesPool[i][0]<<"\t"<<mainlandSpeciesPool[i][1]<<"\t"<<mainlandSpeciesPool[i][2]<<"\t"<<mainlandSpeciesPool[i][3]<<"\t"<<
								mainlandSpeciesPool[i][4]<<"\t"<<mainlandSpeciesPool[i][5]<<"\t"<<mainlandSpeciesPool[i][6]<<"\t"<<mainlandSpeciesPool[i][7]<<"\t"<<
								mainlandSpeciesPool[i][8]<<"\t"<<mainlandSpeciesPool[i][9]<<"\t"<<mainlandSpeciesPool[i][10]<<"\t"<<mainlandSpeciesPool[i][11]<<"\t"<<
								mainlandSpeciesPool[i][12]<<"\t"<<mainlandSpeciesPool[i][13]<<"\t"<<mainlandSpeciesPool[i][14]<<"\t"<<mainlandSpeciesPool[i][15]<<"\t"<<
								mainlandSpeciesPool[i][16]<<"\t"<<mainlandSpeciesPool[i][17]<<"\t"<<mainlandSpeciesPool[i][18]<<"\t"<<mainlandSpeciesPool[i][19]<<"\t"<<
								mainlandSpeciesPool[i][20]<<"\t"<<mainlandSpeciesPool[i][21]<<"\t"<<mainlandSpeciesPool[i][22]<<"\t"<<mainlandSpeciesPool[i][23]<<"\t"<<
								mainlandSpeciesPool[i][24]<<"\t"<<mainlandSpeciesPool[i][25]<<"\t"<<mainlandSpeciesPool[i][26]<<"\n";
					}	
			exportd5.close();
			while(timeStep<simulationTimesteps)
			{   				
				//Simulation:
				if (colonize=="island") colonization();
				if ((colonize=="forest")&&(timeStep==0)) colonizationForest();
				if((invasivesBegin>0)&&(timeStep>=invasivesBegin)) invasiveColonization();
				if(allometabolic==1) localdynamics(speciesAdult,speciesRecruit,speciesSeed, Ktree,speciesDisp, speciesHabt,speciesHabtMain);
				else localdynamicsNoMTE(speciesAdult,speciesRecruit,speciesSeed, Ktree,speciesDisp, speciesHabt,speciesHabtMain);
				int exportd3Tobewritten= 0;
				if(timeStep==speciesOutputtimestep*esportd3Step) 
				{
					esportd3Step++;
					exportd3Tobewritten= 1;
				}

				timeStep++;				
				islanddynamics(speciesAdult,speciesRecruit,speciesSeed, Ktree,speciesDisp, speciesHabt,speciesHabtMain,TempKelvinBareIsland);

				if(exportd3Tobewritten==1)
				{//Output abundance and alpha richness:
					for(unsigned int i=0; i<patchesPosition.size(); i++) 
					{
						unsigned int alphaTreeRichness=0;
						unsigned int alphaShrubRichness=0;
						unsigned int alphaHerbRichness=0;

						unsigned int alphaTreeRichnessPreSpecies=0;
						unsigned int alphaShrubRichnessPreSpecies=0;
						unsigned int alphaHerbRichnessPreSpecies=0;

						float cellAbundTreeA=0;
						float cellAbundShrubA=0;
						float cellAbundHerbA=0;
						float cellAbundTreeR=0;
						float cellAbundShrubR=0;
						float cellAbundHerbR=0;
						float cellAbundTreeS=0;
						float cellAbundShrubS=0;
						float cellAbundHerbS=0;

						unsigned int alphaTreeRichnessEndemic=0;
						unsigned int alphaShrubRichnessEndemic=0;
						unsigned int alphaHerbRichnessEndemic=0;

						unsigned int alphaTreeRichnessI=0;
						unsigned int alphaShrubRichnessI=0;
						unsigned int alphaHerbRichnessI=0;
						float cellAbundTreeAI=0;
						float cellAbundShrubAI=0;
						float cellAbundHerbAI=0;
						float cellAbundTreeRI=0;
						float cellAbundShrubRI=0;
						float cellAbundHerbRI=0;
						float cellAbundTreeSI=0;
						float cellAbundShrubSI=0;
						float cellAbundHerbSI=0;
						
						for(map<vector<double>, Matrix<float> >::iterator itr = speciesAdult.begin(); itr != speciesAdult.end(); ++itr) {
							const vector<double>& parameters = itr->first; // key (== parameterSet)
							Matrix<float>& abundances = itr->second; // value (== abundances)
							if (timeStep>=speciesAbundDistriBegin){ //lowlandDegradation-speciesOutputtimestep*2
								exportd9<<replicate<<"\t"<<timeStep<<"\t"
								<<patchesPosition[i].first<<"\t"<<patchesPosition[i].second<<"\t"<<parameters[0]<<"\t"<<parameters[25]<<"\t"<<parameters[17]<<"\t"<<parameters[18]<<
								"\t"<<parameters[16]<<"\t"<<parameters[15]<<
								"\t"<<abundances(patchesPosition[i])<<"\t"<<speciesRecruit[parameters](patchesPosition[i])<<"\t"<<speciesSeed[parameters](patchesPosition[i])<<"\n";	
							}
							if (abundances(patchesPosition[i])>0) 
							{
						
								if (parameters[16]>100000) {
									cellAbundTreeA+=abundances(patchesPosition[i]);
									if ((parameters[25]==0)||(parameters[25]==2)||(parameters[25]==3)||(parameters[25]==4)||(parameters[25]==6)) alphaTreeRichness++;
									if ((parameters[25]==1)||(parameters[25]==5)) alphaTreeRichnessPreSpecies++;
									if((parameters[25]==2)||(parameters[25]==3)||(parameters[25]==6)) alphaTreeRichnessEndemic++;
								}
								if ((parameters[16]>1000) && (parameters[16]<=100000)) {
									cellAbundShrubA+=abundances(patchesPosition[i]);
									if ((parameters[25]==0)||(parameters[25]==2)||(parameters[25]==3)||(parameters[25]==4)||(parameters[25]==6)) alphaShrubRichness++;
									if ((parameters[25]==1)||(parameters[25]==5)) alphaShrubRichnessPreSpecies++;
									if((parameters[25]==2)||(parameters[25]==3)||(parameters[25]==6)) alphaShrubRichnessEndemic++;
								}
								if (parameters[16]<1000){
									cellAbundHerbA+=abundances(patchesPosition[i]);
									if(parameters[15]==0){
										if ((parameters[25]==0)||(parameters[25]==2)||(parameters[25]==3)||(parameters[25]==4)||(parameters[25]==6))alphaHerbRichness++;
										if ((parameters[25]==1)||(parameters[25]==5)) alphaHerbRichnessPreSpecies++;
										if((parameters[25]==2)||(parameters[25]==3)||(parameters[25]==6)) alphaShrubRichnessEndemic++;
									}
								}
							}
							if (speciesRecruit[parameters](patchesPosition[i])>0) 
							{
						
								if (parameters[16]>100000) {
									cellAbundTreeR+=speciesRecruit[parameters](patchesPosition[i]);
								}
								if ((parameters[16]>1000) && (parameters[16]<=100000)) {
									cellAbundShrubR+=speciesRecruit[parameters](patchesPosition[i]);
								}
								if (parameters[16]<1000){
									cellAbundHerbR+=speciesRecruit[parameters](patchesPosition[i]);
								}
							}
							if (speciesSeed[parameters](patchesPosition[i])>0) 
							{
						
								if (parameters[16]>100000) {
									cellAbundTreeS+=speciesSeed[parameters](patchesPosition[i]);
								}
								if ((parameters[16]>1000) && (parameters[16]<=100000)) {
									cellAbundShrubS+=speciesSeed[parameters](patchesPosition[i]);
								}
								if (parameters[16]<1000){
									cellAbundHerbS+=speciesSeed[parameters](patchesPosition[i]);
									if(parameters[15]==1){
										if ((parameters[25]==0)||(parameters[25]==2)||(parameters[25]==3)||(parameters[25]==4)||(parameters[25]==6))alphaHerbRichness++;
										if ((parameters[25]==1)||(parameters[25]==5)) alphaHerbRichnessPreSpecies++;
										if((parameters[25]==2)||(parameters[25]==3)||(parameters[25]==6)) alphaShrubRichnessEndemic++;
									}
								}
							}
							//Invasive
							if (parameters[25]==4){
								if (abundances(patchesPosition[i])>0) 
								{
						
									if (parameters[16]>100000) {
										cellAbundTreeAI+=abundances(patchesPosition[i]);
										alphaTreeRichnessI++;
									}
									if ((parameters[16]>1000) && (parameters[16]<=100000)) {
										cellAbundShrubAI+=abundances(patchesPosition[i]);
										alphaShrubRichnessI++;
									}
									if (parameters[16]<1000){
										cellAbundHerbAI+=abundances(patchesPosition[i]);
										alphaHerbRichnessI++;
									}
								}
								if (speciesRecruit[parameters](patchesPosition[i])>0) 
								{
						
									if (parameters[16]>100000) {
										cellAbundTreeRI+=speciesRecruit[parameters](patchesPosition[i]);
									}
									if ((parameters[16]>1000) && (parameters[16]<=100000)) {
										cellAbundShrubRI+=speciesRecruit[parameters](patchesPosition[i]);
									}
									if (parameters[16]<1000){
										cellAbundHerbRI+=speciesRecruit[parameters](patchesPosition[i]);
									}
								}
								if (speciesSeed[parameters](patchesPosition[i])>0) 
								{
						
									if (parameters[16]>100000) {
										cellAbundTreeSI+=speciesSeed[parameters](patchesPosition[i]);
									}
									if ((parameters[16]>1000) && (parameters[16]<=100000)) {
										cellAbundShrubSI+=speciesSeed[parameters](patchesPosition[i]);
									}
									if (parameters[16]<1000){
										cellAbundHerbSI+=speciesSeed[parameters](patchesPosition[i]);
									}
								}
							}
						}

						exportd2<<replicate<<"\t"<<timeStep<<"\t"<<(patchesPosition[i].first*1000) + patchesPosition[i].second
								<<"\t"<<alphaTreeRichness<<"\t"<<alphaShrubRichness<<"\t"<<alphaHerbRichness
								<<"\t"<<alphaTreeRichnessEndemic<<"\t"<<alphaShrubRichnessEndemic<<"\t"<<alphaHerbRichnessEndemic<<"\t"<<Ktree(patchesPosition[i])
								<<"\t"<<cellAbundTreeA<<"\t"<<cellAbundShrubA<<"\t"<<cellAbundHerbA<<
								"\t"<<cellAbundTreeR<<"\t"<<cellAbundShrubR<<"\t"<<cellAbundHerbR<<
								"\t"<<cellAbundTreeS<<"\t"<<cellAbundShrubS<<"\t"<<cellAbundHerbS<<"\n";
						if(invasivesBegin>0){
							exportd8<<replicate<<"\t"<<timeStep<<"\t"<<(patchesPosition[i].first*1000) + patchesPosition[i].second
								<<"\t"<<alphaTreeRichnessI<<"\t"<<alphaShrubRichnessI<<"\t"<<alphaHerbRichnessI
								<<"\t"<<cellAbundTreeAI<<"\t"<<cellAbundShrubAI<<"\t"<<cellAbundHerbAI<<
								"\t"<<cellAbundTreeRI<<"\t"<<cellAbundShrubRI<<"\t"<<cellAbundHerbRI<<
								"\t"<<cellAbundTreeSI<<"\t"<<cellAbundShrubSI<<"\t"<<cellAbundHerbSI<<"\n";
						}
					}

				}
				//Output Gamma richness
				int invasiveTrees=0;
				int invasiveShrubs=0;
				int invasiveHerbs=0;	
				int invasiveTreesAdult=0;
				int invasiveShrubsAdult=0;
				int invasiveHerbsAdult=0;
				int endemicCladoTreesAdult=0;
				int endemicCladoShrubsAdult=0;
				int endemicCladoHerbsAdult=0;
				int endemicCladoTreesRecruit=0;
				int endemicCladoShrubsRecruit=0;
				int endemicCladoHerbsRecruit=0;
				int endemicCladoTreesSeed=0;
				int endemicCladoShrubsSeed=0;
				int endemicCladoHerbsSeed=0;
				int endemicCladoTrees=0;
				int endemicCladoShrubs=0;
				int endemicCladoHerbs=0;	
				int endemicAnaTreesAdult=0;
				int endemicAnaShrubsAdult=0;
				int endemicAnaHerbsAdult=0;
				int endemicAnaTreesRecruit=0;
				int endemicAnaShrubsRecruit=0;
				int endemicAnaHerbsRecruit=0;
				int endemicAnaTreesSeed=0;
				int endemicAnaShrubsSeed=0;
				int endemicAnaHerbsSeed=0;
				int endemicAnaTrees=0;
				int endemicAnaShrubs=0;
				int endemicAnaHerbs=0;	
				int GammaTreeRichness=0;
				int GammaShrubRichness=0;
				int GammaHerbRichness=0;
				int GammaTreeRichnessAdult=0;
				int GammaShrubRichnessAdult=0;
				int GammaHerbRichnessAdult=0;
				int GammaTreeRichnessRecruit=0;
				int GammaShrubRichnessRecruit=0;
				int GammaHerbRichnessRecruit=0;
				int GammaTreeRichnessSeed=0;
				int GammaShrubRichnessSeed=0;
				int GammaHerbRichnessSeed=0;

				for(map<vector<double>, Matrix<float> >::iterator itr = speciesAdult.begin(); itr != speciesAdult.end(); ++itr) 
				{	
					const vector<double>& parameters = itr->first; // key (== parameterSet)
					Matrix<float>& abundances = itr->second; // value (== abundances) 
					//bool empty = matrixIsEmpty(abundances);
					if ((parameters[16]>100000)&&((parameters[25]==0)||(parameters[25]==2)||(parameters[25]==3)||(parameters[25]==4)||(parameters[25]==6))) GammaTreeRichness++;
					if (((parameters[16]>1000) && (parameters[16]<=100000))&&((parameters[25]==0)||(parameters[25]==2)||(parameters[25]==3)||(parameters[25]==4)||(parameters[25]==6))) GammaShrubRichness++;
					if ((parameters[16]<1000)&&((parameters[25]==0)||(parameters[25]==2)||(parameters[25]==3)||(parameters[25]==4)||(parameters[25]==6))) GammaHerbRichness++;	
					if (!matrixIsEmpty(abundances))
					{
						if ((parameters[16]>100000)&&((parameters[25]==0)||(parameters[25]==2)||(parameters[25]==3)||(parameters[25]==4)||(parameters[25]==6))) GammaTreeRichnessAdult++;
						if (((parameters[16]>1000) && (parameters[16]<=100000))&&((parameters[25]==0)||(parameters[25]==2)||(parameters[25]==3)||(parameters[25]==4)||(parameters[25]==6))) GammaShrubRichnessAdult++;
						if ((parameters[16]<1000)&&((parameters[25]==0)||(parameters[25]==2)||(parameters[25]==3)||(parameters[25]==4)||(parameters[25]==6))) GammaHerbRichnessAdult++;	
					}
					if (!matrixIsEmpty(speciesRecruit[parameters]))
					{
						if ((parameters[16]>100000)&&((parameters[25]==0)||(parameters[25]==2)||(parameters[25]==3)||(parameters[25]==4)||(parameters[25]==6))) GammaTreeRichnessRecruit++;
						if (((parameters[16]>1000) && (parameters[16]<=100000))&&((parameters[25]==0)||(parameters[25]==2)||(parameters[25]==3)||(parameters[25]==4)||(parameters[25]==6))) GammaShrubRichnessRecruit++;
						if ((parameters[16]<1000)&&((parameters[25]==0)||(parameters[25]==2)||(parameters[25]==3)||(parameters[25]==4)||(parameters[25]==6))) GammaHerbRichnessRecruit++;	
					}
					if (!matrixIsEmpty(speciesSeed[parameters]))
					{
						if ((parameters[16]>100000)&&((parameters[25]==0)||(parameters[25]==2)||(parameters[25]==3)||(parameters[25]==4)||(parameters[25]==6))) GammaTreeRichnessSeed++;
						if (((parameters[16]>1000) && (parameters[16]<=100000))&&((parameters[25]==0)||(parameters[25]==2)||(parameters[25]==3)||(parameters[25]==4)||(parameters[25]==6))) GammaShrubRichnessSeed++;
						if ((parameters[16]<1000)&&((parameters[25]==0)||(parameters[25]==2)||(parameters[25]==3)||(parameters[25]==4)||(parameters[25]==6))) GammaHerbRichnessSeed++;	
					}
					unsigned int cellnumber =0;
				    
					if((parameters[25]==2)||(parameters[25]==1))
					{
					   exportd11<<replicate<<"\t"<<timeStep<<"\t"
						<<parameters[0]<<"\t"<<parameters[15]<<"\t"<<parameters[19]<<"\t"<<parameters[9]<<
						"\t"<<parameters[23]<<"\t"<<parameters[16]<<"\t"<<parameters[3]<<
						"\t"<<parameters[6]<<"\t"<<parameters[7]<<"\t"<<parameters[14]<<"\t"
						<<parameters[24]<<"\t"<<parameters[10]<<
						"\t"<<parameters[17]*1000+parameters[16]<<"\t"<<parameters[18]<<"\t"<<parameters[17]<<"\n";	
					}

					float sumabundA =0;		
					float sumabundR =0;
					float sumabundS =0;
					
					if(exportd3Tobewritten==1)
					{
						for(unsigned int i=0; i<patchesPosition.size(); i++)
						{
							bool counted = false;
							if(abundances(patchesPosition[i])>0) {
								if(!counted) {
									cellnumber++;	
									counted=true;
								}
								sumabundA += abundances(patchesPosition[i]);
							}

							if(speciesRecruit[parameters](patchesPosition[i])>0) {
								if(!counted) {
									cellnumber++;	
									counted=true;
								}	
								sumabundR += speciesRecruit[parameters](patchesPosition[i]);
							}
							if(speciesSeed[parameters](patchesPosition[i])>0) {
								if(!counted) {
									cellnumber++;	
									counted=true;
								}	
								sumabundS += speciesSeed[parameters](patchesPosition[i]);
							}
						}
					}	
					//counting invasive:
					if (parameters[25]==4){
						if (parameters[16]>100000) invasiveTrees++;
						if ((parameters[16]>1000) && (parameters[16]<=100000)) invasiveShrubs++;
						if (parameters[16]<1000) invasiveHerbs++;
						if (!matrixIsEmpty(abundances))
						{
							if (parameters[16]>100000) invasiveTreesAdult++;
							if ((parameters[16]>1000) && (parameters[16]<=100000)) invasiveShrubsAdult++;
							if (parameters[16]<1000) invasiveHerbsAdult++;	
						}
					}
					//counting cladogenetic endemics:
					if((parameters[25]==2)||(parameters[25]==6)) 
					{
						if (parameters[16]>100000) endemicCladoTrees++;
						if ((parameters[16]>1000) && (parameters[16]<=100000)) endemicCladoShrubs++;
						if (parameters[16]<1000) endemicCladoHerbs++;
						if (!matrixIsEmpty(abundances))
						{
							if (parameters[16]>100000) endemicCladoTreesAdult++;
							if ((parameters[16]>1000) && (parameters[16]<=100000)) endemicCladoShrubsAdult++;
							if (parameters[16]<1000) endemicCladoHerbsAdult++;	
						}
						if (!matrixIsEmpty(speciesRecruit[parameters]))
						{
							if (parameters[16]>100000) endemicCladoTreesRecruit++;
							if ((parameters[16]>1000) && (parameters[16]<=100000)) endemicCladoShrubsRecruit++;
							if (parameters[16]<1000) endemicCladoHerbsRecruit++;	
						}
						if (!matrixIsEmpty(speciesSeed[parameters]))
						{
							if (parameters[16]>100000) endemicCladoTreesSeed++;
							if ((parameters[16]>1000) && (parameters[16]<=100000)) endemicCladoShrubsSeed++;
							if (parameters[16]<1000) endemicCladoHerbsSeed++;	
						}
						if(exportd3Tobewritten==1)
						{
							exportd4<<replicate<<"\t"<<timeStep<<"\t"<<parameters[17]<<
							"\t"<<parameters[0]<<"\t"<<parameters[25]<<"\t"<<parameters[15]<<"\t"<<parameters[23]<<"\t"<<parameters[16]<<
							"\t"<<parameters[3]<<
							"\t"<<parameters[6]<<"\t"<<parameters[7]<<"\t"<<parameters[14]<<"\t"
							<<cellnumber<<"\t"<<parameters[10]<<"\t"<<parameters[19]<<"\t"<<parameters[9]<<
							"\t"<<parameters[17]*1000+parameters[16]<<"\t"<<parameters[18]<<"\t"<<parameters[24]<<"\n";	
						}
					
					}		
					//counting anagenetic endemics:
					if(parameters[25]==3) 
					{
						if (parameters[16]>100000) endemicAnaTrees++;
						if ((parameters[16]>1000) && (parameters[16]<=100000)) endemicAnaShrubs++;
						if (parameters[16]<1000) endemicAnaHerbs++;
						if (!matrixIsEmpty(abundances))
						{
							if (parameters[16]>100000) endemicAnaTreesAdult++;
							if ((parameters[16]>1000) && (parameters[16]<=100000)) endemicAnaShrubsAdult++;
							if (parameters[16]<1000) endemicAnaHerbsAdult++;	
						}
						if (!matrixIsEmpty(speciesRecruit[parameters]))
						{
							if (parameters[16]>100000) endemicAnaTreesRecruit++;
							if ((parameters[16]>1000) && (parameters[16]<=100000)) endemicAnaShrubsRecruit++;
							if (parameters[16]<1000) endemicAnaHerbsRecruit++;	
						}
						if (!matrixIsEmpty(speciesSeed[parameters]))
						{
							if (parameters[16]>100000) endemicAnaTreesSeed++;
							if ((parameters[16]>1000) && (parameters[16]<=100000)) endemicAnaShrubsSeed++;
							if (parameters[16]<1000) endemicAnaHerbsSeed++;	
						}
						if(exportd3Tobewritten==1)
						{
							exportd4<<replicate<<"\t"<<timeStep<<"\t"<<parameters[17]<<
							"\t"<<parameters[0]<<"\t"<<parameters[25]<<"\t"<<parameters[15]<<"\t"<<parameters[23]<<"\t"<<parameters[16]<<
							"\t"<<parameters[3]<<
							"\t"<<parameters[6]<<"\t"<<parameters[7]<<"\t"<<parameters[14]<<"\t"
							<<cellnumber<<"\t"<<parameters[10]<<"\t"<<parameters[19]<<"\t"<<parameters[9]<<
							"\t"<<parameters[17]*1000+parameters[16]<<"\t"<<parameters[18]<<"\t"<<parameters[24]<<"\n";	
						}
					
					}
					if(exportd3Tobewritten==1)
					{
						if (sumabundA!=0) sumabundA=sumabundA/cellnumber;
						if (sumabundR!=0) sumabundR=sumabundR/cellnumber;
						if (sumabundS!=0) sumabundS=sumabundS/cellnumber;

						exportd3<<replicate<<"\t"<<timeStep<<"\t"
						<<parameters[0]<<"\t"<<parameters[25]<<"\t"<<parameters[15]<<"\t"<<parameters[19]<<"\t"<<parameters[9]<<
						"\t"<<parameters[23]<<"\t"<<parameters[16]<<"\t"<<parameters[3]<<
						"\t"<<parameters[6]<<"\t"<<parameters[7]<<"\t"<<parameters[14]<<"\t"
						<<cellnumber<<"\t"<<sumabundA<<"\t"<<sumabundR<<"\t"<<sumabundS<<"\t"<<parameters[24]<<"\t"<<parameters[10]<<
						"\t"<<parameters[17]*1000+parameters[16]<<"\t"<<parameters[18]<<"\t"<<parameters[17]<<"\n";	
					}
				}
						
				exportd<<replicate<<"\t"<<timeStep<<"\t"<<patchesPosition.size()<<"\t"<<GammaTreeRichness<<"\t"<<GammaShrubRichness<<"\t"<<GammaHerbRichness<<"\t"<<
					GammaTreeRichnessAdult<<"\t"<<GammaShrubRichnessAdult<<"\t"<<GammaHerbRichnessAdult<<"\t"<<Col
					<<"\t"<<DemExt<<"\t"<<RandomExt<<"\t"<<VulcanicExt<<"\t"<<ErosionExt<<"\t"<<SpecClado<<"\t"<<SpecAna<<"\t"  
					<<endemicCladoTrees<<"\t"<<endemicCladoShrubs<<"\t"<<endemicCladoHerbs<<"\t"
					<<endemicCladoTreesAdult<<"\t"<<endemicCladoShrubsAdult<<"\t"<<endemicCladoHerbsAdult<<"\t"
					<<endemicAnaTrees<<"\t"<<endemicAnaShrubs<<"\t"<<endemicAnaHerbs<<"\t"
					<<endemicAnaTreesAdult<<"\t"<<endemicAnaShrubsAdult<<"\t"<<endemicAnaHerbsAdult<<"\n";
			
				if(invasivesBegin>0){
					exportd10<<replicate<<"\t"<<timeStep<<"\t"<<patchesPosition.size()<<"\t"<<Col
					<<"\t"<<DemExt<<"\t"<<RandomExt<<"\t"<<VulcanicExt<<"\t"<<ErosionExt<<"\t"<<SpecClado<<"\t"<<SpecAna<<"\t"  
					<<invasiveTrees<<"\t"<<invasiveShrubs<<"\t"<<invasiveHerbs<<"\t"
					<<invasiveTreesAdult<<"\t"<<invasiveShrubsAdult<<"\t"<<invasiveHerbsAdult<<"\n";
				}

				exportd6<<replicate<<"\t"<<timeStep<<"\t"<<patchesPosition.size()<<"\t"<<
					GammaTreeRichnessRecruit<<"\t"<<GammaShrubRichnessRecruit<<"\t"<<GammaHerbRichnessRecruit<<"\t"<<
					GammaTreeRichnessSeed<<"\t"<<GammaShrubRichnessSeed<<"\t"<<GammaHerbRichnessSeed<<"\t"
					<<endemicCladoTreesRecruit<<"\t"<<endemicCladoShrubsRecruit<<"\t"<<endemicCladoHerbsRecruit<<"\t"
					<<endemicCladoTreesSeed<<"\t"<<endemicCladoShrubsSeed<<"\t"<<endemicCladoHerbsSeed<<"\t"
					<<endemicAnaTreesRecruit<<"\t"<<endemicAnaShrubsRecruit<<"\t"<<endemicAnaHerbsRecruit<<"\t"
					<<endemicAnaTreesSeed<<"\t"<<endemicAnaShrubsSeed<<"\t"<<endemicAnaHerbsSeed<<"\n";

				if((vulcanicDisturbance==1)||(randomDisturbance==1)){
					exportd7<<replicate<<"\t"<<timeStep<<"\t"<<patchesPosition.size()<<"\t"<<disturbOccurrence<<"\t"<<				
						disturbRadius<<"\t"<<epicenterX<<"\t"<<epicenterY<<"\n";
					disturbRadius=0;
					disturbOccurrence=0;
					epicenterX=0;
					epicenterY=0;
				}
            
			}
		
			clock_t duration = clock() - start;
			cout << "Replicate " << replicate << " : " << duration << " ms" << endl;
			replicate++;			
			exportd.close();
			exportd2.close();
			exportd3.close();
			exportd4.close();		
			exportd6.close();
			exportd7.close();
		}
	}

}

void Simulation::initialization()
{
		/**********************************/
		/*   Initialization               */  
		/**********************************/
	PI = 3.141592;

    ///////////////////////// Trees:
	//read parameter values fomr hte parameter file:
	std::ifstream input(parameterRangesTree.c_str());
	std::string line;

	std::istringstream linestream;

	// ignore the first line (header)
	std::getline(input, line);

	while(std::getline(input, line)) {
		linestream.clear();
		linestream.str(line);

		double min, max;
		unsigned int useLog;
		float allometabolConst;
		std::string paramName;

		linestream >> paramName >> min >> max >> useLog >> allometabolConst;
		if(useLog) {
            useLogForParam.push_back(true);
			minValueTree.push_back(log(min));
			maxValueTree.push_back(log(max));
			allometabolConstTree.push_back(allometabolConst);
		}
		else {
			// do not use log 
			useLogForParam.push_back(false);
			minValueTree.push_back(min);
			maxValueTree.push_back(max);
			allometabolConstTree.push_back(allometabolConst);
		}
	}

	input.close();

	/////////////////////////Shrubs:
	//read parameter values fomr the parameter file:
	std::ifstream input2(parameterRangesShrub.c_str());
	std::string line2;

	std::istringstream linestream2;

	// ignore the first line (header)
	std::getline(input2, line2);

	while(std::getline(input2, line2)) {
		linestream2.clear();
		linestream2.str(line2);

		double min, max;
		unsigned int useLog;
		float allometabolConst;
		std::string paramName;

		linestream2 >> paramName >> min >> max >> useLog >> allometabolConst;
		if(useLog) {
			useLogForParam.push_back(true);
			minValueShrub.push_back(log(min));
			maxValueShrub.push_back(log(max));
			allometabolConstShrub.push_back(allometabolConst);
		}
		else {
			// do not use log
			useLogForParam.push_back(false);
			minValueShrub.push_back(min);
			maxValueShrub.push_back(max);
			allometabolConstShrub.push_back(allometabolConst);
		}
	}

	input2.close();

	///////////////////// Herbs:
	//read parameter values fomr hte parameter file:
	std::ifstream input3(parameterRangesHerb.c_str());
	std::string line3;

	std::istringstream linestream3;

	// ignore the first line (header)
	std::getline(input3, line3);

	while(std::getline(input3, line3)) {
		linestream3.clear();
		linestream3.str(line3);

		double min, max;
		unsigned int useLog;
		float allometabolConst;
		std::string paramName;

		linestream3 >> paramName >> min >> max >> useLog >> allometabolConst;
		if(useLog) {
			useLogForParam.push_back(true);
			minValueHerb.push_back(log(min));
			maxValueHerb.push_back(log(max));
			allometabolConstHerb.push_back(allometabolConst);
		}
		else {
			// do not use log 
			useLogForParam.push_back(false);
			minValueHerb.push_back(min);
			maxValueHerb.push_back(max);
			allometabolConstHerb.push_back(allometabolConst);
		}
	}

	input3.close();
      /////////////////////////////////////////

    //init of interactions
	timeStep = 0;

	//
	growthStep=1;
	erosionStep=1;
	erosionStep2=maxIslandRadius;

	/*// initialize State*/ 
	//-99: water; 1: island; 0: mainland
	State=gisdata2;

	//For disturbance output:
	disturbRadius=0;
	disturbOccurrence=0;
	epicenterX=0;
	epicenterY=0;

	//Initialize our Carrying Capacity and Temperature matrices
	Ktree = Matrix<float> (NX,NY);  //set up as 0 automatically
	TempKelvin = Matrix<float> (NY,NY);//set up as 0 automatically
	TempKelvinBareIsland = Matrix<float> (NY,NY);//set up as 0 automatically
	lowlandDegradationMatrix = Matrix<float> (NY,NY);
	for (unsigned int r=0;r<NY;++r)
	{
		for (unsigned int c=0;c<NY;++c)
		{
			if (altitudinalTemp=="on"){          //15-15-
				TempKelvin(r,c)= temperature - (State(r,c)-1); 
				TempKelvinBareIsland(r,c)= temperature - (State(r,c)-1); 
			}
			else{
				TempKelvin(r,c)= temperature; 
				TempKelvinBareIsland(r,c)= temperature; 
			}
			if((State(r,c)>0)&&(State(r,c)<= (int)(maxIslandRadius-1)*0.5)) 
			{
				float Degradation = (rand()% degradationImpact)*0.01; 
				lowlandDegradationMatrix(r,c) = Degradation;
			}
		}
	}

	 

	for(unsigned int i=0; i<patchesPosition.size(); i++) 
	{
		Ktree(patchesPosition[i])= maxK;
		
	}
	//Central coordinates: only important for concentric habitat suitability:
	centralx= ((unsigned int)NY-1)*0.5;
	centraly= ((unsigned int)NY-1)*0.5;

	speciesAdult.clear();
	speciesRecruit.clear();
	speciesSeed.clear();
	speciesDisp.clear();
	speciesHabt.clear();
		
	//Mutant maps:
	speciesAdultMutant.clear();
	speciesRecruitMutant.clear();
	speciesSeedMutant.clear();
	speciesDispMutant.clear();
    speciesHabtMutant.clear();

	//Gene Flow map:
	immGeneflow.clear();

	//Initializing Mainland:
	//One Mainland per simulation for all its replicates:
	if(replicate==1)
	{
		mainlandSpeciesPool.clear();
		speciesHabtMain.clear();
		speciesDispMain.clear();
		speciesInitialHabtMain.clear();

				for (unsigned int i=0; i< mainlandSpeciesNumber; i++)
				{
					//Saving the parameter values into a vector, which will be used as a map key:
					vector<double> parameterSet;
					//first, set a PFT: 1=Tree, 2=Shrub, 3=Herb:
					float PFT = (rand()% 1000 + 1 )/(float)1000.0;

					if (PFT <=pftTree)
					{
						writeParameterSet(parameterSet,minValueTree,maxValueTree,allometabolConstTree);
					}
					if ((PFT >pftTree) && (PFT <=(pftTree+pftShrub)))
					{
						writeParameterSet(parameterSet,minValueShrub,maxValueShrub,allometabolConstShrub);
						if (allometabolic==0) {
							parameterSet[8]=parameterSet[8]/maxValueShrub[13];
							parameterSet[22]=parameterSet[8]*0.5;
						}
					}
					if ((PFT>(pftTree+pftShrub)) && (PFT <=1))
					{
						writeParameterSet(parameterSet,minValueHerb,maxValueHerb,allometabolConstHerb);
						if (allometabolic==0) {
							parameterSet[8]=parameterSet[8]/maxValueHerb[13];
							parameterSet[22]=parameterSet[8]*0.5;
						}
					}
					parameterSet[17]=i; //lineage
					//include the species into the vector of vectors (vector of spp):
					mainlandSpeciesPool.push_back(parameterSet);
					speciesHabtMain[parameterSet] = Matrix<float>(NY,NY);
					createDispersalKernel(parameterSet,speciesDispMain);
					habitatSuitability(parameterSet,speciesHabtMain);
					//cout << "speciesHabtMain (centralx,centraly) " << speciesHabtMain[parameterSet](centralx,centraly) <<endl;
					//cout << "speciesHabtMain(centralx,centraly)-1 " << speciesHabtMain[parameterSet](centralx-1,centraly) <<endl;
					//cout << "speciesHabtMain(centralx,centraly)-2 " << speciesHabtMain[parameterSet](centralx-2,centraly) <<endl;
					//cout << "speciesHabtMain(centralx,centraly)-3 " << speciesHabtMain[parameterSet](centralx-3,centraly) <<endl;
					//cout << "speciesHabtMain(centralx,centraly)-4 " << speciesHabtMain[parameterSet](centralx-4,centraly) <<endl;
					//cout << "speciesHabtMain(centralx,centraly)-5 " << speciesHabtMain[parameterSet](centralx-5,centraly) <<endl;
					//cout << "Disp" << speciesDispMain[parameterSet] <<endl;
					//cout << "mainlandSpeciesPool Size:  " << mainlandSpeciesPool.size() <<endl;
					//cout << "mainlandSpecies i:  " << mainlandSpeciesPool[i][0] << "\t"<< mainlandSpeciesPool[i][1] <<"\t"<< mainlandSpeciesPool[i][2] <<
					//	 "\t"<< mainlandSpeciesPool[i][3] <<"\t"<< mainlandSpeciesPool[i][4] <<"\t"<< mainlandSpeciesPool[i][5] <<endl;
				}
				if(readMainlandSpeciespool=="read"){   
				    cout<<"speciesDispMain size"<< speciesDispMain.size()<< endl;
					//read parameter values fomr hte parameter file:
					std::ifstream input2("MainlandSpeciesPoolRead.txt");
					std::string line2;

					std::istringstream linestream2;

					// ignore the first line (header)
					std::getline(input2, line2);
					mainlandSpeciesPool.clear();
					speciesHabtMain.clear();
					speciesDispMain.clear();
					vector<double> parameterSet2;
					for (unsigned int i=0; i< mainlandSpeciesNumber; i++)
					{
					while(std::getline(input2, line2)) {
						linestream2.clear();
						linestream2.str(line2);

						double Fenology,Rmax,M,Gamma,GenTime,MutRate,alpha,P,Kpftfactor,randomTrend,islandSide,Growth,Mrecruit,Germ,Mseed,Annual,Mass,Lineage,MotherspID,
							   optimalAltitude,MassSeedling,MassSeed,KpftfactorSeedling,generalist,verticalgeneralist,speciesType,yearOfSpec;  

						linestream2 >> Fenology>>Rmax>>M>>Gamma>>GenTime>>MutRate>>alpha>>P>>Kpftfactor>>randomTrend>>islandSide>>Growth>>Mrecruit>>Germ>>Mseed>>Annual>>Mass>>
								   Lineage>>MotherspID>>optimalAltitude>>MassSeedling>>MassSeed>>KpftfactorSeedling>>generalist>>verticalgeneralist>>speciesType>>yearOfSpec;
						parameterSet2.clear();
						parameterSet2.push_back(Fenology);
						parameterSet2.push_back(Rmax);
						parameterSet2.push_back(M);
						parameterSet2.push_back(Gamma);
						parameterSet2.push_back(GenTime);
						parameterSet2.push_back(MutRate);
						parameterSet2.push_back(alpha);
						parameterSet2.push_back(P-ldd);
						parameterSet2.push_back(Kpftfactor);
						parameterSet2.push_back(randomTrend);
						parameterSet2.push_back(islandSide);
						parameterSet2.push_back(Growth);
						parameterSet2.push_back(Mrecruit);
						parameterSet2.push_back(Germ);
						parameterSet2.push_back(Mseed);
						parameterSet2.push_back(Annual);
						parameterSet2.push_back(Mass);
						parameterSet2.push_back(Lineage);
						parameterSet2.push_back(MotherspID);
						parameterSet2.push_back(optimalAltitude);
						parameterSet2.push_back(MassSeedling);
						parameterSet2.push_back(MassSeed);
						parameterSet2.push_back(KpftfactorSeedling);
						parameterSet2.push_back(generalist);
						parameterSet2.push_back(verticalgeneralist);
						parameterSet2.push_back(speciesType);
						parameterSet2.push_back(yearOfSpec);
						mainlandSpeciesPool.push_back(parameterSet2);
						speciesHabtMain[parameterSet2] = Matrix<float>(NY,NY);
						createDispersalKernel(parameterSet2,speciesDispMain);
						habitatSuitability(parameterSet2,speciesHabtMain);
					}
					    
					input2.close();
				}
                cout<<"speciesDispMain size"<< speciesDispMain.size()<< endl;
			  }
		cout<<"speciesInitialHabtMain size"<< speciesInitialHabtMain.size()<< endl;
		speciesInitialHabtMain=speciesHabtMain;
		cout<<"speciesInitialHabtMain size"<< speciesInitialHabtMain.size()<< endl;
		if(invasivesBegin>0) {
			       	//read parameter values fomr hte parameter file:
					std::ifstream input3("MainlandSpeciesPoolInvasive.txt");
					std::string line3;

					std::istringstream linestream3;

					// ignore the first line (header)
					std::getline(input3, line3);
					invasiveSpeciesPool.clear();
					vector<double> parameterSet3;
					for (unsigned int i=0; i< invasiveSpeciesNumber; i++)
					{
					while(std::getline(input3, line3)) {
						linestream3.clear();
						linestream3.str(line3);

						double Fenology,Rmax,M,Gamma,GenTime,MutRate,alpha,P,Kpftfactor,randomTrend,islandSide,Growth,Mrecruit,Germ,Mseed,Annual,Mass,Lineage,MotherspID,
							   optimalAltitude,MassSeedling,MassSeed,KpftfactorSeedling,generalist,verticalgeneralist,speciesType,yearOfSpec;  

						linestream3 >> Fenology>>Rmax>>M>>Gamma>>GenTime>>MutRate>>alpha>>P>>Kpftfactor>>randomTrend>>islandSide>>Growth>>Mrecruit>>Germ>>Mseed>>Annual>>Mass>>
								   Lineage>>MotherspID>>optimalAltitude>>MassSeedling>>MassSeed>>KpftfactorSeedling>>generalist>>verticalgeneralist>>speciesType>>yearOfSpec;
						parameterSet3.clear();
						parameterSet3.push_back(Fenology);
						parameterSet3.push_back(Rmax);
						parameterSet3.push_back(M);
						parameterSet3.push_back(Gamma);
						parameterSet3.push_back(mainlandSpeciesNumber*2); //identifyer as invasive (ex-GenTime)
						parameterSet3.push_back(MutRate);
						parameterSet3.push_back(alpha);
						parameterSet3.push_back(P);
						parameterSet3.push_back(Kpftfactor);
						parameterSet3.push_back(randomTrend);
						parameterSet3.push_back(islandSide);
						parameterSet3.push_back(Growth);
						parameterSet3.push_back(Mrecruit);
						parameterSet3.push_back(Germ);
						parameterSet3.push_back(Mseed);
						parameterSet3.push_back(Annual);
						parameterSet3.push_back(Mass);
						parameterSet3.push_back(Lineage);
						parameterSet3.push_back(-1); //invasive (another identifyer!)
						parameterSet3.push_back(optimalAltitude);
						parameterSet3.push_back(MassSeedling);
						parameterSet3.push_back(MassSeed);
						parameterSet3.push_back(KpftfactorSeedling);
						parameterSet3.push_back(generalist);
						parameterSet3.push_back(verticalgeneralist);
						parameterSet3.push_back(4); //invasive 
						parameterSet3.push_back(yearOfSpec);
						invasiveSpeciesPool.push_back(parameterSet3);
						speciesHabtMain[parameterSet3] = Matrix<float>(NY,NY);
						createDispersalKernel(parameterSet3,speciesDispMain);
						habitatSuitability(parameterSet3,speciesHabtMain);
					}					    
					input3.close();
				}
		   }
	}
	else {
		speciesHabtMain=speciesInitialHabtMain;
			if(invasivesBegin>0) {
			    //read parameter values fomr hte parameter file:
				std::ifstream input3("MainlandSpeciesPoolInvasive.txt");
				std::string line3;

				std::istringstream linestream3;

				// ignore the first line (header)
				std::getline(input3, line3);
				invasiveSpeciesPool.clear();
				vector<double> parameterSet3;
				for (unsigned int i=0; i< invasiveSpeciesNumber; i++) //for simplicity, number of invasive species pool is the same as of mainland
				{
					while(std::getline(input3, line3)) {
					linestream3.clear();
					linestream3.str(line3);

					double Fenology,Rmax,M,Gamma,GenTime,MutRate,alpha,P,Kpftfactor,randomTrend,islandSide,Growth,Mrecruit,Germ,Mseed,Annual,Mass,Lineage,MotherspID,
							   optimalAltitude,MassSeedling,MassSeed,KpftfactorSeedling,generalist,verticalgeneralist,speciesType,yearOfSpec;  

					linestream3 >> Fenology>>Rmax>>M>>Gamma>>GenTime>>MutRate>>alpha>>P>>Kpftfactor>>randomTrend>>islandSide>>Growth>>Mrecruit>>Germ>>Mseed>>Annual>>Mass>>
								Lineage>>MotherspID>>optimalAltitude>>MassSeedling>>MassSeed>>KpftfactorSeedling>>generalist>>verticalgeneralist>>speciesType>>yearOfSpec;
					parameterSet3.clear();
					parameterSet3.push_back(Fenology);
					parameterSet3.push_back(Rmax);
					parameterSet3.push_back(M);
					parameterSet3.push_back(Gamma);
					parameterSet3.push_back(mainlandSpeciesNumber*2);
					parameterSet3.push_back(MutRate);
					parameterSet3.push_back(alpha);
					parameterSet3.push_back(P);
					parameterSet3.push_back(Kpftfactor);
					parameterSet3.push_back(randomTrend);
					parameterSet3.push_back(islandSide);
					parameterSet3.push_back(Growth);
					parameterSet3.push_back(Mrecruit);
					parameterSet3.push_back(Germ);
					parameterSet3.push_back(Mseed);
					parameterSet3.push_back(Annual);
					parameterSet3.push_back(Mass);
					parameterSet3.push_back(Lineage);
					parameterSet3.push_back(-1);//invasive
					parameterSet3.push_back(optimalAltitude);
					parameterSet3.push_back(MassSeedling);
					parameterSet3.push_back(MassSeed);
					parameterSet3.push_back(KpftfactorSeedling);
					parameterSet3.push_back(generalist);
					parameterSet3.push_back(verticalgeneralist);
					parameterSet3.push_back(4);//invasive
					parameterSet3.push_back(yearOfSpec);
					invasiveSpeciesPool.push_back(parameterSet3);
					speciesHabtMain[parameterSet3] = Matrix<float>(NY,NY);
					createDispersalKernel(parameterSet3,speciesDispMain);
					habitatSuitability(parameterSet3,speciesHabtMain);
				}					    
				input3.close();
			}
		}
	}
}
void Simulation::colonizationForest()
{
	//distruting seeds of each species of the mainland:
	speciesInfoCol.clear();
	speciesDispCol.clear();
	speciesHabtCol.clear();
	if(singleSpecies=="true"){
			//Saving the parameter values into a vector, which will be used as a map key:
			vector<double> parameterSet;
			parameterSet= mainlandSpeciesPool[currentSpecies];
			addSpeciesIfNew(parameterSet, speciesInfoCol, speciesDispCol, speciesHabtCol,speciesDispMain, speciesHabtMain);			
	}
	if(singleSpecies=="false"){
		 for (unsigned int i=0; i<mainlandSpeciesNumber; i++)
		 {		 
			//Saving the parameter values into a vector, which will be used as a map key:
			vector<double> parameterSet;
			parameterSet= mainlandSpeciesPool[i];
			addSpeciesIfNew(parameterSet, speciesInfoCol, speciesDispCol, speciesHabtCol,speciesDispMain, speciesHabtMain);
		 }
	}
	    
		//Dispersal into the island:
	    dispersalfromMainland(speciesInfoCol, speciesDispCol, Ktree, speciesHabtCol);
		//Take the species that could not colonize out of the colonizer map
		removeNoncolonizingSpecies(speciesInfoCol, speciesDispCol, speciesHabtCol);
		//increment the islander's map
		unsigned int Initialrichness =speciesAdult.size();
		Col=0;
		if(speciesInfoCol.size()>0)
		{
          addSpeciesIfcolonizing(speciesInfoCol,speciesAdult, speciesRecruit, speciesSeed, Ktree,speciesDispCol,speciesDisp,speciesHabtCol,speciesHabt);
		}
        if(speciesAdult.size()>Initialrichness) Col=speciesAdult.size()-Initialrichness;
}
void Simulation::colonization()
{
	//Picking random "Species" whithin the mainland species pool that will be given the oportunity
	//to disperse into the islands
	//unsigned int migrationEvents = 10;
	speciesInfoCol.clear();
	speciesDispCol.clear();
	speciesHabtCol.clear();
	if(singleSpecies=="true"){
			//Saving the parameter values into a vector, which will be used as a map key:
			vector<double> parameterSet;
			parameterSet= mainlandSpeciesPool[currentSpecies];
			addSpeciesIfNew(parameterSet, speciesInfoCol, speciesDispCol, speciesHabtCol,speciesDispMain, speciesHabtMain);			
	}
	if(singleSpecies=="false"){
		 for (unsigned int i=0; i< migrationEvents; i++)
		 {		 
			//Saving the parameter values into a vector, which will be used as a map key:
			vector<double> parameterSet;		
			int randomSpfromMainLandPool = rand()%( mainlandSpeciesNumber);
			parameterSet= mainlandSpeciesPool[randomSpfromMainLandPool];
			addSpeciesIfNew(parameterSet, speciesInfoCol, speciesDispCol, speciesHabtCol,speciesDispMain, speciesHabtMain);	
		 }
	}
	    
		//Dispersal into the island:
	    dispersalfromMainland(speciesInfoCol, speciesDispCol, Ktree, speciesHabtCol);
		//Take the species that could not colonize out of the colonizer map
		removeNoncolonizingSpecies(speciesInfoCol, speciesDispCol, speciesHabtCol);
		//increment the islander's map
		//cout  << "speciesInfoCol size: " <<speciesInfoCol.size() <<endl;
		unsigned int Initialrichness =speciesAdult.size();
		Col=0;
		if(speciesInfoCol.size()>0)
		{
          addSpeciesIfcolonizing(speciesInfoCol,speciesAdult, speciesRecruit, speciesSeed, Ktree,speciesDispCol,speciesDisp,speciesHabtCol,speciesHabt);
		}
        if(speciesAdult.size()>Initialrichness) Col=speciesAdult.size()-Initialrichness;
}
void Simulation::invasiveColonization()
{
	//Picking random "Species" whithin the mainland species pool that will be given the oportunity
	//to disperse into the islands
	speciesInfoCol.clear();
	speciesDispCol.clear();
	speciesHabtCol.clear();
     for (unsigned int i=0; i< invasivesperTimestep; i++)
	 {		 
		//Saving the parameter values into a vector, which will be used as a map key:
		vector<double> parameterSet;
		int randomSpfromMainLandPool = rand()%(invasiveSpeciesNumber);
		parameterSet= invasiveSpeciesPool[randomSpfromMainLandPool];
		addSpeciesIfNew(parameterSet, speciesInfoCol, speciesDispCol, speciesHabtCol,speciesDispMain, speciesHabtMain);

	 }
	    
		//Dispersal into the area:
	    dispersalfromMainland(speciesInfoCol, speciesDispCol, Ktree, speciesHabtCol);
		//Take the species that could not colonize out of the colonizer map
		removeNoncolonizingSpecies(speciesInfoCol, speciesDispCol, speciesHabtCol);

		unsigned int Initialrichness =speciesAdult.size();
		Col=0;
		if(speciesInfoCol.size()>0)
		{
          addSpeciesIfcolonizing(speciesInfoCol,speciesAdult, speciesRecruit, speciesSeed, Ktree,speciesDispCol,speciesDisp,speciesHabtCol,speciesHabt);
		}
        if(speciesAdult.size()>Initialrichness) Col=speciesAdult.size()-Initialrichness;
}
void Simulation::writeParameterSet(vector<double>& parameterSet,vector<double>& minValue,vector<double>& maxValue,vector<double>& allometabolConst)
{
	if(allometabolic==0){
		//Relative Fenology:
		double Fenology = rand()% 32000;

		//Rmax:
		double Rmax = rand()%( (int)((maxValue[1]-minValue[1]) *1000));
		Rmax = (minValue[1]*1000) + Rmax;
		Rmax = Rmax*0.001;
		//M:
		double M = rand()%( (int)((maxValue[2]-minValue[2]) *1000));
		M = (minValue[2]*1000) + M;
		M = M*0.001;
		//Gamma version 1:
		double Gamma = rand()%( (int)((maxValue[3]-minValue[3])*1000));
		Gamma = (minValue[3]*1000) + Gamma;
		Gamma = Gamma*0.001;
		//GenTime:
		double GenTime = rand()%( (int)((maxValue[4]-minValue[4]) *1000));
		GenTime = (minValue[4]*1000) + GenTime;
		GenTime = GenTime*0.001;
		//MutRate:
		double MutRate = rand()%( (int)((maxValue[5]-minValue[5]) *1000000000000000));
		MutRate = (minValue[5]*1000000000000000) + MutRate;
		MutRate = MutRate*0.000000000000001;
				
		//Alpha (mean dispersal distance - in cells):
		double alpha = rand()%( (int)((maxValue[6]-minValue[6]) *1000000));
		alpha = (minValue[6]*1000000) + alpha;
		alpha = alpha*0.000001;
		//cout << "alpha: " << alpha <<endl;

		//double alpha= (rand()% 50)*0.001 + 0.00001;
		//P (for 2Dt dispersal kernel):
		double P = rand()%( (int)((maxValue[7]-minValue[7]) *10000));
		P = (minValue[7]*10000) + P;
		P = P*0.0001;
		double Kpftfactor = 1.0;
		double KpftfactorSeedling = Kpftfactor;		

		//suitability on the optimal altitudinal cell layer:
		double randomTrend = (rand()% 50 +1)*0.01 + 0.5;
		double islandSide= (rand()% 4) +1;

		//growth probability
		double Growth = rand()%( (int)((maxValue[8]-minValue[8]) *1000));
		Growth = (minValue[8]*1000) + Growth;
		Growth = Growth*0.001;
		//Recruit M:
		double Mrecruit = rand()%( (int)((maxValue[9]-minValue[9]) *1000));
		Mrecruit = (minValue[9]*1000) + Mrecruit;
		Mrecruit = Mrecruit*0.001;
		//Germination probability
		double Germ = rand()%( (int)((maxValue[10]-minValue[10]) *1000));
		Germ = (minValue[10]*1000) + Germ;
		Germ = Germ*0.001;
		//Seed M:
		double Mseed = rand()%( (int)((maxValue[11]-minValue[11]) *1000));
		Mseed = (minValue[11]*1000) + Mseed;
		Mseed = Mseed*0.001;
		double Annual = rand()% 2;	
		double Lineage = 0;
		double MotherspID = 0;

		//For hsuitability:
		double optimalAltitude = (rand()% maxOptaltitude) + minOptaltitude;

		//specialist: 1; generalist: 4:
		double generalist = (rand()% 4) +1;
		//Vertical specialist: 1; generalist: 4:
		double verticalgeneralist = (rand()% verticalpreference) +1;
	
		double Mass, MassSeedling, MassSeed;
		//allometabolConst[0] gives the PFT (1=tree, 2=shrub, 3=herb):
		if(allometabolConst[0]==1) 
		{
			int randomnumber = rand() % 100 +1;
		Mass = randomnumber*100000;
		}
		if(allometabolConst[0]==2)
		{
		int randomnumber = rand() % 100 +1;
		Mass = randomnumber*1000;
		}
		if(allometabolConst[0]==3)
		{
		int randomnumber = rand() % 100 +1;
		Mass = randomnumber*10;
		Annual = rand()% 2;
		}
		if (Mass>55000)  {
			MassSeedling = Mass*0.2;
			MassSeed = (exp((rand()% 80)*0.1))*0.01; 			
		}
		else {
			MassSeedling = Mass*0.2;
			MassSeed= pow(Mass*0.0000046,0.5);//0.00043*Mass;			
		}

		parameterSet.push_back(Fenology);
		parameterSet.push_back(Rmax);
		parameterSet.push_back(M);
		parameterSet.push_back(Gamma);
		parameterSet.push_back(GenTime);
		parameterSet.push_back(MutRate);
		parameterSet.push_back(alpha);
		parameterSet.push_back(P);
		parameterSet.push_back(Kpftfactor);
		parameterSet.push_back(randomTrend);
		parameterSet.push_back(islandSide);
		parameterSet.push_back(Growth);
		parameterSet.push_back(Mrecruit);
		parameterSet.push_back(Germ);
		parameterSet.push_back(Mseed);
		parameterSet.push_back(Annual);
		parameterSet.push_back(Mass); // key =Mass for the metabolic version
		parameterSet.push_back(Lineage);
		parameterSet.push_back(MotherspID);
		parameterSet.push_back(optimalAltitude);
		parameterSet.push_back(MassSeedling); // key =MassSeedling for hte metabolic version
		parameterSet.push_back(MassSeed); // key =MassSeed for hte metabolic version
		parameterSet.push_back(KpftfactorSeedling);
		parameterSet.push_back(generalist);
		parameterSet.push_back(verticalgeneralist);// key =Mass for hte metabolic version
		parameterSet.push_back(0);//0= species from the initial species pool
		parameterSet.push_back(0); //time step of initialization

	}

	else
	{
		double Mass, MassSeedling, MassSeed;
		//Annual (1) or not (0):
		double Annual = 0;
				
		//allometabolConst[0] gives the PFT (1=tree, 2=shrub, 3=herb):
		if(allometabolConst[0]==1) 
		{
			int randomnumber = rand() % 100 +1;
		Mass = randomnumber*100000;
		}
		if(allometabolConst[0]==2)
		{
		int randomnumber = rand() % 100 +1;
		Mass = randomnumber*1000;
		}
		if(allometabolConst[0]==3)
		{
		int randomnumber = rand() % 100 +1;
		Mass = randomnumber*10;
		Annual = rand()% 2;
		}
		//Relative Fenology:
		double Fenology = rand()% 32000;
		
		//For hsuitability and metabolic rates:
		double optimalAltitude = (rand()% (maxOptaltitude-minOptaltitude)) + minOptaltitude;
		//suitability on the optimal altitudinal cell layer:
		double randomTrend = (rand()% 101)*0.01;	
		//optimal mountain/island side:
		double islandSide= (rand()% 4) +1;
		//island side amplitude: specialist: 1; generalist: 4:
		double generalist = (rand()% 4) +1;
		//Vertical specialist: 1; generalist: 4:
		double verticalgeneralist = (rand()% verticalpreference) +1;
		//Alpha (mean dispersal distance - in cells):
		double alpha= (rand()% (int) (maxValue[6]*1000))*0.001 + minValue[6];
		//P (for 2Dt dispersal kernel):
		double P = rand()%( (int)((maxValue[7]-minValue[7]) *100));
		P = (minValue[7]*100) + P;
		P = P*0.01;
		
		
		//Rmax:  
		double Rmax = (1/pow(Mass,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConst[1];
		//M:
		double M = (1/pow(Mass,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConst[2];
		//GenTime:
		double GenTime = rand()%( (int)((maxValue[4]-minValue[4]) *1000));
		GenTime = (minValue[4]*1000) + GenTime;
		GenTime = GenTime*0.001;
		//MutRate:
		double MutRate;
		MutRate =  (1/pow(Mass,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConst[5];
		
		double Kpftfactor, KpftfactorSeedling;
		//From Niklas (1994) Plant Allometry:
		if(allometabolConst[0]==1) //allometabolConst[0]==1 => tree
		{  
		 double diameter = (1/63.0)*pow(Mass,0.4);
         Kpftfactor = (pow(diameter,2))*0.01;
		 KpftfactorSeedling = (pow(diameter*0.5,2))*0.01;
		}
		else
		{
		 double diameter = (1/46.41589)*pow(Mass,0.3333);
		 Kpftfactor = (pow(diameter,2))*0.01;
		 KpftfactorSeedling = (pow(diameter*0.5,2))*0.01;
		}				
		//with the same metabolic constant than adults, but lower body mass:
		//Recruit M:
		double Mrecruit = 0;
		//Seed M:
		double Mseed =0;
		//growth probability
		double Growth =0;
		//Germination probability
		double Germ = 0;
		if (Mass>55000)  {
			MassSeedling = Mass*0.2;
			MassSeed = (exp((rand()% 80)*0.1))*0.01; 
			Mrecruit =  (1/pow(MassSeedling,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConst[9];
			Mseed =  (1/pow(MassSeed,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConst[11];//tree
			Growth =  (1/pow(MassSeedling,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConst[8];
			Germ =  (1/pow(MassSeed,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConst[10];
		}
		else {
			MassSeedling = Mass*0.2;
			MassSeed= pow(Mass*0.0000046,0.5);//0.00043*Mass;
			Mrecruit =  (1/pow(MassSeedling,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConst[9];
			Mseed =  (1/pow(MassSeed,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConst[11];//shrub
			Growth =  (1/pow(MassSeedling,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConst[8];
			Germ =  (1/pow(MassSeed,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConst[10];
		}

	
		double Lineage = 0;
		double MotherspID = 0;

		//Gamma: Version 1 //Can this be scaled with tamperature and mass?????
		double Gamma = rand()%( (int)((maxValue[3]-minValue[3])*1000));
		Gamma = minValue[3]*1000 + Gamma;
		Gamma = Gamma*0.001;
		if (Gamma>0) Gamma=rand()% 1000 + 1;
		    

		parameterSet.push_back(Fenology);
		parameterSet.push_back(Rmax);
		parameterSet.push_back(M);
		parameterSet.push_back(Gamma);
		parameterSet.push_back(GenTime);
		parameterSet.push_back(MutRate);
		parameterSet.push_back(alpha);
		parameterSet.push_back(P);
		parameterSet.push_back(Kpftfactor);
		parameterSet.push_back(randomTrend);
		parameterSet.push_back(islandSide);
		parameterSet.push_back(Growth);
		parameterSet.push_back(Mrecruit);
		parameterSet.push_back(Germ);
		parameterSet.push_back(Mseed);
		parameterSet.push_back(Annual);
		parameterSet.push_back(Mass);
		parameterSet.push_back(Lineage);
		parameterSet.push_back(MotherspID);
		parameterSet.push_back(optimalAltitude);
		parameterSet.push_back(MassSeedling);
		parameterSet.push_back(MassSeed);
		parameterSet.push_back(KpftfactorSeedling);
		parameterSet.push_back(generalist);
		parameterSet.push_back(verticalgeneralist);
		parameterSet.push_back(0); //0= species from the initial species pool
		parameterSet.push_back(0); //time step of initialization
	}
}
void Simulation::addSpeciesIfNew(vector<double>& parameterSet, map<vector<double>, Matrix<float> >& spInfo,
	                             map<vector<double>, Matrix<float> >& spDisp,
								 map<vector<double>, Matrix<float> >& spHabt, map<vector<double>, Matrix<float> >& spDispMain,
								 map<vector<double>, Matrix<float> >& spHabtMain) 
{
	// look if sp contains parameterSet
	bool contained = spInfo.find(parameterSet) != spInfo.end();

	if(!contained) {
		spInfo[parameterSet] = Matrix<float>(NY,NY);
		spDisp[parameterSet] = spDispMain[parameterSet];
		spHabt[parameterSet] = spHabtMain[parameterSet];
	}
}
void Simulation::addRandomSpeciesIfNew(vector<double>& parameterSet, map<vector<double>, Matrix<float> >& spInfo,
	                             map<vector<double>, Matrix<float> >& spDisp,
								 map<vector<double>, Matrix<float> >& spHabt) 
{
	// look if sp contains parameterSet
	bool contained = spInfo.find(parameterSet) != spInfo.end();
	if(!contained) {
		spInfo[parameterSet] = Matrix<float>(NY,NY);
		spHabt[parameterSet] = Matrix<float>(NY,NY);
        createDispersalKernel(parameterSet,spDisp);
		habitatSuitability(parameterSet,spHabt);
	}
}
void Simulation::habitatSuitability(vector<double>& parameterSet, map<vector<double>, Matrix<float> >& spHabt)
{
	vector<pair<unsigned int,unsigned int> > positionsToProcess;
	Matrix<float>& suitability = spHabt[parameterSet];	
	//cout <<"TempKelvin(centralx,centraly) during suitability update: "<< TempKelvin(centralx,centraly) << endl;
	//Autocorrelated random variant:
	if (habSuitability=="Random")
	{
		for (unsigned int x=1;x<NY-1;++x)
		{
			for (unsigned int y=1;y<NY-1;++y)
			{
				positionsToProcess.push_back(std::pair<int,int>(x,y));				
				suitability(x,y) = 0; // 0 -> not initialized
				//uncorrelated suitability:
	//			double randomSuitability = ((rand()% 100) +1)/100.0;
	//			spHabt[parameterSet](r,c) = randomSuitability;
			}
		}
		while(positionsToProcess.size() > 0) {
			// select a random index (between 0 and coords.size() - 1)
			int randomCoord = rand() % positionsToProcess.size();
			unsigned int randomX = positionsToProcess[randomCoord].first;
			unsigned int randomY = positionsToProcess[randomCoord].second;
			float randomSuitability = (float) ((rand()% 100 +1)*0.01);
			for (unsigned int x=randomX-1;x<=randomX+1;++x)
			{
				for (unsigned int y=randomY-1;y<=randomY+1;++y)
				{
					if(suitability(x,y)==0) 
					{
						// cell was not processed yet
						suitability(x,y) = randomSuitability;

						// remove from positionsToProcess if not in the edge
						if((x>0)&&(x<NY-1)&&(y>0)&&(y<NY-1)) {

							// get an interator on the cell coords in positionsToProcess
							vector<pair<unsigned int,unsigned int> >::iterator itr
								= std::find(positionsToProcess.begin(), positionsToProcess.end(), std::pair<unsigned int, unsigned int>(x,y));

							positionsToProcess.erase(itr);
						}
					}				
				}
			}		     

		}
	}

	//cout << "habSuitability " << habSuitability <<endl;
	//Autocorrelated concentric variant:
	checkedSuitability = Matrix<int> (NY,NY);
	positionsToProcess.clear();
	
	if (habSuitability=="Concentric")
	{
			//decrease in suitability depending on the distance:
			double decreaseSuitperDist = parameterSet[9]/parameterSet[24]; 
			int stepToProcess = 0;

			if (dynamicGeology=="yes") {
				if(growthStep<=(maxIslandRadius+1)) stepToProcess =growthStep-1;
				if(erosionStep>=2) stepToProcess =growthStep-erosionStep;
			}			
			else stepToProcess = mountainRadius; //mountain radius

			if((int)parameterSet[23]==1){
				if((int)parameterSet[10]==1){
				//updating to the current island size:
					int counting =0;
					for (unsigned int x=centralx-stepToProcess;x<=centralx;++x)
					{				
						for (unsigned int y=centraly-stepToProcess+counting;y<=centraly+stepToProcess-counting;++y)
						{
							//modular distance to the optimal altitude:
							int tempPref = 0;
							if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
							else tempPref = (int)parameterSet[19];
							int modulDist=0;
							if(fixedIslandaltitude=="no"){
								 modulDist= (x-(centralx-stepToProcess))- tempPref;
							}
							else modulDist = 1 - tempPref; //1= lowland
							
							if(modulDist<0) modulDist = -modulDist;	
							suitability(x,y) = (float) (parameterSet[9]) - (float) decreaseSuitperDist*modulDist;
							if (suitability(x,y)<0.0001) suitability(x,y)=0;

						}
						counting++;
					}
				}
				if((int)parameterSet[10]==2){
					//updating to the current island size:
					int counting =0;
					for (unsigned int x=centralx+stepToProcess;x>=centralx;--x)
					{
						
						for (unsigned int y=centraly+stepToProcess-counting;y>=centraly-stepToProcess+counting;--y)
						{
							//modular distance to the optimal altitude:
							int tempPref = 0;
							if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
							else tempPref = (int)parameterSet[19];
							int modulDist=0;
							if(fixedIslandaltitude=="no"){
								 modulDist= ((centralx+stepToProcess)-x)- tempPref;
							}
							else modulDist = 1 - tempPref; //1= lowland

							if(modulDist<0) modulDist = -modulDist;	
							suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
							if (suitability(x,y)<0.0001) suitability(x,y)=0;
							
						}
						counting++;
					}
				}
				if((int)parameterSet[10]==3){
					//updating to the current island size:
					int counting =0;
					for (unsigned int y=centraly-stepToProcess;y<=centraly;++y)
					{
						
						for (unsigned int x=centralx-stepToProcess+counting;x<=centralx+stepToProcess-counting;++x)
						{
							//modular distance to the optimal altitude:
							int tempPref = 0;
							if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
							else tempPref = (int)parameterSet[19];
							int modulDist=0;
							if(fixedIslandaltitude=="no"){
								 modulDist= (y-(centraly-stepToProcess))- tempPref;
							}
							else modulDist = 1 - tempPref; //1= lowland
							
							if(modulDist<0) modulDist = -modulDist;	
							suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
							if (suitability(x,y)<0.0001) suitability(x,y)=0;

						}
						counting++;
					}
				}
				if((int)parameterSet[10]==4){
					//updating to the current island size:
					int counting =0;
					for (unsigned int y=centraly+stepToProcess;y>=centraly;--y)
					{
						
						for (unsigned int x=centralx+stepToProcess-counting;x>=centralx-stepToProcess+counting;--x)
						{
							//modular distance to the optimal altitude:
							int tempPref = 0;
							if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
							else tempPref = (int)parameterSet[19];
							int modulDist=0;
							if(fixedIslandaltitude=="no"){
								 modulDist= ((centraly+stepToProcess)-y)- tempPref;
							}
							else modulDist = 1 - tempPref; //1= lowland

							if(modulDist<0) modulDist = -modulDist;	
							suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
							if (suitability(x,y)<0.0001) suitability(x,y)=0;
						}
						counting++;
					}
				}
			}


			if((int)parameterSet[23]==2){
				int neighbour = rand()% 1; //0 to the right, 1 to the left
				if((int)parameterSet[10]==1){

				//updating to the current island size:
					int counting =0;
					for (unsigned int x=centralx-stepToProcess;x<=centralx;++x)
					{						
				
						for (unsigned int y=centraly-stepToProcess+counting;y<=centraly+stepToProcess-counting;++y)
						{
							//modular distance to the optimal altitude:
							int tempPref = 0;
							if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
							else tempPref = (int)parameterSet[19];
							int modulDist=0;
							if(fixedIslandaltitude=="no"){
								 modulDist= (x-(centralx-stepToProcess))- tempPref;
							}
							else modulDist = 1 - tempPref; //1= lowland

							if(modulDist<0) modulDist = -modulDist;	
							suitability(x,y) = (float) (parameterSet[9]) - (float) decreaseSuitperDist*modulDist;
							if (suitability(x,y)<0.0001) suitability(x,y)=0;
						}
						counting++;
					}
					counting =0;
					if (neighbour==0){
						for (unsigned int y=centraly+stepToProcess;y>=centraly;--y)
						{
							
							for (unsigned int x=centralx+stepToProcess-counting;x>=centralx-stepToProcess+counting;--x)
							{
								//modular distance to the optimal altitude:
								int tempPref = 0;
								if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
								else tempPref = (int)parameterSet[19];
								int modulDist=0;
								if(fixedIslandaltitude=="no"){
									 modulDist= ((centraly+stepToProcess)-y)- tempPref;
								}
								else modulDist = 1 - tempPref; //1= lowland 

								if(modulDist<0) modulDist = -modulDist;	
								suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
								if (suitability(x,y)<0.0001) suitability(x,y)=0;
							}
							counting++;
						}
					}
					else{
						for (unsigned int y=centraly-stepToProcess;y<=centraly;++y)
						{
							
							for (unsigned int x=centralx-stepToProcess+counting;x<=centralx+stepToProcess-counting;++x)
							{
								//modular distance to the optimal altitude:
								int tempPref = 0;
								if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
								else tempPref = (int)parameterSet[19];
								int modulDist=0;
								if(fixedIslandaltitude=="no"){
									 modulDist= (y-(centraly-stepToProcess))- tempPref;
								}
								else modulDist = 1 - tempPref; //1= lowland 

								if(modulDist<0) modulDist = -modulDist;	
								suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
								if (suitability(x,y)<0.0001) suitability(x,y)=0;
							}
							counting++;
						}
					}
				}
				if((int)parameterSet[10]==2){
					//updating to the current island size:
					int counting =0;
					for (unsigned int x=centralx+stepToProcess;x>=centralx;--x)
					{
						
						for (unsigned int y=centraly+stepToProcess-counting;y>=centraly-stepToProcess+counting;--y)
						{
							//modular distance to the optimal altitude:
							int tempPref = 0;
							if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
							else tempPref = (int)parameterSet[19];
							int modulDist=0;
							if(fixedIslandaltitude=="no"){
									modulDist= ((centralx+stepToProcess)-x)- tempPref;
							}
							else modulDist = 1 - tempPref; //1= lowland 

							if(modulDist<0) modulDist = -modulDist;	
							suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
							if (suitability(x,y)<0.0001) suitability(x,y)=0;
						}
						counting++;
					}
					counting =0;
					if (neighbour==0){
						for (unsigned int y=centraly-stepToProcess;y<=centraly;++y)
						{
							
							for (unsigned int x=centralx-stepToProcess+counting;x<=centralx+stepToProcess-counting;++x)
							{
								//modular distance to the optimal altitude:
								int tempPref = 0;
								if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
								else tempPref = (int)parameterSet[19];
								int modulDist=0;
								if(fixedIslandaltitude=="no"){
										modulDist=(y-(centraly-stepToProcess))- tempPref;
								}
								else modulDist = 1 - tempPref; //1= lowland 
								
								if(modulDist<0) modulDist = -modulDist;	
								suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
								if (suitability(x,y)<0.0001) suitability(x,y)=0;
							}
							counting++;
						}
					}
					else{
						for (unsigned int y=centraly+stepToProcess;y>=centraly;--y)
						{
							
							for (unsigned int x=centralx+stepToProcess-counting;x>=centralx-stepToProcess+counting;--x)
							{
								//modular distance to the optimal altitude:
								int tempPref = 0;
								if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
								else tempPref = (int)parameterSet[19];
								int modulDist=0;
								if(fixedIslandaltitude=="no"){
										modulDist=((centraly+stepToProcess)-y)- tempPref;
								}
								else modulDist = 1 - tempPref; //1= lowland 

								if(modulDist<0) modulDist = -modulDist;	
								suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
								if (suitability(x,y)<0.0001) suitability(x,y)=0;
							}
							counting++;
						}
					}

				}
				if((int)parameterSet[10]==3){
					//updating to the current island size:
					int counting =0;
					for (unsigned int y=centraly-stepToProcess;y<=centraly;++y)
					{
						
						for (unsigned int x=centralx-stepToProcess+counting;x<=centralx+stepToProcess-counting;++x)
						{
							//modular distance to the optimal altitude:
							int tempPref = 0;
							if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
							else tempPref = (int)parameterSet[19];
							int modulDist=0;
							if(fixedIslandaltitude=="no"){
									modulDist=(y-(centraly-stepToProcess))- tempPref;
							}
							else modulDist = 1 - tempPref; //1= lowland 

							if(modulDist<0) modulDist = -modulDist;	
							suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
							if (suitability(x,y)<0.0001) suitability(x,y)=0;
						}
						counting++;
					}
					counting=0; 
					if (neighbour==0){
						for (unsigned int x=centralx-stepToProcess;x<=centralx;++x)
						{
							
				
							for (unsigned int y=centraly-stepToProcess+counting;y<=centraly+stepToProcess-counting;++y)
							{
								//modular distance to the optimal altitude:
								int tempPref = 0;
								if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
								else tempPref = (int)parameterSet[19];
								int modulDist=0;
								if(fixedIslandaltitude=="no"){
										modulDist=(x-(centralx-stepToProcess))- tempPref;
								}
								else modulDist = 1 - tempPref; //1= lowland

								if(modulDist<0) modulDist = -modulDist;	
								suitability(x,y) = (float) (parameterSet[9]) - (float) decreaseSuitperDist*modulDist;
								if (suitability(x,y)<0.0001) suitability(x,y)=0;
							}
							counting++;
						}
					}
					else{
						for (unsigned int x=centralx+stepToProcess;x>=centralx;--x)
						{
								
							for (unsigned int y=centraly+stepToProcess-counting;y>=centraly-stepToProcess+counting;--y)
							{
								//modular distance to the optimal altitude:
								int tempPref = 0;
								if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
								else tempPref = (int)parameterSet[19];
								int modulDist=0;
								if(fixedIslandaltitude=="no"){
										modulDist=((centralx+stepToProcess)-x)- tempPref;
								}
								else modulDist = 1 - tempPref; //1= lowland
								
								if(modulDist<0) modulDist = -modulDist;
								suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
								if (suitability(x,y)<0.0001) suitability(x,y)=0;
							}
							counting++;
						}
					}
				}
				if((int)parameterSet[10]==4){
					//updating to the current island size:
					int counting =0;
					for (unsigned int y=centraly+stepToProcess;y>=centraly;--y)
					{
						
						for (unsigned int x=centralx+stepToProcess-counting;x>=centralx-stepToProcess+counting;--x)
						{
							//modular distance to the optimal altitude:
							int tempPref = 0;
							if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
							else tempPref = (int)parameterSet[19];
							int modulDist=0;
							if(fixedIslandaltitude=="no"){
									modulDist=((centraly+stepToProcess)-y)- tempPref;
							}
							else modulDist = 1 - tempPref; //1= lowland

							if(modulDist<0) modulDist = -modulDist;	
							suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
							if (suitability(x,y)<0.0001) suitability(x,y)=0;
						}
						counting++;
					}
					counting=0; 
					if (neighbour==0){
						for (unsigned int x=centralx+stepToProcess;x>=centralx;--x)
						{
							
							for (unsigned int y=centraly+stepToProcess-counting;y>=centraly-stepToProcess+counting;--y)
							{
								//modular distance to the optimal altitude:
								int tempPref = 0;
								if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
								else tempPref = (int)parameterSet[19];
								int modulDist=0;
								if(fixedIslandaltitude=="no"){
										modulDist=((centralx+stepToProcess)-x)- tempPref;
								}
								else modulDist = 1 - tempPref; //1= lowland

								if(modulDist<0) modulDist = -modulDist;	
								suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
								if (suitability(x,y)<0.0001) suitability(x,y)=0;
							}
							counting++;
						}
					}
					else{
						for (unsigned int x=centralx-stepToProcess;x<=centralx;++x)
						{
							
				
							for (unsigned int y=centraly-stepToProcess+counting;y<=centraly+stepToProcess-counting;++y)
							{
								//modular distance to the optimal altitude:
								int tempPref = 0;
								if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
								else tempPref = (int)parameterSet[19];
								int modulDist=0;
								if(fixedIslandaltitude=="no"){
										modulDist=(x-(centralx-stepToProcess))- tempPref;
								}
								else modulDist = 1 - tempPref; //1= lowland

								if(modulDist<0) modulDist = -modulDist;	
								suitability(x,y) = (float) (parameterSet[9]) - (float) decreaseSuitperDist*modulDist;
								if (suitability(x,y)<0.0001) suitability(x,y)=0;
							}
							counting++;
						}
					}
				}
			}	


			if((int)parameterSet[23]==3){
				if((int)parameterSet[10]==1){
				//updating to the current island size:
					int counting =0;
					for (unsigned int x=centralx-stepToProcess;x<=centralx;++x)
					{
						
				
						for (unsigned int y=centraly-stepToProcess+counting;y<=centraly+stepToProcess-counting;++y)
						{
							//modular distance to the optimal altitude:
							int tempPref = 0;
							if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
							else tempPref = (int)parameterSet[19];
							int modulDist=0;
							if(fixedIslandaltitude=="no"){
									modulDist=(x-(centralx-stepToProcess))- tempPref;
							}
							else modulDist = 1 - tempPref; //1= lowland

							if(modulDist<0) modulDist = -modulDist;	
							suitability(x,y) = (float) (parameterSet[9]) - (float) decreaseSuitperDist*modulDist;
							if (suitability(x,y)<0.0001) suitability(x,y)=0;
						}
						counting++;
					}
					counting =0;
					for (unsigned int y=centraly+stepToProcess;y>=centraly;--y)
					{
						
						for (unsigned int x=centralx+stepToProcess-counting;x>=centralx-stepToProcess+counting;--x)
						{
							//modular distance to the optimal altitude:
							int tempPref = 0;
							if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
							else tempPref = (int)parameterSet[19];
							int modulDist=0;
							if(fixedIslandaltitude=="no"){
									modulDist=((centraly+stepToProcess)-y)- tempPref;
							}
							else modulDist = 1 - tempPref; //1= lowland

							if(modulDist<0) modulDist = -modulDist;	
							suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
							if (suitability(x,y)<0.0001) suitability(x,y)=0;
						}
						counting++;
					}
					counting=0; 
					for (unsigned int y=centraly-stepToProcess;y<=centraly;++y)
					{
						
						for (unsigned int x=centralx-stepToProcess+counting;x<=centralx+stepToProcess-counting;++x)
						{
							//modular distance to the optimal altitude:
							int tempPref = 0;
							if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
							else tempPref = (int)parameterSet[19];
							int modulDist=0;
							if(fixedIslandaltitude=="no"){
									modulDist=(y-(centraly-stepToProcess))- tempPref;
							}
							else modulDist = 1 - tempPref; //1= lowland

							if(modulDist<0) modulDist = -modulDist;	
							suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
							if (suitability(x,y)<0.0001) suitability(x,y)=0;
						}
						counting++;
					}
				}
				if((int)parameterSet[10]==2){
					//updating to the current island size:
					int counting =0;
					for (unsigned int x=centralx+stepToProcess;x>=centralx;--x)
					{
						
						for (unsigned int y=centraly+stepToProcess-counting;y>=centraly-stepToProcess+counting;--y)
						{
							//modular distance to the optimal altitude:
							int tempPref = 0;
							if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
							else tempPref = (int)parameterSet[19];
							int modulDist=0;
							if(fixedIslandaltitude=="no"){
									modulDist=((centralx+stepToProcess)-x)- tempPref;
							}
							else modulDist = 1 - tempPref; //1= lowland

							if(modulDist<0) modulDist = -modulDist;	
							suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
							if (suitability(x,y)<0.0001) suitability(x,y)=0;
						}
						counting++;
					}
					counting =0;
					for (unsigned int y=centraly-stepToProcess;y<=centraly;++y)
					{
						
						for (unsigned int x=centralx-stepToProcess+counting;x<=centralx+stepToProcess-counting;++x)
						{
							//modular distance to the optimal altitude:
							int tempPref = 0;
							if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
							else tempPref = (int)parameterSet[19];
							int modulDist=0;
							if(fixedIslandaltitude=="no"){
									modulDist=(y-(centraly-stepToProcess))- tempPref;
							}
							else modulDist = 1 - tempPref; //1= lowland

							if(modulDist<0) modulDist = -modulDist;	
							suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
							if (suitability(x,y)<0.0001) suitability(x,y)=0;
						}
						counting++;
					}
					counting=0; 
					for (unsigned int y=centraly+stepToProcess;y>=centraly;--y)
					{
						
						for (unsigned int x=centralx+stepToProcess-counting;x>=centralx-stepToProcess+counting;--x)
						{
							//modular distance to the optimal altitude:
							int tempPref = 0;
							if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
							else tempPref = (int)parameterSet[19];
							int modulDist=0;
							if(fixedIslandaltitude=="no"){
									modulDist=((centraly+stepToProcess)-y)- tempPref;
							}
							else modulDist = 1 - tempPref; //1= lowland

							if(modulDist<0) modulDist = -modulDist;	
							suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
							if (suitability(x,y)<0.0001) suitability(x,y)=0;
						}
						counting++;
					}
				}
				if((int)parameterSet[10]==3){
					//updating to the current island size:
					int counting =0;
					for (unsigned int y=centraly-stepToProcess;y<=centraly;++y)
					{
						
						for (unsigned int x=centralx-stepToProcess+counting;x<=centralx+stepToProcess-counting;++x)
						{
							//modular distance to the optimal altitude:
							int tempPref = 0;
							if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
							else tempPref = (int)parameterSet[19];
							int modulDist=0;
							if(fixedIslandaltitude=="no"){
									modulDist=(y-(centraly-stepToProcess))- tempPref;
							}
							else modulDist = 1 - tempPref; //1= lowland

							if(modulDist<0) modulDist = -modulDist;	
							suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
							if (suitability(x,y)<0.0001) suitability(x,y)=0;
						}
						counting++;
					}
					counting=0; 			
					for (unsigned int x=centralx-stepToProcess;x<=centralx;++x)
					{
						
				
						for (unsigned int y=centraly-stepToProcess+counting;y<=centraly+stepToProcess-counting;++y)
						{
							//modular distance to the optimal altitude:
							int tempPref = 0;
							if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
							else tempPref = (int)parameterSet[19];
							int modulDist=0;
							if(fixedIslandaltitude=="no"){
									modulDist=(x-(centralx-stepToProcess))- tempPref;
							}
							else modulDist = 1 - tempPref; //1= lowland

							if(modulDist<0) modulDist = -modulDist;	
							suitability(x,y) = (float) (parameterSet[9]) - (float) decreaseSuitperDist*modulDist;
							if (suitability(x,y)<0.0001) suitability(x,y)=0;
						}
						counting++;
					}
					counting=0; 
					for (unsigned int x=centralx+stepToProcess;x>=centralx;--x)
					{
						
						for (unsigned int y=centraly+stepToProcess-counting;y>=centraly-stepToProcess+counting;--y)
						{
							//modular distance to the optimal altitude:
							int tempPref = 0;
							if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
							else tempPref = (int)parameterSet[19];
							int modulDist=0;
							if(fixedIslandaltitude=="no"){
									modulDist=((centralx+stepToProcess)-x)- tempPref;
							}
							else modulDist = 1 - tempPref; //1= lowland

							if(modulDist<0) modulDist = -modulDist;	
							suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
							if (suitability(x,y)<0.0001) suitability(x,y)=0;
						}
						counting++;
					}
				}
				if((int)parameterSet[10]==4){
					//updating to the current island size:
					int counting =0;
					for (unsigned int y=centraly+stepToProcess;y>=centraly;--y)
					{
						
						for (unsigned int x=centralx+stepToProcess-counting;x>=centralx-stepToProcess+counting;--x)
						{
							//modular distance to the optimal altitude:
							int tempPref = 0;
							if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
							else tempPref = (int)parameterSet[19];
							int modulDist=0;
							if(fixedIslandaltitude=="no"){
									modulDist=((centraly+stepToProcess)-y)- tempPref;
							}
							else modulDist = 1 - tempPref; //1= lowland

							if(modulDist<0) modulDist = -modulDist;	
							suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
							if (suitability(x,y)<0.0001) suitability(x,y)=0;
						}
						counting++;
					}
					counting=0; 				
					for (unsigned int x=centralx+stepToProcess;x>=centralx;--x)
					{
						
						for (unsigned int y=centraly+stepToProcess-counting;y>=centraly-stepToProcess+counting;--y)
						{
							//modular distance to the optimal altitude:
							int tempPref = 0;
							if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
							else tempPref = (int)parameterSet[19];
							int modulDist=0;
							if(fixedIslandaltitude=="no"){
									modulDist=((centralx+stepToProcess)-x)- tempPref;
							}
							else modulDist = 1 - tempPref; //1= lowland

							if(modulDist<0) modulDist = -modulDist;	
							suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
							if (suitability(x,y)<0.0001) suitability(x,y)=0;
						}
						counting++;
					}
					counting=0; 
					for (unsigned int x=centralx-stepToProcess;x<=centralx;++x)
					{
						
				
						for (unsigned int y=centraly-stepToProcess+counting;y<=centraly+stepToProcess-counting;++y)
						{
							//modular distance to the optimal altitude:
							int tempPref = 0;
							if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
							else tempPref = (int)parameterSet[19];
							int modulDist=0;
							if(fixedIslandaltitude=="no"){
									modulDist=(x-(centralx-stepToProcess))- tempPref;
							}
							else modulDist = 1 - tempPref; //1= lowland

							if(modulDist<0) modulDist = -modulDist;	
							suitability(x,y) = (float) (parameterSet[9]) - (float) decreaseSuitperDist*modulDist;
							if (suitability(x,y)<0.0001) suitability(x,y)=0;
						}
						counting++;
					}
				}
			}	


			if((int)parameterSet[23]==4){
			//updating to the current island size:
				int counting =0;
				for (unsigned int x=centralx-stepToProcess;x<=centralx;++x)
				{				
				
					for (unsigned int y=centraly-stepToProcess+counting;y<=centraly+stepToProcess-counting;++y)
					{
						//modular distance to the optimal altitude:
						int tempPref = 0;
						if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
						else tempPref = (int)parameterSet[19];
						int modulDist=0;
						if(fixedIslandaltitude=="no"){
								modulDist=(x-(centralx-stepToProcess))- tempPref;
						}
						else modulDist = 1 - tempPref; //1= lowland

						if(modulDist<0) modulDist = -modulDist;	
						suitability(x,y) = (float) (parameterSet[9]) - (float) decreaseSuitperDist*modulDist;
						if (suitability(x,y)<0.0001) suitability(x,y)=0;
					}
					counting++;
				}
				counting =0;				
				for (unsigned int x=centralx+stepToProcess;x>=centralx;--x)
				{
					
					for (unsigned int y=centraly+stepToProcess-counting;y>=centraly-stepToProcess+counting;--y)
					{
						//modular distance to the optimal altitude:
						int tempPref = 0;
						if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
						else tempPref = (int)parameterSet[19];
						int modulDist=0;
						if(fixedIslandaltitude=="no"){
								modulDist=((centralx+stepToProcess)-x)- tempPref;
						}
						else modulDist = 1 - tempPref; //1= lowland

						if(modulDist<0) modulDist = -modulDist;	
						suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
						if (suitability(x,y)<0.0001) suitability(x,y)=0;
					}
					counting++;
				}
				counting =0;				
				for (unsigned int y=centraly-stepToProcess;y<=centraly;++y)
				{
					for (unsigned int x=centralx-stepToProcess+counting;x<=centralx+stepToProcess-counting;++x)
					{
						//modular distance to the optimal altitude:
						int tempPref = 0;
						if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
						else tempPref = (int)parameterSet[19];
						int modulDist=0;
						if(fixedIslandaltitude=="no"){
								modulDist=(y-(centraly-stepToProcess))-tempPref;
						}
						else modulDist = 1 - tempPref; //1= lowland

						if(modulDist<0) modulDist = -modulDist;	
						suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
						if (suitability(x,y)<0.0001) suitability(x,y)=0;
					}
					counting++;
				}
				counting=0; 				
				for (unsigned int y=centraly+stepToProcess;y>=centraly;--y)
				{	
					for (unsigned int x=centralx+stepToProcess-counting;x>=centralx-stepToProcess+counting;--x)
					{
						//modular distance to the optimal altitude:
						int tempPref = 0;
						if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y));
						else tempPref = (int)parameterSet[19];
						int modulDist=0;
						if(fixedIslandaltitude=="no"){
								modulDist=((centraly+stepToProcess)-y)- tempPref;
						}
						else modulDist = 1 - tempPref; //1= lowland

						if(modulDist<0) modulDist = -modulDist;
						suitability(x,y) = (float) (parameterSet[9]) - (float)decreaseSuitperDist*modulDist;
						if (suitability(x,y)<0.0001) suitability(x,y)=0;
					}
					counting++;
				}
				
			}
			//applying degration onto lowlands:
			if((lowlandDegradation>0)&&(timeStep>= lowlandDegradation)){
				for (unsigned int r=0;r<NY;++r)
				{
					for (unsigned int c=0;c<NX;++c)
					{
						if((State(c,r)>0)&&(State(c,r)<= (int) (maxIslandRadius-1)*0.5 )) 
						{
							suitability(c,r) = suitability(c,r) -  suitability(c,r)*lowlandDegradationMatrix(c,r);
						}
					}
				}	
			}
	}

	if (habSuitability=="Plane")
	{
		//decrease in suitability depending on the distance:
		double decreaseSuitperDist = parameterSet[9]/parameterSet[24]; 
		for (unsigned int y=0;y<NY;++y)
			{
				for (unsigned int x=0;x<NX;++x)
				{	
						if (State(x,y)>=1)
						{
						//modular distance to the optimal altitude:
						int tempPref = 0;
						//this will make tempPref != optimal altitude only in case of feedbacks to temperature (e.g. with trees as eco engineers):
						if (altitudinalTemp=="on") tempPref = (int)parameterSet[19] - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y)); 
						else tempPref = (int)parameterSet[19];
						int modulDist = State(x,y)- tempPref;
						/*if (altitudinalTemp=="on") tempPref = (temperature - (int)parameterSet[19]) - (int) (TempKelvinBareIsland(x,y) - TempKelvin(x,y)); 
						else tempPref =  (temperature - (int)parameterSet[19]);
						int modulDist = TempKelvin(x,y)- tempPref;*/
						if(modulDist<0) modulDist = -modulDist;	
						suitability(x,y) = (float) (parameterSet[9]) - (float) decreaseSuitperDist*modulDist;
						if (suitability(x,y)<0.0001) suitability(x,y)=0;
						}
			}
		}
	}
   
	//cout << "GrowthStep  "<< growthStep << endl;
	/*cout << "Optimal Altitude "<< parameterSet[19] << endl;
	cout << "Max suitability  "<< parameterSet[9] << endl;
	cout << suitability(centralx-3,centraly-3)  <<"  " <<suitability(centralx-2,centraly-3) <<"  " <<suitability(centralx-1,centraly-3)<<"  " <<suitability(centralx,centraly-3) <<"  " <<suitability(centralx+1,centraly-3) <<"  " <<suitability(centralx+2,centraly-3)<<"  " <<suitability(centralx+3,centraly-3)   <<endl;
	cout << suitability(centralx-3,centraly-2)  <<"  " <<suitability(centralx-2,centraly-2) <<"  " <<suitability(centralx-1,centraly-2)<<"  " <<suitability(centralx,centraly-2) <<"  " <<suitability(centralx+1,centraly-2) <<"  " <<suitability(centralx+2,centraly-2)<<"  " <<suitability(centralx+3,centraly-2)   <<endl;
	cout << suitability(centralx-3,centraly-1)  <<"  " <<suitability(centralx-2,centraly-1) <<"  " <<suitability(centralx-1,centraly-1)<<"  " <<suitability(centralx,centraly-1) <<"  " <<suitability(centralx+1,centraly-1) <<"  " <<suitability(centralx+2,centraly-1)<<"  " <<suitability(centralx+3,centraly-1)   <<endl;
	cout << suitability(centralx-3,centraly)  <<"  " <<suitability(centralx-2,centraly) <<"  " <<suitability(centralx-1,centraly)<<"  " <<suitability(centralx,centraly) <<"  " <<suitability(centralx+1,centraly) <<"  " <<suitability(centralx+2,centraly)<<"  " <<suitability(centralx+3,centraly)   <<endl;
	cout << suitability(centralx-3,centraly+1)  <<"  " <<suitability(centralx-2,centraly+1) <<"  " <<suitability(centralx-1,centraly+1)<<"  " <<suitability(centralx,centraly+1) <<"  " <<suitability(centralx+1,centraly+1) <<"  " <<suitability(centralx+2,centraly+1)<<"  " <<suitability(centralx+3,centraly+1)   <<endl;
	cout << suitability(centralx-3,centraly+2)  <<"  " <<suitability(centralx-2,centraly+2) <<"  " <<suitability(centralx-1,centraly+2)<<"  " <<suitability(centralx,centraly+2) <<"  " <<suitability(centralx+1,centraly+2) <<"  " <<suitability(centralx+2,centraly+2)<<"  " <<suitability(centralx+3,centraly+2)   <<endl;
	cout << suitability(centralx-3,centraly+3)  <<"  " <<suitability(centralx-2,centraly+3) <<"  " <<suitability(centralx-1,centraly+3)<<"  " <<suitability(centralx,centraly+3) <<"  " <<suitability(centralx+1,centraly+3) <<"  " <<suitability(centralx+2,centraly+3)<<"  " <<suitability(centralx+3,centraly+3)   <<endl;
			
	cout << "Optimal Altitude "<< parameterSet[19] << endl;
	cout << "Max suitability  "<< parameterSet[9] << endl;
	cout << "Temp amplitude  "<< parameterSet[24] << endl;
	cout << suitability(2,2)  <<"  " <<suitability(3,2) <<"  " <<suitability(4,2)<<endl;
	cout << suitability(2,3)  <<"  " <<suitability(3,3) <<"  " <<suitability(4,3) <<endl;
	cout << suitability(2,4)  <<"  " <<suitability(3,4) <<"  " <<suitability(4,4)<<endl;
	cout << suitability(2,5)  <<"  " <<suitability(3,5) <<"  " <<suitability(4,5) <<endl;
	cout << suitability(2,6)  <<"  " <<suitability(3,6) <<"  " <<suitability(4,6) <<endl;
	cout << suitability(2,7)  <<"  " <<suitability(3,7) <<"  " <<suitability(4,7) <<endl;
	cout << suitability(2,8)  <<"  " <<suitability(3,8) <<"  " <<suitability(4,8)<<endl;
	cout << suitability(2,9)  <<"  " <<suitability(3,9) <<"  " <<suitability(4,9)<<endl;
	cout << suitability(2,10)  <<"  " <<suitability(3,10) <<"  " <<suitability(4,10)<<endl;
	cout << suitability(2,11)  <<"  " <<suitability(3,11) <<"  " <<suitability(4,11)<<endl;
	cout << suitability(2,12)  <<"  " <<suitability(3,12) <<"  " <<suitability(4,12)<<endl;
	cout << suitability(2,13)  <<"  " <<suitability(3,13) <<"  " <<suitability(4,13)<<endl;*/

}
void Simulation::createDispersalKernel(vector<double>& parameterSet, map<vector<double>, Matrix<float> >& spDisp)
{
   unsigned int x = kernelRadius + 1;
   unsigned int y = NY*2;
   double sum=0;
   bool isNegExp = (dispKernel=="NegExp");

   // create a new dispersal Matrix

   spDisp[parameterSet] = Matrix<float>(x,y);

   for(unsigned int i=0;i<x; i++) {
		for(unsigned int j=0;j<x; j++) {
		   double r = sqrt((double)(i*i + j*j));
			 
		   double dispersal;

		   if (isNegExp) {
			   dispersal=dispersalNegExpFunction(parameterSet[6],r);
		   }
		   else {
			   dispersal=dispersal2DtFunction(parameterSet[6],parameterSet[7],r);
		   }

		   if(i==0) {
			   if(j==0) {
					// center
				   sum += dispersal;
			   }
			   else {
					// y axis
				   sum += 2*dispersal;
			   }
		   }
		   else {
			   if(j==0) {
					// x axis
				   sum += 2*dispersal;
			   }
			   else {
					// normal point
				   sum += 4* dispersal;
			   }
		   }

		   if(j < y) {
				spDisp[parameterSet](i,j) = (float) dispersal;
		   }
		}
   }
   //cout << "Sum: " << sum <<endl;
  //cout << "float spDisp[parameterSet](kernelRadius-6,NY-6) : " << (float) spDisp[parameterSet](kernelRadius-6,NY-6) <<endl;
   //double sum2;
   //cout <<"spDisp[parameters](18,18):" << spDisp[parameterSet](18,18)<< endl;
   for(unsigned int i=0;i<x; i++) {
		for(unsigned int j=0;j<y; j++) {
			 spDisp[parameterSet](i,j)=spDisp[parameterSet](i,j)/ (float)sum;
			 //sum2 += spDisp[parameterSet](i,j);
		}
   }
  //cout <<"parameterSet[6]:" << parameterSet[6]<< endl;
  //cout <<"parameterSet[7]:" << parameterSet[7]<< endl;
  //cout <<"spDisp[parameters](18,18):" << spDisp[parameterSet](18,18)<< endl;
// cout << "spDisp[parameterSet](floor((kernelRadius+1)*0.5)-1,floor((NY+1)*0.5))-1) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5)-1,floor((NY+1)*0.5))-1 <<endl;
// cout << "spDisp[parameterSet](floor((kernelRadius+1)*0.5),floor((NY+1)*0.5))) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5),floor((NY+1)*0.5)) <<endl;
  //cout << "spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,floor((NY+1)*0.5))+1) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,floor((NY+1)*0.5))+1 <<endl;   
/* cout << "spDisp[parameterSet](0,0) : " << spDisp[parameterSet](0,0) <<endl;
	 cout << "spDisp[parameterSet](1,1) : " << spDisp[parameterSet](1,1) <<endl;
     cout << "spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,0) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,0) <<endl;
     cout << "spDisp[parameterSet](kernelRadius,NY) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,1) <<endl;
     cout << "spDisp[parameterSet](kernelRadius,NY) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,2) <<endl;
     cout << "spDisp[parameterSet](kernelRadius,NY) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,3) <<endl;
     cout << "spDisp[parameterSet](kernelRadius,NY) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,4) <<endl;
	 cout << "spDisp[parameterSet](kernelRadius,NY) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,5) <<endl;
	 cout << "spDisp[parameterSet](kernelRadius,NY) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,6) <<endl;
	 cout << "spDisp[parameterSet](kernelRadius,NY) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,7) <<endl;
	 cout << "spDisp[parameterSet](kernelRadius,NY) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,8) <<endl;
	 cout << "spDisp[parameterSet](kernelRadius,NY) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,9) <<endl;
	 cout << "spDisp[parameterSet](kernelRadius,NY) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,10) <<endl;
	 cout << "spDisp[parameterSet](kernelRadius,NY) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,11) <<endl;
	 cout << "spDisp[parameterSet](kernelRadius,NY) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,12) <<endl;
	 cout << "spDisp[parameterSet](kernelRadius,NY) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,13) <<endl;
	 cout << "spDisp[parameterSet](kernelRadius,NY) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,14) <<endl;
	 cout << "spDisp[parameterSet](kernelRadius,NY) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,15) <<endl;
	 cout << "spDisp[parameterSet](kernelRadius,NY) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,16) <<endl;
	 cout << "spDisp[parameterSet](kernelRadius,NY) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,17) <<endl;
	 cout << "spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,18) : " << spDisp[parameterSet](floor((kernelRadius+1)*0.5)+1,18) <<endl;
     cout << "spDisp[parameterSet](kernelRadius,NY) : " << spDisp[parameterSet](kernelRadius,NY) <<endl;*/
  //cout << "Sum2 : " << sum2 <<endl;
}
double Simulation::dispersalNegExpFunction(double alpha, double r)
{
   double N = 2*PI*(alpha*alpha);
   double p = (1/N )* exp(- r/alpha);
   if(p<0) p =0;
   return p;
}
double Simulation::dispersal2DtFunction(double alpha, double P, double r)
{
	double u = pow(alpha, 2)*gamma(P);

   double p = P/(pow(PI*u*(1+(pow(r,2)/u)),(P+1))); 
   if(p<0) p =0;
   return p;
}
double Simulation::gamma(double x)
{
    int i,k,m;
    double ga,gr,r,z;

    static double g[] = {
        1.0,
        0.5772156649015329,
       -0.6558780715202538,
       -0.420026350340952e-1,
        0.1665386113822915,
       -0.421977345555443e-1,
       -0.9621971527877e-2,
        0.7218943246663e-2,
       -0.11651675918591e-2,
       -0.2152416741149e-3,
        0.1280502823882e-3,
       -0.201348547807e-4,
       -0.12504934821e-5,
        0.1133027232e-5,
       -0.2056338417e-6,
        0.6116095e-8,
        0.50020075e-8,
       -0.11812746e-8,
        0.1043427e-9,
        0.77823e-11,
       -0.36968e-11,
        0.51e-12,
       -0.206e-13,
       -0.54e-14,
        0.14e-14};

    if (x > 171.0) return 1e308;    // This value is an overflow flag.
    if (x == (int)x) {
        if (x > 0.0) {
            ga = 1.0;               // use factorial
            for (i=2;i<x;i++) {
               ga *= i;
            }
         }
         else
            ga = 1e308;
     }
     else {
        if (fabs(x) > 1.0) {
            z = fabs(x);
            m = (int)z;
            r = 1.0;
            for (k=1;k<=m;k++) {
                r *= (z-k);
            }
            z -= m;
        }
        else
            z = x;
        gr = g[24];
        for (k=23;k>=0;k--) {
            gr = gr*z+g[k];
        }
        ga = 1.0/(gr*z);
        if (fabs(x) > 1.0) {
            ga *= r;
            if (x < 0.0) {
                ga = -M_PI/(x*ga*sin(M_PI*x));
            }
        }
    }
    return ga;
}

void Simulation::dispersalfromMainland(map<vector<double>, Matrix<float> >& spInfo, map<vector<double>, Matrix<float> >& spDisp, Matrix<float>& K, map<vector<double>, Matrix<float> >& spHabt)
{
 	map<vector<double>, Matrix<float> >::iterator itr = spInfo.begin();
    
	if(colonize=="island"){
		while(itr != spInfo.end()) {
			const vector<double>& parameters = itr->first; // key (== parameterSet)
			Matrix<float>& abundances = itr->second; // value (== abundances)
			Matrix<double> incomingDispersalUnits (NY,NY);  
       
			for (unsigned int i=0; i<patchesPositionMainLand.size(); i++ )
			{
			    //Set first a random number of dispersal units (can be local defined):			
				int dispersalUnits = (rand() % 10000)*(rand() % 100) + 10;
				int is = patchesPositionMainLand[i].first;
				int iz = patchesPositionMainLand[i].second;
			    //Use a dispersal kernel to disperse the dispersal units (use a dispersal funtion)
				for(unsigned int j=0; j<patchesPosition.size(); j++)
				{
					int ib = patchesPosition[j].first;
					int ia = patchesPosition[j].second;
					incomingDispersalUnits(ib,ia) += spDisp[parameters](abs(ib-is),abs(ia-iz))*dispersalUnits;
				}
			}

			//use the incoming dispersal units to draw a random recruitment
			for(unsigned int i=0; i<patchesPosition.size(); i++)
			{
			   if(incomingDispersalUnits(patchesPosition[i])>0)
			   {
				   float seeds = Zufall.poisson(incomingDispersalUnits(patchesPosition[i]));

				   if(seeds>0)
				   {
					   abundances(patchesPosition[i])= seeds;
				   }				   	   	   	   	   	   	   	   	   	 
			   }
			}

			++itr;
		}
	}
	if((colonize=="forest")&&(timeStep==0)){
		while(itr != spInfo.end()) {
			const vector<double>& parameters = itr->first; // key (== parameterSet)
			Matrix<float>& abundances = itr->second; // value (== abundances)
     
			for(unsigned int j=0; j<patchesPosition.size(); j++)
			{
				int ib = patchesPosition[j].first;
				int ia = patchesPosition[j].second;
				abundances(ib,ia)= (float) Zufall.poisson(100*(rand()% 400));					
			}

			++itr;
		}
	}

	if((colonize=="forest")&&(timeStep==invasivesBegin)){ // invasion happens only one time step
		while(itr != spInfo.end()) {

			const vector<double>& parameters = itr->first; // key (== parameterSet)
			Matrix<float>& abundances = itr->second; // value (== abundances)
     		
			for(unsigned int i=0; i<randomInvasivecenters; i++)
			{
			   int randomCoordInvasivecenters = rand ()% patchesPosition.size(); 
			   float NumberofInvasives = rand()% 1000 + 1; //Invasion pressure
               abundances(patchesPosition[randomCoordInvasivecenters])= NumberofInvasives;
			 }

			++itr;
		}
	}
}
bool Simulation::matrixIsEmpty(Matrix<float>& matrix) {
	for(unsigned int i = 0; i < matrix.getWidth(); ++i) {
		for(unsigned int j = 0; j < matrix.getHeight(); ++j) {
			//only integer individuals!
			if((matrix(i,j) >= 1.0)&&(State(i,j)>=1)) {
				return false;
			}
		}
	}
	return true;
}

void Simulation::removeNoncolonizingSpecies(map<vector<double>, Matrix<float> >& spInfo, map<vector<double>, Matrix<float> >& spDisp,
								 map<vector<double>, Matrix<float> >& spHabt){
	// we need a "special" loop (with "while") when an element of the map may be erased
	map<vector<double>, Matrix<float> >::iterator itr = spInfo.begin();
	map<vector<double>, Matrix<float> >::iterator it = spDisp.begin();
	map<vector<double>, Matrix<float> >::iterator itH = spHabt.begin();

	while(itr != spInfo.end()) {

		const vector<double>& parameters = itr->first; // key (== parameterSet)
		Matrix<float>& abundances = itr->second; // value (== abundances)
		
		if(matrixIsEmpty(abundances)) {
			// remove it
			itr = spInfo.erase(itr);
		
		}
		else {
			++itr;
		}
	}

	while(it != spDisp.end()) {

		const vector<double>& parameters = it->first; // key (== parameterSet)
		Matrix<float>& abundances = it->second; // value (== abundances)
		
		bool contained = spInfo.find(parameters) != spInfo.end();
		if(!contained) {
			// remove it
			it= spDisp.erase(it);
			
		}
		else {
			++it;
		}
	}

	while(itH != spHabt.end()) {

		const vector<double>& parameters = itH->first; // key (== parameterSet)
		Matrix<float>& abundances = itH->second; // value (== abundances)
		
		bool contained = spInfo.find(parameters) != spInfo.end();
		if(!contained) {
			// remove it
			itH= spHabt.erase(itH);
			
		}
		else {
			++itH;
		}
	}
}

void Simulation::removeExtinctSpecies(map<vector<double>, Matrix<float> >& spAdult, 
	map<vector<double>, Matrix<float> >& spRecruit,	map<vector<double>, Matrix<float> >& spSeed,
	map<vector<double>, Matrix<float> >& spDisp, map<vector<double>, Matrix<float> >& spHabt){
	// we need a "special" loop (with "while") when an element of the map may be erased
	map<vector<double>, Matrix<float> >::iterator itr = spSeed.begin();
	map<vector<double>, Matrix<float> >::iterator itrr = spRecruit.begin();
	map<vector<double>, Matrix<float> >::iterator itra = spAdult.begin();
	map<vector<double>, Matrix<float> >::iterator it = spDisp.begin();
	map<vector<double>, Matrix<float> >::iterator itH = spHabt.begin();
	while(itr != spSeed.end()) {
		const vector<double>& parameters = itr->first; // key (== parameterSet)
		Matrix<float>& abundances = itr->second; // value (== abundances)		
		if((matrixIsEmpty(abundances))&&(matrixIsEmpty(spRecruit[parameters]))&&(matrixIsEmpty(spAdult[parameters]))) {
			// remove it
			if ((parameters[25]==1)||(parameters[25]==5)) prespeciesExtinct++;
			itr = spSeed.erase(itr);			
		}
		else {
			++itr;
		}
	}

	while(itrr != spRecruit.end()) {
		const vector<double>& parameters = itrr->first; // key (== parameterSet)
		Matrix<float>& abundances = itrr->second; // value (== abundances)		
		bool contained = spSeed.find(parameters) != spSeed.end();
		if(!contained) {
			// remove it
			itrr= spRecruit.erase(itrr);			
		}
		else {
			++itrr;
		}
	}

	while(itra != spAdult.end()) {
		const vector<double>& parameters = itra->first; // key (== parameterSet)
		Matrix<float>& abundances = itra->second; // value (== abundances)		
		bool contained = spSeed.find(parameters) != spSeed.end();
		if(!contained) {
			// remove it
			itra= spAdult.erase(itra);			
		}
		else {
			++itra;
		}
	}

	while(it != spDisp.end()) {
		const vector<double>& parameters = it->first; // key (== parameterSet)
		Matrix<float>& abundances = it->second; // value (== abundances)		
		bool contained = spSeed.find(parameters) != spSeed.end();
		if(!contained) {
			// remove it
			it= spDisp.erase(it);			
		}
		else {
			++it;
		}
	}

	while(itH != spHabt.end()) {
		const vector<double>& parameters = itH->first; // key (== parameterSet)
		Matrix<float>& abundances = itH->second; // value (== abundances)		
		bool contained = spSeed.find(parameters) != spSeed.end();
		if(!contained) {
			// remove it
			itH= spHabt.erase(itH);			
		}
		else {
			++itH;
		}
	}
}
void Simulation::addSpeciesIfcolonizing(map<vector<double>, Matrix<float> >& spInfoCol, map<vector<double>, Matrix<float> >& spAdult,
				map<vector<double>, Matrix<float> >& spRecruit, map<vector<double>, Matrix<float> >& spSeed, 
	           Matrix<float>& K,map<vector<double>, Matrix<float> >& spDispCol, 
			   map<vector<double>, Matrix<float> >& spDisp,
			   map<vector<double>, Matrix<float> >& spHabtCol,map<vector<double>, Matrix<float> >& spHabt)
{
		for(map<vector<double>, Matrix<float> >::iterator itr = spInfoCol.begin(); itr != spInfoCol.end(); ++itr) {
			const vector<double>& parameters = itr->first; // key (== parameterSet)
			Matrix<float>& abundances = itr->second; // value (== abundances)
			bool contained = spSeed.find(parameters) != spSeed.end();
			if(!contained) {//create maps
				spSeed[parameters] = abundances;
				spRecruit[parameters] = Matrix<float>(NY,NY);
				spAdult[parameters] = Matrix<float>(NY,NY);
				spDisp[parameters] = spDispCol[parameters];
				spHabt[parameters] = spHabtCol[parameters];

				if (anaspeciation=="on"){
					//for the gene flow map affecting anagenesis
					double TimestepforAnagenesis= timeStep +10000;
					TimestepforAnagenesis= timeStep + (pow(parameters[16],0.25)*exp(E/(Boltzmann*(temperature-parameters[19])))*allometabolConstTree[12]);
					immGeneflow[parameters] = Matrix<float>(2,2);
					immGeneflow[parameters] (0,0)=TimestepforAnagenesis;
					immGeneflow[parameters] (0,1)=0;
				}
			}
			else{ //increment seedbank
                 for (unsigned int i=0; i<patchesPosition.size(); i++){
					 spSeed[parameters] (patchesPosition[i]) += abundances(patchesPosition[i]);
					 if (anaspeciation=="on")  immGeneflow[parameters](0,1) += abundances(patchesPosition[i]);
				 }
			}
		}
}

void Simulation::localdynamics(map<vector<double>, Matrix<float> >& spAdult, map<vector<double>, Matrix<float> >& spRecruit,
					map<vector<double>, Matrix<float> >& spSeed, Matrix<float>& K,
	               map<vector<double>, Matrix<float> >& spDisp, map<vector<double>, Matrix<float> >& spHabt, map<vector<double>, Matrix<float> >& spHabtMain)
{
	Matrix<float> AreaCount (NY,NY);

	//Calculating tree cover for temperature feedback:
	if (treeTemperatureFeedback=="with") {		
		TempKelvin=TempKelvinBareIsland;
		//cout <<"TempKelvin(centralx,centraly) equal to bare temperatures: "<< TempKelvin(centralx,centraly) << endl;
		float totalTrees=0;
		for(map<vector<double>, Matrix<float> >::iterator itera = spAdult.begin(); itera != spAdult.end(); ++itera) 
		{
			const vector<double>& parameters = itera->first; // key (== parameterSet)
			Matrix<float>& abundances = itera->second; 
			
			if(parameters[16]>100000)
			{
				for (unsigned int i=0; i<patchesPosition.size(); i++ )
				{			
					totalTrees+= spAdult[parameters](patchesPosition[i])/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*0.0000199);
				}
			}

		}
		 double treeCover=totalTrees/(maxK*patchesPosition.size());
		 if(treeCover>1) treeCover=1;
		 if(treeCover<0) treeCover=0;
		 for (int ig=0; ig<(int)NY; ig++ )
		 {
			for (int ih=0; ih<(int)NY; ih++ )
			{
				if(State(ig,ih)>0)
				{
					TempKelvin(ig,ih) = TempKelvin(ig,ih) - 5*treeCover;
				}
			}
		}
    
	}
	//upgrading states of the population:
	//cout <<"TempKelvin(centralx,centraly) after tree cover update: "<< TempKelvin(centralx,centraly) << endl;
	for(map<vector<double>, Matrix<float> >::iterator itr = spSeed.begin(); itr != spSeed.end(); ++itr) 
	{
		
		const vector<double>& parameters = itr->first; // key (== parameterSet)
		Matrix<float>& abundances = itr->second; // value (== abundances)
		//update habitat suitability after tree cover change:
		vector<double> parameterSet = parameters;
		if (treeTemperatureFeedback=="with"){
			habitatSuitability(parameterSet,spHabt);
			if(parameters[4]==mainlandSpeciesNumber*2) habitatSuitability(parameterSet,spHabtMain);
		}

		//For mainland-island differentiation:
		
		float suminitialSeeds =0;
		float sumrecruitingSeeds =0;
		float suminitialRecruits =0;
		float summaturingRecruits =0;

		for (unsigned int i=0; i<patchesPosition.size(); i++ )
		{			
			//for mainland-island differentiation:
			if (anaspeciation=="on"){
				suminitialSeeds+=spSeed[parameters](patchesPosition[i]);
				suminitialRecruits+=spRecruit[parameters](patchesPosition[i]);
			}

			//Updating first the Adults (stronger competitors):
			//Growth to adulthood:
			int adults;			                                              
			if (metabolicRateToCorrect=="none") adults = Zufall.binomial(1 - exp(-(1/pow(parameters[20],0.25))*exp(-(E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[8]),(int)spRecruit[parameters](patchesPosition[i]));
			if ((metabolicRateToCorrect=="bgr")||(metabolicRateToCorrect=="all")) adults = Zufall.binomial(1 - exp(-parameters[11]),(int)spRecruit[parameters](patchesPosition[i]));

			float currentAdults=spAdult[parameters](patchesPosition[i])+adults;

			//for mainland-island differentiation:
			if (anaspeciation=="on") summaturingRecruits+= adults;

			//For updating K:
			float areaRecruitafterGrowth,areaTotalbeforeGrowth, correctedK;
			  if (metabolicRateToCorrect=="none"){
				areaRecruitafterGrowth=(spRecruit[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13]))) - (adults*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13])));
				areaTotalbeforeGrowth= (spRecruit[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13])))+spAdult[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13]));	
				correctedK = (K(patchesPosition[i])*spHabt[parameters](patchesPosition[i]))/(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13]));
			 }
			  if ((metabolicRateToCorrect=="all")||(metabolicRateToCorrect=="brmax")){
				areaRecruitafterGrowth=(spRecruit[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13]))) - (adults*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13])));
				areaTotalbeforeGrowth= (spRecruit[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13])))+spAdult[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13]));	
				correctedK = (K(patchesPosition[i])*spHabt[parameters](patchesPosition[i]))/(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13]));
			 }
			 if(correctedK>0)                          //
				{                             
					if(currentAdults<=correctedK) spAdult[parameters](patchesPosition[i]) = currentAdults;
					if(currentAdults>correctedK)  spAdult[parameters](patchesPosition[i]) = correctedK;
			}

			//For updating K:
			float areaTotalafterGrowth; 
			if (metabolicRateToCorrect=="none")areaTotalafterGrowth= areaRecruitafterGrowth + spAdult[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13]));
			if ((metabolicRateToCorrect=="all")||(metabolicRateToCorrect=="brmax"))areaTotalafterGrowth= areaRecruitafterGrowth + spAdult[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13]));
			
			AreaCount(patchesPosition[i])+=areaTotalafterGrowth;
			//updating K
			if(AreaCount(patchesPosition[i])>maxK) K(patchesPosition[i])=0;
			else{
				K(patchesPosition[i])=K(patchesPosition[i]) + (areaTotalbeforeGrowth-areaTotalafterGrowth);
			}

			if(K(patchesPosition[i])>maxK) K(patchesPosition[i])=maxK;
			if(K(patchesPosition[i])<0) K(patchesPosition[i])=0;
			
			//Seedling mortality after growth for those who did not grow:
			float survivingSeedlings;
			if (metabolicRateToCorrect=="none")   survivingSeedlings = (float) Zufall.binomial(1-(1 - exp(-((1/pow(parameters[20],0.25))*exp(-(E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[9]))),(int)spRecruit[parameters](patchesPosition[i]));
			if ((metabolicRateToCorrect=="bm")||(metabolicRateToCorrect=="all")||(metabolicRateToCorrect=="brmax&bm")) survivingSeedlings = (float) Zufall.binomial(1-(1 - exp(-parameters[12])),(int)spRecruit[parameters](patchesPosition[i]));
		
			float areaDeadRecruits,areaRecruitsbeforeRecruitment;
			if (metabolicRateToCorrect=="none")  areaDeadRecruits = (spRecruit[parameters](patchesPosition[i])-survivingSeedlings)*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13]));
			if ((metabolicRateToCorrect=="all")||(metabolicRateToCorrect=="brmax"))  areaDeadRecruits = (spRecruit[parameters](patchesPosition[i])-survivingSeedlings)*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13]));
			spRecruit[parameters](patchesPosition[i]) = survivingSeedlings;
			if (metabolicRateToCorrect=="none")  areaRecruitsbeforeRecruitment = spRecruit[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13]));
			if ((metabolicRateToCorrect=="all")||(metabolicRateToCorrect=="brmax"))  areaRecruitsbeforeRecruitment = spRecruit[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13]));
			//for updating K again after mortality and before recruitment:
            AreaCount(patchesPosition[i])-=areaDeadRecruits;		

			//updating K
			if(AreaCount(patchesPosition[i])>maxK) K(patchesPosition[i])=0;
			else{
				K(patchesPosition[i])=K(patchesPosition[i]) + areaDeadRecruits;
			}
			if(K(patchesPosition[i])>maxK) K(patchesPosition[i])=maxK;
			if(K(patchesPosition[i])<0) K(patchesPosition[i])=0;

			//Germination:
			int recruits;
			if (metabolicRateToCorrect=="none")  recruits = Zufall.binomial(1 - exp(-(1/pow(parameters[21],0.25))*exp(-(E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[10]),(int)abundances(patchesPosition[i])); 
			if ((metabolicRateToCorrect=="bge")||(metabolicRateToCorrect=="all"))  recruits = Zufall.binomial(1 - exp(-parameters[13]),(int)abundances(patchesPosition[i])); 
			
			//updating seedbank:
			abundances(patchesPosition[i]) -= recruits;
			//Seed mortality:
			if (metabolicRateToCorrect=="none")  abundances(patchesPosition[i]) = (float) Zufall.binomial(1-(1 - exp(-((1/pow(parameters[21],0.25))*exp(-(E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[11]))),((int)abundances(patchesPosition[i])));   
			if ((metabolicRateToCorrect=="bm")||(metabolicRateToCorrect=="all")||(metabolicRateToCorrect=="brmax&bm")) abundances(patchesPosition[i]) = (float) Zufall.binomial(1-(1 - exp(-parameters[14])),((int)abundances(patchesPosition[i])));    
			
			//for anagenesis:
			if (anaspeciation=="on") sumrecruitingSeeds+= recruits;

			// Then updating Recruits:	

			float correctedK2; 
			if (metabolicRateToCorrect=="none")	correctedK2= (K(patchesPosition[i])*spHabt[parameters](patchesPosition[i]))/(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13]));///(float)parameters[22]
			if ((metabolicRateToCorrect=="all")||(metabolicRateToCorrect=="brmax"))	correctedK2= (K(patchesPosition[i])*spHabt[parameters](patchesPosition[i]))/(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13]));///(float)parameters[22]
			if(correctedK2>0)
						{
							if(recruits<=correctedK2) spRecruit[parameters](patchesPosition[i]) += recruits;
							if(recruits>correctedK2)  spRecruit[parameters](patchesPosition[i]) = correctedK2;
					}
			float areaRecruitsafterRecruitment; 
			if (metabolicRateToCorrect=="none")	areaRecruitsafterRecruitment = spRecruit[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13]));
			if ((metabolicRateToCorrect=="all")||(metabolicRateToCorrect=="brmax"))	areaRecruitsafterRecruitment = spRecruit[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13]));
			//for updating K again after recruitment:
            AreaCount(patchesPosition[i])+= areaRecruitsafterRecruitment-areaRecruitsbeforeRecruitment;		
			//updating K
			if(AreaCount(patchesPosition[i])>maxK) K(patchesPosition[i])=0;
			else{
				K(patchesPosition[i])=K(patchesPosition[i]) - (areaRecruitsafterRecruitment-areaRecruitsbeforeRecruitment);
			}
			if(K(patchesPosition[i])>maxK) K(patchesPosition[i])=maxK;
			if(K(patchesPosition[i])<0) K(patchesPosition[i])=0;

		}

		//Annuals: just loop for the recruits to update them to adults:
        if((parameters[16]<=1000)&&(parameters[15]==1)){
			for (unsigned int i=0; i<patchesPosition.size(); i++ )
			{
				float survivingSeedlings =0;
				//growth to adulthood for all those who survive:
				int adults;
				if (metabolicRateToCorrect=="none") adults = Zufall.binomial(1- (1 - exp(-(1/pow(parameters[20],0.25))*exp(-(E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[9])),(int)spRecruit[parameters](patchesPosition[i]));
				if ((metabolicRateToCorrect=="bgr")||(metabolicRateToCorrect=="all")) adults = Zufall.binomial(1- (1 - exp(-parameters[12])),(int)spRecruit[parameters](patchesPosition[i]));
				//For updating K:
				float areaRecruitTminus1=0;
				if (metabolicRateToCorrect=="none")  areaRecruitTminus1 = spRecruit[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13]));
				if ((metabolicRateToCorrect=="all")||(metabolicRateToCorrect=="brmax"))  areaRecruitTminus1 = spRecruit[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13]));
				float areaAdultTminus1; 
				if (metabolicRateToCorrect=="none") areaAdultTminus1=  spAdult[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13]));
				if ((metabolicRateToCorrect=="all")||(metabolicRateToCorrect=="brmax")) areaAdultTminus1=  spAdult[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13]));
				
				spRecruit[parameters](patchesPosition[i]) = 0;

				//for anagenesis:
			    if (anaspeciation=="on") summaturingRecruits+= adults;

				//And updating Adults:
				float currentAdults=spAdult[parameters](patchesPosition[i])+adults;
				
				float correctedK3; 
				if (metabolicRateToCorrect=="none") correctedK3= (K(patchesPosition[i])*spHabt[parameters](patchesPosition[i]))/(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13]));///(float)parameters[8]
				if ((metabolicRateToCorrect=="all")||(metabolicRateToCorrect=="brmax")) correctedK3= (K(patchesPosition[i])*spHabt[parameters](patchesPosition[i]))/(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13]));///(float)parameters[8]
				if(correctedK3>0)
				{
					if(currentAdults<=correctedK3) spAdult[parameters](patchesPosition[i]) = currentAdults;
					if(currentAdults>correctedK3)  spAdult[parameters](patchesPosition[i]) = correctedK3;
				}

				//AbundCount(patchesPosition[i])+= ((spAdult[parameters](patchesPosition[i])*(float) parameters[8])- abundAdultTminus1) -abundRecruitTminus1;
				if (metabolicRateToCorrect=="none") AreaCount(patchesPosition[i])+= ((spAdult[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13])))- areaAdultTminus1) -areaRecruitTminus1;
				if ((metabolicRateToCorrect=="all")||(metabolicRateToCorrect=="brmax")) AreaCount(patchesPosition[i])+= ((spAdult[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13])))- areaAdultTminus1) -areaRecruitTminus1;
				//updating K
				if(AreaCount(patchesPosition[i])>maxK) K(patchesPosition[i])=0;
				else{
					if (metabolicRateToCorrect=="none")K(patchesPosition[i])=K(patchesPosition[i]) - (((spAdult[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13])))- areaAdultTminus1) -areaRecruitTminus1);
					if ((metabolicRateToCorrect=="all")||(metabolicRateToCorrect=="brmax"))K(patchesPosition[i])=K(patchesPosition[i]) - (((spAdult[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13])))- areaAdultTminus1) -areaRecruitTminus1);
				}

				if(K(patchesPosition[i])>maxK) K(patchesPosition[i])=maxK;
				if(K(patchesPosition[i])<0) K(patchesPosition[i])=0;
				//cout <<"K After updating : "<< K(patchesPosition[i])<< endl;
			}
		}
		
		if (anaspeciation=="on"){
			bool contained = immGeneflow.find(parameters) != immGeneflow.end();
			if((contained)&&(immGeneflow[parameters](0,1)>=1)) {
				//updating gene flow map for anagenesis:   suminitialSeeds,sumrecruitingSeeds,suminitialRecruits,summaturingRecruits
					double immigrantrecruiting = (immGeneflow[parameters](0,1)/suminitialSeeds)*sumrecruitingSeeds;
					double immigrantmaturing = (immigrantrecruiting/suminitialRecruits)*summaturingRecruits;
					double GenerationTime = (pow(parameters[16],0.25)*exp(E/(Boltzmann*(temperature-parameters[19])))*allometabolConstTree[12])*0.0001;
					immGeneflow[parameters](0,0) =immGeneflow[parameters](0,0) +10*GenerationTime*immigrantmaturing;
					double maxanatime = timeStep + (pow(parameters[16],0.25)*exp(E/(Boltzmann*(temperature-parameters[19])))*allometabolConstTree[12]);
					if (immGeneflow[parameters](0,0)>maxanatime) immGeneflow[parameters](0,0) =maxanatime; //time for anagenesis cannot exceed the number of generations for it
					immGeneflow[parameters](0,1) =0;
			}
		}
	}	
	
	//updating Carrying capacity
	Matrix<float> areaCount (NY,NY);
	for(map<vector<double>, Matrix<float> >::iterator itr = spSeed.begin(); itr != spSeed.end(); ++itr) 
	{
		const vector<double>& parameters = itr->first; // key (== parameterSet)
		Matrix<float>& abundances = itr->second; 
		for(unsigned int i=0; i<patchesPosition.size(); i++) 
	    {
			int c = patchesPosition[i].first;
			int r = patchesPosition[i].second;		
			if(spAdult[parameters](c,r)<0) spAdult[parameters](c,r)=0;	
			if(spRecruit[parameters](c,r)<0) spRecruit[parameters](c,r)=0;	
			if(spSeed[parameters](c,r)<0) spSeed[parameters](c,r)=0;	
			float AreaTotal; 
			if (metabolicRateToCorrect=="none")AreaTotal =spAdult[parameters](c,r)*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(c,r))))*allometabolConstTree[13])) +spRecruit[parameters](c,r)*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*TempKelvin(c,r))))*allometabolConstTree[13]));
			if ((metabolicRateToCorrect=="all")||(metabolicRateToCorrect=="brmax"))AreaTotal =spAdult[parameters](c,r)*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13])) +spRecruit[parameters](c,r)*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13]));
			areaCount(c,r)+=AreaTotal;
		}
	}
	for(unsigned int i=0; i<patchesPosition.size(); i++) 
    {
		int c = patchesPosition[i].first;
		int r = patchesPosition[i].second;		
		K(c,r) = maxK-areaCount(c,r);
		if(K(c,r)<0) K(c,r) =0;	
	}
	//Mutant maps:
	speciesAdultMutant.clear();
	speciesRecruitMutant.clear();
	speciesSeedMutant.clear();
	speciesDispMutant.clear();
    speciesHabtMutant.clear();
	
	//reproduction and speciation:
	for(map<vector<double>, Matrix<float> >::iterator itr = spAdult.begin(); itr != spAdult.end(); ++itr) {
		const vector<double>& parameters = itr->first; // key (== parameterSet)
		Matrix<float>& abundances = itr->second; // value (== abundances)
		 
		//Matrix<float> incomingDispersalUnits (NY,NY);  
		 
		//Using matrix addition of library armadillo (a lot faster):	
		mat D(NY*2 +1,NY*2 +1);
		D.fill(0.0);
		mat A(NY,NY);
		A.fill(0.0);
		for (int co=0; co<(int)NY*2 +1; co++ )
		 {
		 	for (int ro=0; ro<(int)NY*2+1; ro++ )
		 	{						
				D(co,ro) =spDisp[parameters]( abs((int)NY +1 - co) , abs((int)NY+1 - ro) );					
			}
		 } 
		
		//D.print("D:");
		/*cout << "central lines:" << endl;
		cout << D(floor(NY*0.5)-1,floor(NY*0.5)-1)<<"  "<<D(floor(NY*0.5),floor(NY*0.5)-1)<<"  "<<D(floor(NY*0.5)+1,floor(NY*0.5)-1)<<"  "<<D(floor(NY*0.5)+2,floor(NY*0.5)-1)<<"  "<<D(floor(NY*0.5)+3,floor(NY*0.5)-1)<<endl;
		cout << D(floor(NY*0.5)-1,floor(NY*0.5))<<"  "<<D(floor(NY*0.5),floor(NY*0.5))<<"  "<<D(floor(NY*0.5)+1,floor(NY*0.5))<<"  "<<D(floor(NY*0.5)+2,floor(NY*0.5))<<"  "<<D(floor(NY*0.5)+3,floor(NY*0.5))<<endl;
		cout << D(floor(NY*0.5)-1,floor(NY*0.5)+1)<<"  "<<D(floor(NY*0.5),floor(NY*0.5)+1)<<"  "<<D(floor(NY*0.5)+1,floor(NY*0.5)+1)<<"  "<<D(floor(NY*0.5)+2,floor(NY*0.5)+1)<<"  "<<D(floor(NY*0.5)+3,floor(NY*0.5)+1)<<endl;
        cout << D(floor(NY*0.5)-1,floor(NY*0.5)+2)<<"  "<<D(floor(NY*0.5),floor(NY*0.5)+2)<<"  "<<D(floor(NY*0.5)+1,floor(NY*0.5)+2)<<"  "<<D(floor(NY*0.5)+2,floor(NY*0.5)+2)<<"  "<<D(floor(NY*0.5)+3,floor(NY*0.5)+2)<<endl;
		cout << D(floor(NY*0.5)-1,floor(NY*0.5)+3)<<"  "<<D(floor(NY*0.5),floor(NY*0.5)+3)<<"  "<<D(floor(NY*0.5)+1,floor(NY*0.5)+3)<<"  "<<D(floor(NY*0.5)+2,floor(NY*0.5)+3)<<"  "<<D(floor(NY*0.5)+3,floor(NY*0.5)+3)<<endl;
		*/
         
		for (unsigned int i=0; i<patchesPosition.size(); i++ )
		{
    		// Beverton-Holt Modell extende for Allee effects
			float spTotalExploitedArea; 
			double Ki;
			if (metabolicRateToCorrect=="none") {
				spTotalExploitedArea= (spRecruit[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13])))+spAdult[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13]));	
				Ki= ((K(patchesPosition[i])+spTotalExploitedArea)*spHabt[parameters](patchesPosition[i]))/(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13]));
			}
			if ((metabolicRateToCorrect=="all")||(metabolicRateToCorrect=="brmax")){
				spTotalExploitedArea= (spRecruit[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13])))+spAdult[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13]));	
				Ki= ((K(patchesPosition[i])+spTotalExploitedArea)*spHabt[parameters](patchesPosition[i]))/(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13]));
			}
			/*cout <<"Timestep: "<< timeStep<< endl;
			cout <<"X: "<< patchesPosition[i].first<< endl;
			cout <<"Y: "<< patchesPosition[i].second<< endl;
			cout <<"Mass: "<< parameters[16]<< endl;
			cout <<"Adult area: "<< (maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13]))<< endl;
			cout <<"spAdult[parameters](patchesPosition[i]): "<< spAdult[parameters](patchesPosition[i])<< endl;
			cout <<"spTotalExploitedArea: "<< spTotalExploitedArea<< endl;
            cout <<"spadultsExploitedArea /: "<< spAdult[parameters](patchesPosition[i])*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13]))<< endl;
			cout <<"K(patchesPosition[i]): "<< K(patchesPosition[i])<< endl;
			cout <<"spHabt[parameters](patchesPosition[i]): "<< spHabt[parameters](patchesPosition[i])<< endl;
            cout <<"Ki: "<< Ki<< endl;*/

			double C = 0;
			if (metabolicRateToCorrect=="none") {
				if(parameters[3]<0) C= parameters[3]*(((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13])*spHabt[parameters](patchesPosition[i]));
				else C= parameters[3];
			}
			if ((metabolicRateToCorrect=="all")||(metabolicRateToCorrect=="brmax")){
				if(parameters[3]<0) C= parameters[3]*(((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13])*spHabt[parameters](patchesPosition[i]));
				else C= parameters[3];
			}
			//making minimum value of C being -K
			if(C<-Ki) C=-Ki;
			if (C==Ki) C= C-1;
			//Beverton-Holt carrying capacity for nonsprouters:
			double R, M, k, c;
			if (metabolicRateToCorrect=="none") {
				R = (1/pow(parameters[16],0.25))*exp(-(E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[1];
				M = (1/pow(parameters[16],0.25))*exp(-(E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[2];
			}
			if (metabolicRateToCorrect=="brmax") {
				R = parameters[1];
				M = (1/pow(parameters[16],0.25))*exp(-(E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[2];
			}
			if (metabolicRateToCorrect=="bm") {
				R = (1/pow(parameters[16],0.25))*exp(-(E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[1];
				M = parameters[2];
			}
			if ((metabolicRateToCorrect=="all")||(metabolicRateToCorrect=="brmax&bm")) {
				R = parameters[1];
				M = parameters[2];
			}
			k = 4*(R-M)/(M*pow((Ki-C),2));
			c = C + sqrt((R-M)/(M*k));	
			
			   			
			//Beverton-Holt:
			double dispersalUnits; 
			if (treeline =="on"){	
				if(R>0.15) 
					dispersalUnits= (abundances(patchesPosition[i])*R)/(1+ k*(pow((abundances(patchesPosition[i]) - c),2)));
				else {
					dispersalUnits= 0;
					if (timeStep<=1) abundances(patchesPosition[i])=0;

				}
			}
			else dispersalUnits= (abundances(patchesPosition[i])*R)/(1+ k*(pow((abundances(patchesPosition[i]) - c),2)));
			//Ricker:
			//dispersalUnits = abundances(patchesPosition[i])*pow(R,(1-(abundances(patchesPosition[i])/K(patchesPosition[i])))*((abundances(patchesPosition[i])-C)/K(patchesPosition[i]))*(4/(1-(C/K(patchesPosition[i]))*(1-(C/K(patchesPosition[i]))))));
			//cout <<"dispersalUnits: "<< dispersalUnits<< endl;

			int is = patchesPosition[i].first;
			int iz = patchesPosition[i].second;			
			//Use a dispersal kernel to disperse the dispersal units (use a dispersal funtion)           
			if(dispersalUnits>0){		
				
				//Using matrix addition of library armadillo (a lot faster):	
				A=A+  D.submat(NY+1-is,NY+1-iz,NY+(NY-is),NY+(NY-iz))*dispersalUnits;
				//A.print("A:");

				//Using loop over sink cells during dispersal(slower):
				//Absorbing boundaries:
				/*if (boundaries=="absorbing"){
					for (int ib=0; ib<(int)NY; ib++ )
					{
						for (int ia=0; ia<(int)NY; ia++ )
						{					
							incomingDispersalUnits(ib,ia) += (float) (spDisp[parameters](abs(ib-is),abs(ia-iz))*dispersalUnits);
						}
					}
				}
				//periodic boundaries:
				if (boundaries=="periodic"){
					for (int ib=(int)(is-maxIslandRadius); ib<(int)(is+maxIslandRadius); ib++ )
					{
						for (int ia=(int)(iz-maxIslandRadius); ia<(int)(iz+maxIslandRadius); ia++ )
						{
							int id, ie;
							id=ib;
							ie=ia;
							if(ib<1) id=NY-1+ib;
							if(ib>((int)NY-1))	id=ib-((int)NY-1);
							if(ia<1) ie=(int)NY-1+ia;
							if(ia>((int)NY-1)) ie=ia-((int)NY-1);
							incomingDispersalUnits(id,ie) += (float) (spDisp[parameters](abs(ib-is),abs(ia-iz))*dispersalUnits);
						}
					}
				}*/				
			}			
			
			 /*if((parameters[16]>1000000)&&(parameters[19]>22)&&(parameters[24]<4)&&(timeStep>=150)&&(patchesPosition[i].second<=5)){
					cout <<" "<< endl;
					cout <<"Timestep: "<< timeStep<< endl;
					cout <<"TempKelvin(patchesPosition[i]): "<< TempKelvin(patchesPosition[i])<< endl;
					//cout <<"State(patchesPosition[i]): "<< State(patchesPosition[i])<< endl;
					cout <<"optimal altitude : "<<  parameters[19]<< endl;
					cout <<"altitude amplitude: "<<  parameters[24]<< endl;
					cout <<"Mass: "<<  parameters[16]<< endl;
					cout <<"R: "<< R<< endl;
					cout <<"M: "<< M<< endl;
					//cout <<"Ki: "<< Ki<< endl;
					//cout <<"spTotalExploitedArea: "<< spTotalExploitedArea<< endl;
					//cout <<"spHabt[parameters](patchesPosition[i]): "<< spHabt[parameters](patchesPosition[i])<< endl;
					//cout <<"K(patchesPosition[i]): "<< K(patchesPosition[i])<< endl;
					cout <<"abundances(patchesPosition[i]): "<< abundances(patchesPosition[i])<< endl;
					cout <<"dispersalUnits: "<< dispersalUnits<< endl;
					cout <<"incomingDispersalUnits: "<< incomingDispersalUnits(patchesPosition[i])<< endl;
			 }*/
		}
		

		//use the incoming dispersal units to draw a random recruitment
		for(unsigned int i=0; i<patchesPosition.size(); i++)
		{
			int is = patchesPosition[i].first;
			int iz = patchesPosition[i].second;

			float mutants = 0;
            if (cladospeciation=="on"){
			//Speciation: random recruitment of mutants
				
				//Using matrix addition of library armadillo:					
				if (metabolicRateToCorrect=="all") mutants = (float) Zufall.poisson(A(is,iz)*(1 - exp(-parameters[5])));// OBS: below is the mutantion rate varying on the local sink cell, not on the souce cell!!
				else mutants = (float) Zufall.poisson(A(is,iz)*(1 - exp(-(1/pow(parameters[16],0.25))*exp(-(E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[5])));
				
				//Using loop over sink cells during dispersal(slower):
				//if (metabolicRateToCorrect=="all") mutants = (float) Zufall.poisson(incomingDispersalUnits(patchesPosition[i])*parameters[5]);
				//else mutants = (float) Zufall.poisson(incomingDispersalUnits(patchesPosition[i])*(1 - exp(-(1/pow(parameters[16],0.25))*exp(-(E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[5])));
				
				//if at least one mutant of this species recruited (only one mutant per timestep and island cell allowed
				//- all the mutants from the "same" diverging species):
				if(mutants>0)
				{
					vector<double> parameterSetMutant;
						//choosing random parameter values for hte new mutant +- half of the value
						double Annual = 0;
						//Mass:
						double max = (parameters[16]+ parameters[16]*nicheEvolution)*10;
						double min = (parameters[16]- parameters[16]*nicheEvolution)*10;
						double Mass = rand()%( (int)(max-min));
						//cout << "min Mass: "<<min<<"max Mass: " <<max<<"random Mass: " <<Mass<<endl;
						Mass = min + Mass;
						Mass = Mass*0.1;
						if(Mass<10) Mass = 10; //lower bound for herbs
					
						//cout << "mother Mass: "<<parameters[16]<<"corrected mutant Mass: " <<Mass<<endl;
						double MassSeedling, MassSeed;
						//Relative Fenology:
						max = (parameters[0]+ parameters[0]*nicheEvolution)*10;
						min = (parameters[0]- parameters[0]*nicheEvolution)*10;
						double Fenology = rand()%( (int)(max-min));
						Fenology = min + Fenology;
						Fenology = Fenology*0.1;

						//For hsuitability and metabolic rates:
						double optimalAltitude = (rand()%(maxOptaltitude-minOptaltitude)) + minOptaltitude;

						//Rmax:
						double Rmax= (1/pow(Mass,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[1];
						//M:
						double M = (1/pow(Mass,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[2];
						//Gamma:
						double Gamma = 0.001;
						float Evolutionrate = nicheEvolution*100;
						double change =((rand()%((int) Evolutionrate)) +1)*0.01;
						int direction = rand()%2;
						if(direction==0) Gamma = parameters[3] -parameters[3]*change;
						else Gamma = parameters[3] +parameters[3]*change;
				
						//GenTime or invasive:
						double GenTime =0;
						if(parameters[4]<mainlandSpeciesNumber*2) GenTime = rand()% mainlandSpeciesNumber*nicheEvolution;
						else GenTime = mainlandSpeciesNumber*2 + rand()%( 30000) +1;
						//MutRate:
						double MutRate= (1/pow(Mass,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[5];
						
						//Dispersal Ability (mean dispersal distance):
						max = (parameters[6]+ parameters[6]*nicheEvolution) *100000;
						min = (parameters[6]- parameters[6]*nicheEvolution) *100000;
						if(min<1) min=1;
						double Alpha = rand()%( (int)(max-min)+1);
						Alpha = min + Alpha;
						Alpha = Alpha*0.00001;
						if (Alpha<minValueTree[6]) Alpha=minValueTree[6];
						if (Alpha>maxValueTree[6]) Alpha=maxValueTree[6]; 
						//Dispersal Ability (P):

						max = (parameters[7]+ parameters[7]*nicheEvolution) *1000;
						min = (parameters[7]- parameters[7]*nicheEvolution) *1000;
						double P = rand()%( (int)(max-min));
						//cout <<"P: "<< P<< endl;
						//cout <<"nicheEvolution: "<< nicheEvolution<< endl;
						P = min + P;
						//cout <<"min + P: "<< P<< endl;
						//cout <<"P*0.001: "<< P*0.001<< endl;
						P = P*0.001;
						if(P<minValueTree[7]) P=minValueTree[7];
						if(P>maxValueTree[7]) P=maxValueTree[7];
						//K PFTFACTOR:
								//From Niklas (1994) Plant Allometry:
						double Kpftfactor,KpftfactorSeedling;
						if(parameters[8]>0.01) //allometabolConst[0]==1 => tree
						{  
						 double diameter = (1/63.0)*pow(Mass,0.4);
						 Kpftfactor = (pow(diameter,2))*0.01;
						 KpftfactorSeedling = (pow(diameter*0.5,2))*0.01;
						}
						else
						{
						 double diameter = (1/46.41589)*pow(Mass,0.3333);
						 Kpftfactor = (pow(diameter,2))*0.01;
						 KpftfactorSeedling = (pow(diameter*0.5,2))*0.01;
						}										
						double randomTrend = (rand()% 101)*0.01;
						double islandSide= (rand()% 4) +1;
					
						if(Kpftfactor<0.00045) Annual = rand()% 2;	
					
					//with the same metabolic constant than adults, but lower body mass:
					//Recruit M:
						double Mrecruit = 0;
						//Seed M:
						double Mseed =0;
						//Growth:
						double Growth =0;
						//Germ:
						double Germ = 0;
						if (Mass>55000)  {
							MassSeedling = Mass*0.2;
							MassSeed = (exp((rand()% 80)*0.1))*0.01; 
							Mrecruit =  (1/pow(MassSeedling,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[9];
							Mseed =  (1/pow(MassSeed,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[11];//tree
							Growth =  (1/pow(MassSeedling,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[8];
							Germ =  (1/pow(MassSeed,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[10];
						}
						else {
							MassSeedling = Mass*0.2;
							MassSeed= pow(Mass*0.0000046,0.5);//0.00043*Mass;
							Mrecruit =  (1/pow(MassSeedling,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[9];
							Mseed =  (1/pow(MassSeed,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[11];//shrub
							Growth =  (1/pow(MassSeedling,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[8];
							Germ =  (1/pow(MassSeed,0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[10];
						}

						//specialist: 1; generalist: 4:
						double generalist = (rand()% 4) +1;
						//vertical specialist: 1; generalist: verticalpreference:
						double verticalgeneralist = (rand()% verticalpreference) +1;
						
						//metabolic year at which the pre-cladospecies will become a fully recognized new cladospecies
						double yearOfSpec = timeStep + (pow(Mass,0.25)*exp(E/(Boltzmann*tempForFixedMetabolism))*allometabolConstTree[12]);
						double SpType=1;
						if(parameters[25]==4) SpType =5;

						parameterSetMutant.push_back(Fenology);
  						parameterSetMutant.push_back(Rmax);
						parameterSetMutant.push_back(M);
						parameterSetMutant.push_back(Gamma);
						parameterSetMutant.push_back(GenTime);
						parameterSetMutant.push_back(MutRate);
						parameterSetMutant.push_back(Alpha);
						parameterSetMutant.push_back(P);
						parameterSetMutant.push_back(Kpftfactor);
						parameterSetMutant.push_back(randomTrend);
						parameterSetMutant.push_back(islandSide);
						parameterSetMutant.push_back(Growth);
						parameterSetMutant.push_back(Mrecruit);
						parameterSetMutant.push_back(Germ);
						parameterSetMutant.push_back(Mseed);
						parameterSetMutant.push_back(Annual);
						parameterSetMutant.push_back(Mass);
						parameterSetMutant.push_back(parameters[17]);//lineage
						parameterSetMutant.push_back(parameters[17]*1000+parameters[16]); //MotherID 
						parameterSetMutant.push_back(optimalAltitude);
						parameterSetMutant.push_back(MassSeedling);
						parameterSetMutant.push_back(MassSeed);
						parameterSetMutant.push_back(KpftfactorSeedling);
						parameterSetMutant.push_back(generalist);
						parameterSetMutant.push_back(verticalgeneralist);
						parameterSetMutant.push_back(SpType); //1= pre-cladospecies from the initial source pool; 5 =pre-cladospecies from an invasive
						parameterSetMutant.push_back(yearOfSpec); //time step of mutation
					

				     if(parameterSetMutant[16]<=15.0){
							mutants=0;
					}
					else
					{
						addRandomSpeciesIfNew(parameterSetMutant,speciesSeedMutant, speciesDispMutant, speciesHabtMutant);
						speciesSeedMutant[parameterSetMutant](patchesPosition[i]) = mutants;
					 }			

				}

			}
		    //Recruitment:
			float seeds;
			if (treeline =="on"){
				double Gr; 
				if (metabolicRateToCorrect=="none") Gr= (1/pow(parameters[20],0.25))*exp(-(E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[8];
				if ((metabolicRateToCorrect=="bgr")||(metabolicRateToCorrect=="all")) Gr= (1/pow(parameters[20],0.25))*exp(-(E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[8];
				if(Gr>=0.012){
					//Using matrix addition of library armadillo:
					seeds = (float)  Zufall.poisson(A(is,iz)-mutants); 
					//Using loop over sink cells during dispersal(slower):
					//seeds = (float)  Zufall.poisson(incomingDispersalUnits(patchesPosition[i])-mutants);
				}
				else{
					seeds=0;
					spSeed[parameters](patchesPosition[i])=0;
					spRecruit[parameters](patchesPosition[i])=0;
					spAdult[parameters](patchesPosition[i])=0;
				}
			}

			//Using matrix addition of library armadillo:
			else seeds = (float)  Zufall.poisson(A(is,iz)-mutants);
			//Using loop over sink cells during dispersal (slower):
			//else seeds = (float)  Zufall.poisson(incomingDispersalUnits(patchesPosition[i])-mutants);

			//updating seedbank:			
			if(vulnerableSpecialists=="no"){
				if((parameters[24]>=1)&&(spHabt[parameters](patchesPosition[i])>0)) spSeed[parameters](patchesPosition[i])=spSeed[parameters](patchesPosition[i])+seeds;
				if(parameters[24]<=1) spSeed[parameters](patchesPosition[i])=spSeed[parameters](patchesPosition[i])+seeds;
			}
			if(vulnerableSpecialists=="yes"){ //Attention: this makes altitude specialists vulnerable to extinction during island dynamics!!!!!!!!!!!!!!!
				if(spHabt[parameters](patchesPosition[i])>0) spSeed[parameters](patchesPosition[i])=spSeed[parameters](patchesPosition[i])+seeds;
			}

			//Adult mortality:
			float survivors;
			if (metabolicRateToCorrect=="none")survivors = (float) Zufall.binomial(1-(1 - exp(-((1/pow(parameters[16],0.25))*exp(-(E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[2]))),(int)abundances(patchesPosition[i]));
			if ((metabolicRateToCorrect=="bm")||(metabolicRateToCorrect=="all")||(metabolicRateToCorrect=="brmax&bm")) survivors = (float) Zufall.binomial(1-(1 - exp(-parameters[2])),(int)abundances(patchesPosition[i]));
			
			if((parameters[16]<1000)&&(parameters[15]==1))survivors = 0;
			
			//cout <<"survivor adults: "<< survivors<< endl;
			//update abundance and K:
			float areaTimeminus1 ;
			if (metabolicRateToCorrect=="none")areaTimeminus1=abundances(patchesPosition[i])*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13]));
			if ((metabolicRateToCorrect=="all")||(metabolicRateToCorrect=="brmax"))areaTimeminus1=abundances(patchesPosition[i])*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13]));
				
			//double KbeforeRecruit = K(patchesPosition[i]);
			abundances(patchesPosition[i])= survivors;
			
			if (metabolicRateToCorrect=="none") AreaCount(patchesPosition[i])-=(areaTimeminus1-(abundances(patchesPosition[i])*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13]))));
			if ((metabolicRateToCorrect=="all")||(metabolicRateToCorrect=="brmax")) AreaCount(patchesPosition[i])-=(areaTimeminus1-(abundances(patchesPosition[i])*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13]))));
			//updating K:
			if(AreaCount(patchesPosition[i])>maxK) K(patchesPosition[i])=0;
			else{
				if (metabolicRateToCorrect=="none")K(patchesPosition[i])=K(patchesPosition[i]) + (areaTimeminus1-(abundances(patchesPosition[i])*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(patchesPosition[i]))))*allometabolConstTree[13]))));			
				if ((metabolicRateToCorrect=="all")||(metabolicRateToCorrect=="brmax"))K(patchesPosition[i])=K(patchesPosition[i]) + (areaTimeminus1-(abundances(patchesPosition[i])*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*tempForFixedMetabolism)))*allometabolConstTree[13]))));			
			}
			if(K(patchesPosition[i])>maxK) K(patchesPosition[i])=maxK;
            if(K(patchesPosition[i])<0) K(patchesPosition[i])=0;

		}
			
	}
	
	
    addSpeciesIfcolonizing(speciesSeedMutant, spAdult, spRecruit,spSeed, K, speciesDispMutant, spDisp,speciesHabtMutant,spHabt);
	
	//updating species status (from pre-cladospecies to new evovled cladospecies): 1->2
	if (cladospeciation=="on"){
		SpecClado=0; 
		LookforCladoSpecies (spAdult,spRecruit,spSeed,spDisp,spHabt);
	}
	if (anaspeciation=="on"){
		SpecAna=0;
		LookforAnaSpecies (spAdult,spRecruit,spSeed,spDisp,spHabt,immGeneflow);
	}
}
void Simulation::localdynamicsNoMTE(map<vector<double>, Matrix<float> >& spAdult, map<vector<double>, Matrix<float> >& spRecruit,
					map<vector<double>, Matrix<float> >& spSeed, Matrix<float>& K,
	               map<vector<double>, Matrix<float> >& spDisp, map<vector<double>, Matrix<float> >& spHabt, map<vector<double>, Matrix<float> >& spHabtMain)
{
	Matrix<float> AreaCount (NY,NY);
    
	//upgrading states of the population:
	//cout <<"TempKelvin(centralx,centraly) after tree cover update: "<< TempKelvin(centralx,centraly) << endl;
	for(map<vector<double>, Matrix<float> >::iterator itr = spSeed.begin(); itr != spSeed.end(); ++itr) 
	{
		
		const vector<double>& parameters = itr->first; // key (== parameterSet)
		Matrix<float>& abundances = itr->second; // value (== abundances)
		//update habitat suitability after tree cover change:
		vector<double> parameterSet = parameters;
		
		//For anagenesis:
		
		float suminitialSeeds =0;
		float sumrecruitingSeeds =0;
		float suminitialRecruits =0;
		float summaturingRecruits =0;

		for (unsigned int i=0; i<patchesPosition.size(); i++ )
		{			
			//for mainland-island differentiation:
			if (anaspeciation=="on"){
				suminitialSeeds+=spSeed[parameters](patchesPosition[i]);
				suminitialRecruits+=spRecruit[parameters](patchesPosition[i]);
			}

			//Updating first the Adults (stronger competitors):
			//Growth to adulthood:
			int adults;
			adults = Zufall.binomial(1 - exp(-parameters[11]),(int)spRecruit[parameters](patchesPosition[i]));
			float currentAdults=spAdult[parameters](patchesPosition[i])+adults;
			

			//for mainland-island differentiation:
			if (anaspeciation=="on") summaturingRecruits+= adults;

			//For updating K: TO-DO: making a random non-metabolic kfactor for non-metabolic simus  0.0001125
			float areaRecruitafterGrowth,areaTotalbeforeGrowth, correctedK;
			areaRecruitafterGrowth=(spRecruit[parameters](patchesPosition[i])*((pow(parameters[20],0.75))*0.0001125)) - (adults*((pow(parameters[20],0.75))*0.0001125));
			areaTotalbeforeGrowth= (spRecruit[parameters](patchesPosition[i])*((pow(parameters[20],0.75))*0.0001125))+spAdult[parameters](patchesPosition[i])*((pow(parameters[16],0.75))*0.0001125);	
			correctedK = (K(patchesPosition[i])*spHabt[parameters](patchesPosition[i]))/((pow(parameters[16],0.75))*0.0001125);
			if(correctedK>0)                          //
			{                             
				if(currentAdults<=correctedK) spAdult[parameters](patchesPosition[i]) = currentAdults;
				if(currentAdults>correctedK)  spAdult[parameters](patchesPosition[i]) = correctedK;
			}

			//For updating K:
			float areaTotalafterGrowth; 
			areaTotalafterGrowth= areaRecruitafterGrowth + spAdult[parameters](patchesPosition[i])*((pow(parameters[16],0.75))*0.0001125);			
			AreaCount(patchesPosition[i])+=areaTotalafterGrowth;
			//updating K
			if(AreaCount(patchesPosition[i])>maxK) K(patchesPosition[i])=0;
			else K(patchesPosition[i])=K(patchesPosition[i]) + (areaTotalbeforeGrowth-areaTotalafterGrowth);
			if(K(patchesPosition[i])>maxK) K(patchesPosition[i])=maxK;
			if(K(patchesPosition[i])<0) K(patchesPosition[i])=0;
			
			//Seedling mortality after growth for those who did not grow:
			float survivingSeedlings;
			survivingSeedlings = (float) Zufall.binomial((1-(1 - exp(-parameters[12]))),(int)spRecruit[parameters](patchesPosition[i]));	
			
			float areaDeadRecruits,areaRecruitsbeforeRecruitment;
			areaDeadRecruits = (spRecruit[parameters](patchesPosition[i])-survivingSeedlings)*((pow(parameters[20],0.75))*0.0001125);
			spRecruit[parameters](patchesPosition[i]) = survivingSeedlings;
			areaRecruitsbeforeRecruitment = spRecruit[parameters](patchesPosition[i])*((pow(parameters[20],0.75))*0.0001125);
			//for updating K again after mortality and before recruitment:
            AreaCount(patchesPosition[i])-=areaDeadRecruits;		

			//updating K
			if(AreaCount(patchesPosition[i])>maxK) K(patchesPosition[i])=0;
			else K(patchesPosition[i])=K(patchesPosition[i]) + areaDeadRecruits;			
			if(K(patchesPosition[i])>maxK) K(patchesPosition[i])=maxK;
			//Germination:
			int recruits;
			recruits = Zufall.binomial(1 - exp(-parameters[13]),(int)abundances(patchesPosition[i])); 

			//updating seedbank:
			abundances(patchesPosition[i]) -= recruits;
			//Seed mortality:
			abundances(patchesPosition[i]) = (float) Zufall.binomial((1-(1 - exp(-parameters[14]))),((int)abundances(patchesPosition[i])));   
			//for anagenesis:
			if (anaspeciation=="on") sumrecruitingSeeds+= recruits;

			// Then updating Recruits:	
			float correctedK2; 
			correctedK2= K(patchesPosition[i])*spHabt[parameters](patchesPosition[i])/((pow(parameters[20],0.75))*0.0001125);
			if(correctedK2>0)
			{
				if(recruits<=correctedK2) spRecruit[parameters](patchesPosition[i]) += recruits;
				if(recruits>correctedK2)  spRecruit[parameters](patchesPosition[i]) = correctedK2;
			}

			float areaRecruitsafterRecruitment; 
			areaRecruitsafterRecruitment = spRecruit[parameters](patchesPosition[i])*((pow(parameters[20],0.75))*0.0001125);
			//for updating K again after recruitment:
            AreaCount(patchesPosition[i])+= areaRecruitsafterRecruitment-areaRecruitsbeforeRecruitment;		
			//updating K
			if(AreaCount(patchesPosition[i])>maxK) K(patchesPosition[i])=0;
			else K(patchesPosition[i])=K(patchesPosition[i]) + (areaRecruitsbeforeRecruitment-areaRecruitsafterRecruitment);			
			if(K(patchesPosition[i])>maxK) K(patchesPosition[i])=maxK;
			if(K(patchesPosition[i])<0) K(patchesPosition[i])=0;
		}

		//Annuals: just loop for the recruits to update them to adults:
        if((parameters[16]<=1000)&&(parameters[15]==1)){
			for (unsigned int i=0; i<patchesPosition.size(); i++ )
			{
				float survivingSeedlings =0;
				//growth to adulthood for all those who survive:
				int adults;
				adults = Zufall.binomial((1-(1 - exp(-parameters[12]))),(int)spRecruit[parameters](patchesPosition[i]));
				
				//For updating K:
				float areaRecruitTminus1=0;
				areaRecruitTminus1 = spRecruit[parameters](patchesPosition[i])*((pow(parameters[20],0.75))*0.0001125);
				float areaAdultTminus1; 
				areaAdultTminus1=  spAdult[parameters](patchesPosition[i])*((pow(parameters[16],0.75))*0.0001125);
				
				spRecruit[parameters](patchesPosition[i]) = 0;

				//for mainland-island differentiation:
			    if (anaspeciation=="on") summaturingRecruits+= adults;

				//And updating Adults:
				float currentAdults=spAdult[parameters](patchesPosition[i])+adults;
				
				float correctedK3; 
				correctedK3= (K(patchesPosition[i])*spHabt[parameters](patchesPosition[i]))/((pow(parameters[16],0.75))*0.0001125);
				if(correctedK3>0)
				{
					if(currentAdults<=correctedK3) spAdult[parameters](patchesPosition[i]) = currentAdults;
					if(currentAdults>correctedK3)  spAdult[parameters](patchesPosition[i]) = correctedK3;
				}

				AreaCount(patchesPosition[i])+= ((spAdult[parameters](patchesPosition[i])*((pow(parameters[16],0.75))*0.0001125))- areaAdultTminus1) -areaRecruitTminus1;
				//updating K
				if(AreaCount(patchesPosition[i])>maxK) K(patchesPosition[i])=0;
				else K(patchesPosition[i])=K(patchesPosition[i]) - (((spAdult[parameters](patchesPosition[i])*((pow(parameters[16],0.75))*0.0001125))- areaAdultTminus1) -areaRecruitTminus1);
					
				if(K(patchesPosition[i])>maxK) K(patchesPosition[i])=maxK;
				if(K(patchesPosition[i])<0) K(patchesPosition[i])=0;
			}
		}
		
		if (anaspeciation=="on"){
			bool contained = immGeneflow.find(parameters) != immGeneflow.end();
			if((contained)&&(immGeneflow[parameters](0,1)>=1)) {
				//updating gene flow map for anagenesis:   suminitialSeeds,sumrecruitingSeeds,suminitialRecruits,summaturingRecruits
					double immigrantrecruiting = (immGeneflow[parameters](0,1)/suminitialSeeds)*sumrecruitingSeeds;
					double immigrantmaturing = (immigrantrecruiting/suminitialRecruits)*summaturingRecruits;
					double GenerationTime = (pow(parameters[16],0.25)*exp(E/(Boltzmann*(temperature-parameters[19])))*allometabolConstTree[12])*0.0001;
					immGeneflow[parameters](0,0) =immGeneflow[parameters](0,0) +10*GenerationTime*immigrantmaturing;
					double maxanatime = timeStep + (pow(parameters[16],0.25)*exp(E/(Boltzmann*(temperature-parameters[19])))*allometabolConstTree[12]);
					if (immGeneflow[parameters](0,0)>maxanatime) immGeneflow[parameters](0,0) =maxanatime; //time for anagenesis cannot exceed the number of generations for it
					immGeneflow[parameters](0,1) =0;
			}
		}
	}	
	
	//updating Carrying capacity
	Matrix<float> areaCount (NY,NY);
	for(map<vector<double>, Matrix<float> >::iterator itr = spSeed.begin(); itr != spSeed.end(); ++itr) 
	{
		const vector<double>& parameters = itr->first; // key (== parameterSet)
		Matrix<float>& abundances = itr->second; 
		for(unsigned int i=0; i<patchesPosition.size(); i++) 
	    {
			int c = patchesPosition[i].first;
			int r = patchesPosition[i].second;		
			if(spAdult[parameters](c,r)<0) spAdult[parameters](c,r)=0;	
			if(spRecruit[parameters](c,r)<0) spRecruit[parameters](c,r)=0;	
			if(spSeed[parameters](c,r)<0) spSeed[parameters](c,r)=0;	
			float AreaTotal; 
			AreaTotal =spAdult[parameters](c,r)*((pow(parameters[16],0.75))*0.0001125) +spRecruit[parameters](c,r)*((pow(parameters[20],0.75))*0.0001125);
			areaCount(c,r)+=AreaTotal;
		}
	}
	for(unsigned int i=0; i<patchesPosition.size(); i++) 
    {
		int c = patchesPosition[i].first;
		int r = patchesPosition[i].second;		
		K(c,r) = maxK-areaCount(c,r);
		if(K(c,r)<0) K(c,r) =0;	
	}

	//Mutant maps:
	speciesAdultMutant.clear();
	speciesRecruitMutant.clear();
	speciesSeedMutant.clear();
	speciesDispMutant.clear();
    speciesHabtMutant.clear();
	
	//reproduction and speciation:
	for(map<vector<double>, Matrix<float> >::iterator itr = spAdult.begin(); itr != spAdult.end(); ++itr) {
		const vector<double>& parameters = itr->first; // key (== parameterSet)
		Matrix<float>& abundances = itr->second; // value (== abundances)
		 
		//Using matrix addition of library armadillo (a lot faster):	
		mat D(NY*2 +1,NY*2 +1);
		D.fill(0.0);
		mat A(NY,NY);
		A.fill(0.0);
		for (int co=0; co<(int)NY*2 +1; co++ )
		 {
		 	for (int ro=0; ro<(int)NY*2+1; ro++ )
		 	{						
				D(co,ro) =spDisp[parameters]( abs((int)NY +1 - co) , abs((int)NY+1 - ro) );					
			}
		 } 
         
		for (unsigned int i=0; i<patchesPosition.size(); i++ )
		{
    		// Beverton-Holt Modell extende for Allee effects
			float spTotalExploitedArea; 
			double Ki;			
			spTotalExploitedArea= (spRecruit[parameters](patchesPosition[i])*((pow(parameters[20],0.75))*0.0001125))+spAdult[parameters](patchesPosition[i])*((pow(parameters[16],0.75))*0.0001125);	
			Ki= ((K(patchesPosition[i])+spTotalExploitedArea)*spHabt[parameters](patchesPosition[i]))/((pow(parameters[20],0.75))*0.0001125);			
			double C = 0;			
			if(parameters[3]<0) C= parameters[3]*((maxK/((pow(parameters[16],0.75))*0.0001125))*spHabt[parameters](patchesPosition[i]));
			else C= parameters[3];
			
			//making minimum value of C being -K
			if(C<-Ki) C=-Ki;
			if (C==Ki) C= C-1;
			//Beverton-Holt carrying capacity for nonsprouters:
			double R, k, c;
			R = parameters[1];
			k = 4*(parameters[1]-parameters[2])/(parameters[2]*pow((Ki-C),2));
			c = C + sqrt((parameters[1]-parameters[2])/(parameters[2]*k));
					
			//Beverton-Holt:
			double dispersalUnits; 
			dispersalUnits= (abundances(patchesPosition[i])*R)/(1+ k*(pow((abundances(patchesPosition[i]) - c),2)));			
			
			int is = patchesPosition[i].first;
			int iz = patchesPosition[i].second;			
			//Use a dispersal kernel to disperse the dispersal units (use a dispersal funtion)           
			if(dispersalUnits>0){						
				//Using matrix addition of library armadillo (a lot faster):	
				A=A+  D.submat(NY+1-is,NY+1-iz,NY+(NY-is),NY+(NY-iz))*dispersalUnits;
				//A.print("A:");			
			}						
		}
		
		//use the incoming dispersal units to draw a random recruitment
		for(unsigned int i=0; i<patchesPosition.size(); i++)
		{
			int is = patchesPosition[i].first;
			int iz = patchesPosition[i].second;

			float mutants = 0;
            if (cladospeciation=="on"){
				//Speciation: random recruitment of mutants				
				//Using matrix addition of library armadillo:					
				mutants = (float) Zufall.poisson(A(is,iz)*(1 - exp(-parameters[5])));
				
				//if at least one mutant of this species recruited (only one mutant per timestep and island cell allowed
				//- all the mutants from the "same" diverging species):
				if(mutants>0)
				{
					/*cout<<"Mutation rate " << parameters[5] <<endl; 
					cout<<"Mutation probability " << 1 - exp(-parameters[5]) <<endl; 
					cout<<"Numer of mutants" << mutants <<endl; */	

					vector<double> parameterSetMutant;					
						//choosing random parameter values for hte new mutant +- half of the value
						//Relative Fenology:
						double max = (parameters[0]+ parameters[0]*nicheEvolution)*10000 ;
						double min = (parameters[0]- parameters[0]*nicheEvolution)*10000 ;
						double Fenology = rand()%( (int)(max-min));
						Fenology = min + Fenology;
						Fenology = Fenology*0.0001;
						//Rmax:
						max = (parameters[1]+ parameters[1]*nicheEvolution) *10000;
						min = (parameters[1]- parameters[1]*nicheEvolution) *10000;
						double Rmax = rand()%( (int)(max-min));
						Rmax = min + Rmax;
						Rmax = Rmax*0.0001;
						if (Rmax<minValueTree[1]) Rmax=minValueTree[1];
						if (Rmax>maxValueTree[1]) Rmax=maxValueTree[1]; 
						//M:
						max = (parameters[2]+ parameters[2]*nicheEvolution) *10000;
						min = (parameters[2]- parameters[2]*nicheEvolution) *10000;
						double M = rand()%( (int)(max-min));
						M = min + M;
						M = M*0.0001;
						if (M<minValueTree[2]) M=minValueTree[2];
						if (M>maxValueTree[2]) M=maxValueTree[2]; 
						//Gamma:
						max = (parameters[3]+ parameters[3]*nicheEvolution)*10000;
						min = (parameters[3]- parameters[3]*nicheEvolution)*10000;
						double Gamma;
						if((max-min)==0) Gamma = rand()%( (int) (0.5*10000) );
						else Gamma = rand()%( (int)(max-min));						 
						Gamma = min + Gamma;
						Gamma = Gamma*0.0001;
                        if (Gamma<minValueTree[3]) Gamma=minValueTree[3];
						if (Gamma>maxValueTree[3]) Gamma=maxValueTree[3]; 
						//GenTime:
						max = (parameters[4]+ parameters[4]*nicheEvolution) *10000;
						min = (parameters[4]- parameters[4]*nicheEvolution) *10000;
						double GenTime = rand()%( (int)(max-min));
						GenTime = min + GenTime;
						GenTime = GenTime*0.0001;
						//MutRate:
						double MutRate;
						
						max = (parameters[5]+ parameters[5]*nicheEvolution);
						min = (parameters[5]- parameters[5]*nicheEvolution);
						MutRate = rand()%( (int)((max-min)*1000000000000000));
						/*cout <<"Mutation rate: " << MutRate <<endl;
						cout <<"Mutation rate: " << MutRate <<endl;
						cout <<"Mutation rate: " << MutRate <<endl;*/
						MutRate = min + MutRate;
						MutRate = MutRate*0.000000000000001;
						if (MutRate<minValueTree[5]) MutRate=minValueTree[5];
						if (MutRate>maxValueTree[5]) MutRate=maxValueTree[5]; 

						//Dispersal Ability (mean dispersal distance):
						max = (parameters[6]+ parameters[6]*nicheEvolution);
						min = (parameters[6]- parameters[6]*nicheEvolution);

						double Alpha = rand()%( (int)((max-min)*1000000));
						Alpha = min + Alpha;
						Alpha = Alpha*0.000001;
						if (Alpha<minValueTree[6]) Alpha=minValueTree[6];
						if (Alpha>maxValueTree[6]) Alpha=maxValueTree[6]; 

						//Dispersal Ability (P):
						max = (parameters[7]+ parameters[7]*nicheEvolution) *10000;
						min = (parameters[7]- parameters[7]*nicheEvolution) *10000;
						double P = rand()%( (int)(max-min));
						P = min + P;
						P = P*0.0001;
						if (P<minValueTree[7]) P=minValueTree[7];
						if (P>maxValueTree[7]) P=maxValueTree[7]; 
						//K PFTFACTOR:
						double Kpftfactor  = parameters[8];
						double KpftfactorSeedling = parameters[22];

						double randomTrend = (rand()% 50 +1)*0.01+ 0.5;
						double islandSide= (rand()% 4) +1;

						//Growth:
						max = (parameters[11]+ parameters[11]*nicheEvolution) *10000;
						min = (parameters[11]- parameters[11]*nicheEvolution) *10000;
						double Growth = rand()%( (int)(max-min));
						Growth = min + Growth;
						Growth = Growth*0.0001;
						if (Growth<minValueTree[8]) Growth=minValueTree[8];
						if (Growth>maxValueTree[8]) Growth=maxValueTree[8]; 
						//Mrecruit:
						max = (parameters[12]+ parameters[12]*nicheEvolution) *10000;
						min = (parameters[12]- parameters[12]*nicheEvolution) *10000;
						double Mrecruit = rand()%( (int)(max-min));
						Mrecruit = min + Mrecruit;
						Mrecruit = Mrecruit*0.0001;
						if (Mrecruit<minValueTree[9]) Mrecruit=minValueTree[9];
						if (Mrecruit>maxValueTree[9]) Mrecruit=maxValueTree[9]; 
						//Germ:
						max = (parameters[13]+ parameters[13]*nicheEvolution) *10000;
						min = (parameters[13]- parameters[13]*nicheEvolution) *10000;
						double Germ = rand()%( (int)(max-min));
						Germ = min + Germ;
						Germ = Germ*0.0001;
						if (Germ<minValueTree[10]) Germ=minValueTree[10];
						if (Germ>maxValueTree[10]) Germ=maxValueTree[10];
						//Mseed:
						max = (parameters[14]+ parameters[14]*nicheEvolution) *10000;
						min = (parameters[14]- parameters[14]*nicheEvolution) *10000;
						double Mseed = rand()%( (int)(max-min));
						Mseed = min + Mseed;
						Mseed = Mseed*0.0001;
						if (Mseed<minValueTree[11]) Germ=minValueTree[11];
						if (Mseed>maxValueTree[11]) Germ=maxValueTree[11];

						double Annual = rand()% 2;
						//cout<<"Gamma mutant:   "<< Gamma<< endl;
						//cout<<"Kpftfactor mutant:   "<< Kpftfactor<< endl;
				    
						double optimalAltitude = rand()% (maxIslandRadius +1);
						//specialist: 1; generalist: 4:
						double generalist = (rand()% 4) +1;
						//vertical specialist: 1; generalist: verticalpreference:
						double verticalgeneralist = (rand()% verticalpreference) +1;

						
						double MassSeedling, MassSeed;
						//allometabolConst[0] gives the PFT (1=tree, 2=shrub, 3=herb):
						//Mass:
						 max = (parameters[16]+ parameters[16]*nicheEvolution)*10;
						 min = (parameters[16]- parameters[16]*nicheEvolution)*10;
						double Mass = rand()%( (int)(max-min));
						//cout << "min Mass: "<<min<<"max Mass: " <<max<<"random Mass: " <<Mass<<endl;
						Mass = min + Mass;
						Mass = Mass*0.1;
						if(Mass<10) Mass = 10; //lower bound for herbs
						if (Mass>55000)  {
							MassSeedling = Mass*0.2;
							MassSeed = (exp((rand()% 80)*0.1))*0.01; 			
						}
						else {
							MassSeedling = Mass*0.2;
							MassSeed= pow(Mass*0.0000046,0.5);//0.00043*Mass;			
						}

						double yearOfSpec =timeStep + (pow(Mass,0.25)*exp(E/(Boltzmann*tempForFixedMetabolism))*allometabolConstTree[12]);
						double SpType=1;
						if(parameters[25]==4) SpType =5;

						parameterSetMutant.push_back(Fenology);
  						parameterSetMutant.push_back(Rmax);
						parameterSetMutant.push_back(M);
						parameterSetMutant.push_back(Gamma);
						parameterSetMutant.push_back(GenTime);
						parameterSetMutant.push_back(MutRate);
						parameterSetMutant.push_back(Alpha);
						parameterSetMutant.push_back(P);
						parameterSetMutant.push_back(Kpftfactor);
						parameterSetMutant.push_back(randomTrend);
						parameterSetMutant.push_back(islandSide);
						parameterSetMutant.push_back(Growth);
						parameterSetMutant.push_back(Mrecruit);
						parameterSetMutant.push_back(Germ);
						parameterSetMutant.push_back(Mseed);
						parameterSetMutant.push_back(Annual);
						parameterSetMutant.push_back(Mass); //Dummy key for Mass in the Metabolic version
						parameterSetMutant.push_back(parameters[17]); //lineage
						parameterSetMutant.push_back(parameters[17]*1000+parameters[16]); //MotherID 
						parameterSetMutant.push_back(optimalAltitude); //lineage
						parameterSetMutant.push_back(MassSeedling); //Dummy key =Massseedling for hte metabolic version
						parameterSetMutant.push_back(MassSeed); //Dummy key =Massseed for hte metabolic version
						parameterSetMutant.push_back(KpftfactorSeedling);
						parameterSetMutant.push_back(generalist);
						parameterSetMutant.push_back(verticalgeneralist);//Dummy key =Mass for hte metabolic version
						parameterSetMutant.push_back(SpType); //0= species from the initial species pool
						parameterSetMutant.push_back(yearOfSpec); //time step of mutation					

						if(parameterSetMutant[16]<=15.0) mutants=0;					
						else
						{
							addRandomSpeciesIfNew(parameterSetMutant,speciesSeedMutant, speciesDispMutant, speciesHabtMutant);
							speciesSeedMutant[parameterSetMutant](patchesPosition[i]) = mutants;
						 }			
				}
			}
		    //Recruitment:
			float seeds;
			seeds = (float)  Zufall.poisson(A(is,iz)-mutants);
			
			//updating seedbank:			
			if(vulnerableSpecialists=="no"){
				if((parameters[24]>=1)&&(spHabt[parameters](patchesPosition[i])>0)) spSeed[parameters](patchesPosition[i])=spSeed[parameters](patchesPosition[i])+seeds;
				if(parameters[24]<=1) spSeed[parameters](patchesPosition[i])=spSeed[parameters](patchesPosition[i])+seeds;
			}
			if(vulnerableSpecialists=="yes"){ //Attention: this makes altitude specialists vulnerable to extinction during island dynamics!!!!!!!!!!!!!!!
				if(spHabt[parameters](patchesPosition[i])>0) spSeed[parameters](patchesPosition[i])=spSeed[parameters](patchesPosition[i])+seeds;
			}

			//Adult mortality:
			float survivors;
			survivors = (float) Zufall.binomial((1-(1 - exp(-parameters[2]))),(int)abundances(patchesPosition[i]));
			if((parameters[16]<1000)&&(parameters[15]==1))survivors = 0;			

			//update abundance and K:
			float areaTimeminus1 ;
			areaTimeminus1=abundances(patchesPosition[i])*((pow(parameters[16],0.75))*0.0001125);
				
			//double KbeforeRecruit = K(patchesPosition[i]);
			abundances(patchesPosition[i])= survivors;			
			AreaCount(patchesPosition[i])-=(areaTimeminus1-(abundances(patchesPosition[i])*((pow(parameters[16],0.75))*0.0001125)));
			//updating K:
			if(AreaCount(patchesPosition[i])>maxK) K(patchesPosition[i])=0;
			else K(patchesPosition[i])=K(patchesPosition[i]) + (areaTimeminus1-(abundances(patchesPosition[i])*((pow(parameters[16],0.75))*0.0001125)));					
			if(K(patchesPosition[i])>maxK) K(patchesPosition[i])=maxK;
            if(K(patchesPosition[i])<0) K(patchesPosition[i])=0;
		}			
	}
		
    addSpeciesIfcolonizing(speciesSeedMutant, spAdult, spRecruit,spSeed, K, speciesDispMutant, spDisp,speciesHabtMutant,spHabt);
	
	//updating species status (from pre-cladospecies to new evovled cladospecies): 1->2
	if (cladospeciation=="on"){
		SpecClado=0; 
		LookforCladoSpecies (spAdult,spRecruit,spSeed,spDisp,spHabt);
	}
	if (anaspeciation=="on"){
		SpecAna=0;
		LookforAnaSpecies (spAdult,spRecruit,spSeed,spDisp,spHabt,immGeneflow);
	}
}
void Simulation::LookforCladoSpecies(map<vector<double>, Matrix<float> >& spAdult, 
	map<vector<double>, Matrix<float> >& spRecruit,	map<vector<double>, Matrix<float> >& spSeed,
	map<vector<double>, Matrix<float> >& spDisp, map<vector<double>, Matrix<float> >& spHabt)
{
   // we need a "special" loop (with "while") when an element of the map may be erased
	map<vector<double>, Matrix<float> >::iterator itr = spSeed.begin();
	map<vector<double>, Matrix<float> >::iterator itrr = spRecruit.begin();
	map<vector<double>, Matrix<float> >::iterator itra = spAdult.begin();
	map<vector<double>, Matrix<float> >::iterator it = spDisp.begin();
	map<vector<double>, Matrix<float> >::iterator itH = spHabt.begin();
	
	map<vector<double>, Matrix<float> > newspeciesAdult;
	map<vector<double>, Matrix<float> > newspeciesRecruit;
	map<vector<double>, Matrix<float> > newspeciesSeed;
	map<vector<double>, Matrix<float> > newspeciesDisp;
	map<vector<double>, Matrix<float> > newspeciesHabt;
	
	unsigned int Initialrichness =spSeed.size();
	
	while(itr != spSeed.end()) {

		const vector<double>& parameters = itr->first; // key (== parameterSet)
		Matrix<float>& abundances = itr->second; // value (== abundances)
		
		if((parameters[26]>0)&&((int)parameters[26]==timeStep)) {
			vector<double> parameterSetNewcladoSpecies;
			double SpType=2;
			if(parameters[25]==5) SpType =6;

			parameterSetNewcladoSpecies.push_back(parameters[0]);
  			parameterSetNewcladoSpecies.push_back(parameters[1]);
			parameterSetNewcladoSpecies.push_back(parameters[2]);
			parameterSetNewcladoSpecies.push_back(parameters[3]);
			parameterSetNewcladoSpecies.push_back(parameters[4]);
			parameterSetNewcladoSpecies.push_back(parameters[5]);
			parameterSetNewcladoSpecies.push_back(parameters[6]);
			parameterSetNewcladoSpecies.push_back(parameters[7]);
			parameterSetNewcladoSpecies.push_back(parameters[8]);
			parameterSetNewcladoSpecies.push_back(parameters[9]);
			parameterSetNewcladoSpecies.push_back(parameters[10]);
			parameterSetNewcladoSpecies.push_back(parameters[11]);
			parameterSetNewcladoSpecies.push_back(parameters[12]);
			parameterSetNewcladoSpecies.push_back(parameters[13]);
			parameterSetNewcladoSpecies.push_back(parameters[14]);
			parameterSetNewcladoSpecies.push_back(parameters[15]);
			parameterSetNewcladoSpecies.push_back(parameters[16]);
			parameterSetNewcladoSpecies.push_back(parameters[17]); //lineage
			parameterSetNewcladoSpecies.push_back(parameters[18]); //MotherID 
			parameterSetNewcladoSpecies.push_back(parameters[19]); //lineage
			parameterSetNewcladoSpecies.push_back(parameters[20]); //Dummy key =Mass for hte metabolic version
			parameterSetNewcladoSpecies.push_back(parameters[21]); //Dummy key =Mass for hte metabolic version
			parameterSetNewcladoSpecies.push_back(parameters[22]);
			parameterSetNewcladoSpecies.push_back(parameters[23]);
			parameterSetNewcladoSpecies.push_back(parameters[24]);//Dummy key =Mass for hte metabolic version
			parameterSetNewcladoSpecies.push_back(SpType); //2= new clado species; if 5-<endemic from an invasive
			parameterSetNewcladoSpecies.push_back(parameters[26]); //time step of mutation
			//include it in one col-like map:
			newspeciesAdult[parameterSetNewcladoSpecies] = spAdult[parameters];
			newspeciesRecruit[parameterSetNewcladoSpecies] = spRecruit[parameters];
			newspeciesSeed[parameterSetNewcladoSpecies] = spSeed[parameters];
			newspeciesDisp[parameterSetNewcladoSpecies] = spDisp[parameters];
			newspeciesHabt[parameterSetNewcladoSpecies] = spHabt[parameters];
			// remove it
			itr = spSeed.erase(itr);		
		}
		else {
			++itr;
		}
	}
	//counting number of fully fledged new species:
	if(spSeed.size()<Initialrichness) SpecClado=Initialrichness-spSeed.size();

	//remove it from other old maps:
	while(itrr != spRecruit.end()) {
		const vector<double>& parameters = itrr->first; // key (== parameterSet)
		Matrix<float>& abundances = itrr->second; // value (== abundances)		
		bool contained = spSeed.find(parameters) != spSeed.end();
		if(!contained) {
			// remove it
			itrr= spRecruit.erase(itrr);			
		}
		else {
			++itrr;
		}
	}

	while(itra != spAdult.end()) {
		const vector<double>& parameters = itra->first; // key (== parameterSet)
		Matrix<float>& abundances = itra->second; // value (== abundances)		
		bool contained = spSeed.find(parameters) != spSeed.end();
		if(!contained) {
			// remove it
			itra= spAdult.erase(itra);			
		}
		else {
			++itra;
		}
	}

	while(it != spDisp.end()) {
		const vector<double>& parameters = it->first; // key (== parameterSet)
		Matrix<float>& abundances = it->second; // value (== abundances)		
		bool contained = spSeed.find(parameters) != spSeed.end();
		if(!contained) {
			// remove it
			it= spDisp.erase(it);			
		}
		else {
			++it;
		}
	}

	while(itH != spHabt.end()) {
		const vector<double>& parameters = itH->first; // key (== parameterSet)
		Matrix<float>& abundances = itH->second; // value (== abundances)		
		bool contained = spSeed.find(parameters) != spSeed.end();
		if(!contained) {
			// remove it
			itH= spHabt.erase(itH);			
		}
		else {
			++itH;
		}
	}
		
	//re-add to old maps with the parameter key updated:
	for(map<vector<double>, Matrix<float> >::iterator itrSpec = newspeciesSeed.begin(); itrSpec != newspeciesSeed.end(); ++itrSpec) {
		const vector<double>& parameters = itrSpec->first; // key (== parameterSet)
		Matrix<float>& abundances = itrSpec->second; // value (== abundances)
		spSeed[parameters] = abundances;
		spRecruit[parameters] = newspeciesRecruit[parameters];
		spAdult[parameters] = newspeciesAdult[parameters];
		spDisp[parameters] = newspeciesDisp[parameters];
		spHabt[parameters] = newspeciesHabt[parameters];
	}
}
void Simulation::LookforAnaSpecies(map<vector<double>, Matrix<float> >& spAdult, 
	map<vector<double>, Matrix<float> >& spRecruit,	map<vector<double>, Matrix<float> >& spSeed,
	map<vector<double>, Matrix<float> >& spDisp, map<vector<double>, Matrix<float> >& spHabt, map<vector<double>, Matrix<float> >& immGeneflow)
{

	// we need a "special" loop (with "while") when an element of the map may be erased
	map<vector<double>, Matrix<float> >::iterator itr = spSeed.begin();
	map<vector<double>, Matrix<float> >::iterator itrr = spRecruit.begin();
	map<vector<double>, Matrix<float> >::iterator itra = spAdult.begin();
	map<vector<double>, Matrix<float> >::iterator it = spDisp.begin();
	map<vector<double>, Matrix<float> >::iterator itH = spHabt.begin();
	map<vector<double>, Matrix<float> >::iterator itG = immGeneflow.begin();

	map<vector<double>, Matrix<float> > newspeciesAdult;
	map<vector<double>, Matrix<float> > newspeciesRecruit;
	map<vector<double>, Matrix<float> > newspeciesSeed;
	map<vector<double>, Matrix<float> > newspeciesDisp;
	map<vector<double>, Matrix<float> > newspeciesHabt;
	
	unsigned int Initialrichness =spSeed.size();

	while(itr != spSeed.end()) {

		const vector<double>& parameters = itr->first; // key (== parameterSet)
		Matrix<float>& abundances = itr->second; // value (== abundances)
		bool contained = immGeneflow.find(parameters) != immGeneflow.end();
		if((contained)&&(timeStep>= immGeneflow[parameters](0,0))) {
			vector<double> parameterSetNewanaSpecies;
			double SpType=3;

			parameterSetNewanaSpecies.push_back(parameters[0]-1); //the anagenetic species has a phenological advantage to future mainland immigrants. 
  			parameterSetNewanaSpecies.push_back(parameters[1]);
			parameterSetNewanaSpecies.push_back(parameters[2]);
			parameterSetNewanaSpecies.push_back(parameters[3]);
			parameterSetNewanaSpecies.push_back(parameters[4]);
			parameterSetNewanaSpecies.push_back(parameters[5]);
			parameterSetNewanaSpecies.push_back(parameters[6]);
			parameterSetNewanaSpecies.push_back(parameters[7]);
			parameterSetNewanaSpecies.push_back(parameters[8]);
			parameterSetNewanaSpecies.push_back(parameters[9]);
			parameterSetNewanaSpecies.push_back(parameters[10]);
			parameterSetNewanaSpecies.push_back(parameters[11]);
			parameterSetNewanaSpecies.push_back(parameters[12]);
			parameterSetNewanaSpecies.push_back(parameters[13]);
			parameterSetNewanaSpecies.push_back(parameters[14]);
			parameterSetNewanaSpecies.push_back(parameters[15]);
			parameterSetNewanaSpecies.push_back(parameters[16]);
			parameterSetNewanaSpecies.push_back(parameters[17]); //lineage
			parameterSetNewanaSpecies.push_back(parameters[18]); //MotherID 
			parameterSetNewanaSpecies.push_back(parameters[19]); //lineage
			parameterSetNewanaSpecies.push_back(parameters[20]); // key =Mass for hte metabolic version
			parameterSetNewanaSpecies.push_back(parameters[21]); // key =Mass for hte metabolic version
			parameterSetNewanaSpecies.push_back(parameters[22]);
			parameterSetNewanaSpecies.push_back(parameters[23]);
			parameterSetNewanaSpecies.push_back(parameters[24]);// key =Mass for hte metabolic version
			parameterSetNewanaSpecies.push_back(SpType); //3= new ana species; 
			parameterSetNewanaSpecies.push_back(parameters[26]); //time step of mutation
			//include it in one col-like map:
			newspeciesAdult[parameterSetNewanaSpecies] = spAdult[parameters];
			newspeciesRecruit[parameterSetNewanaSpecies] = spRecruit[parameters];
			newspeciesSeed[parameterSetNewanaSpecies] = spSeed[parameters];
			newspeciesDisp[parameterSetNewanaSpecies] = spDisp[parameters];
			newspeciesHabt[parameterSetNewanaSpecies] = spHabt[parameters];
			// remove it
			itr = spSeed.erase(itr);		
		}
		else {
			++itr;
		}
	}
	//counting number of fully fledged new species:
	if(spSeed.size()<Initialrichness) SpecAna=Initialrichness-spSeed.size();

	//remove it from other old maps:
	while(itrr != spRecruit.end()) {
		const vector<double>& parameters = itrr->first; // key (== parameterSet)
		Matrix<float>& abundances = itrr->second; // value (== abundances)		
		bool contained = spSeed.find(parameters) != spSeed.end();
		if(!contained) {
			// remove it
			itrr= spRecruit.erase(itrr);			
		}
		else {
			++itrr;
		}
	}

	while(itra != spAdult.end()) {
		const vector<double>& parameters = itra->first; // key (== parameterSet)
		Matrix<float>& abundances = itra->second; // value (== abundances)		
		bool contained = spSeed.find(parameters) != spSeed.end();
		if(!contained) {
			// remove it
			itra= spAdult.erase(itra);			
		}
		else {
			++itra;
		}
	}

	while(it != spDisp.end()) {
		const vector<double>& parameters = it->first; // key (== parameterSet)
		Matrix<float>& abundances = it->second; // value (== abundances)		
		bool contained = spSeed.find(parameters) != spSeed.end();
		if(!contained) {
			// remove it
			it= spDisp.erase(it);			
		}
		else {
			++it;
		}
	}

	while(itH != spHabt.end()) {
		const vector<double>& parameters = itH->first; // key (== parameterSet)
		Matrix<float>& abundances = itH->second; // value (== abundances)		
		bool contained = spSeed.find(parameters) != spSeed.end();
		if(!contained) {
			// remove it
			itH= spHabt.erase(itH);			
		}
		else {
			++itH;
		}
	}
	
	while(itG != immGeneflow.end()) {
		const vector<double>& parameters = itG->first; // key (== parameterSet)
		Matrix<float>& abundances = itG->second; // value (== abundances)		
		bool contained = spSeed.find(parameters) != spSeed.end();
		if(!contained) {
			// remove it
			itG= immGeneflow.erase(itG);			
		}
		else {
			++itG;
		}
	}
	
	//re-add to old maps with the parameter key updated:
	for(map<vector<double>, Matrix<float> >::iterator itrSpec = newspeciesSeed.begin(); itrSpec != newspeciesSeed.end(); ++itrSpec) {
		const vector<double>& parameters = itrSpec->first; // key (== parameterSet)
		Matrix<float>& abundances = itrSpec->second; // value (== abundances)
		spSeed[parameters] = abundances;
		spRecruit[parameters] = newspeciesRecruit[parameters];
		spAdult[parameters] = newspeciesAdult[parameters];
		spDisp[parameters] = newspeciesDisp[parameters];
		spHabt[parameters] = newspeciesHabt[parameters];
	}
}
void Simulation::islanddynamics(map<vector<double>, Matrix<float> >& spAdult, map<vector<double>, Matrix<float> >& spRecruit,
				map<vector<double>, Matrix<float> >& spSeed, Matrix<float>& K, 
	            map<vector<double>, Matrix<float> >& spDisp, map<vector<double>, Matrix<float> >& spHabt, map<vector<double>, Matrix<float> >& spHabtMain, Matrix<float>& TempKelvinBareIsland)
{
	//checking carrying capacity:
	/*Matrix<float> abundCount (NY,NY);
	for(map<vector<double>, Matrix<float> >::iterator itr = spSeed.begin(); itr != spSeed.end(); ++itr) 
	{
		const vector<double>& parameters = itr->first; // key (== parameterSet)
		Matrix<float>& abundances = itr->second; 
		for (unsigned int r=0;r<NY;++r)
		{
			for (unsigned int c=0;c<NY;++c)
			{
					float AbundTotal =spAdult[parameters](c,r)*(float)parameters[8] +spRecruit[parameters](c,r)*(float)parameters[22];
					abundCount(c,r)+=AbundTotal;
			}
		}
	}
	for (unsigned int r=0;r<NY;++r)
	{
		for (unsigned int c=0;c<NY;++c)
		{
			if(State(c,r)>=1){
				//K(c,r) = maxK-abundCount(c,r);
				//if(K(c,r)<0) K(c,r) =0;

				cout << "c: " << c <<endl;
				cout << "r: " << r <<endl;
				cout << "abundCount(c,r): " << abundCount(c,r) <<endl;
				cout << "K(c,r): " << K(c,r) <<endl;
			}
		}
	}*/
		
	DemExt=0;
	ErosionExt=0;
	VulcanicExt=0;
	RandomExt=0;	
	
	unsigned int InitialGammaRichness=spAdult.size();
	//cout << "InitialGammaRichness: " << InitialGammaRichness <<endl;
     //these both for counting pre-cladospecies and not consider them in the extinction count
	prespeciesExtinct=0;
    
	// erasing possible extincted species:
	removeExtinctSpecies(spAdult,spRecruit, spSeed, spDisp,spHabt);	

	DemExt= InitialGammaRichness - (spAdult.size()+prespeciesExtinct); //Extinction due to demographic stochasticity

	//If climate change is switched on:
	if ((timeStep>=climateChange)&&(timeStep< (climateChange + temperatureChangeInterval)))
	{
		for (unsigned int r=0;r<NY;++r)
		{
			for (unsigned int c=0;c<NY;++c)
			{
					TempKelvin(c,r)=TempKelvin(c,r)+(temperatureChange/temperatureChangeInterval);
					TempKelvinBareIsland(r,c)=TempKelvinBareIsland(r,c)+(temperatureChange/temperatureChangeInterval);
			}
		}
	}


	//Random disturbance:
	float distprob= rand()% 1000*(float) 0.001;
	if ((randomDisturbance==1) && (distprob<=disturbanceProb))
	{
		disturbOccurrence=1;
		int disturbanceRadius = rand () % (maxIslandRadius*2); // is able to reach the whole islands/mountain
		disturbRadius=disturbanceRadius;
		unsigned int randomCoord = rand () % patchesPosition.size();
		unsigned int randomx= patchesPosition[randomCoord].first;
		epicenterX=randomx;
		unsigned int randomy= patchesPosition[randomCoord].second;
		epicenterY=randomy;
		for(map<vector<double>, Matrix<float> >::iterator itr = spAdult.begin(); itr != spAdult.end(); ++itr) {
			const vector<double>& parameters = itr->first; 
			float disturbance = (rand() % 100)*(float)0.01;
			//disturbance to seeds milder:
			//float distRandom = disturbance*100;
			//if (distRandom<2) distRandom =1;
			//float disturbanceSeeds = (rand() % (int)distRandom)*(float)0.01;		
			for (unsigned int x=randomx-disturbanceRadius;x<=randomx+disturbanceRadius;++x)
			{
				for (unsigned int y=randomy-disturbanceRadius;y<=randomy+disturbanceRadius;++y)
				{
					if (((x>0)&&(x<NY))&&((y>0)&&(y<NY))){
						if(State(y,x)>=1){			
							float AreaTotalbeforeDistbrub =(spAdult[parameters](x,y)*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(x,y))))*allometabolConstTree[13]))) +spRecruit[parameters](x,y)*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*TempKelvin(x,y))))*allometabolConstTree[13]));
							spAdult[parameters](x,y) =spAdult[parameters](x,y) -spAdult[parameters](x,y)*disturbance;
							spRecruit[parameters](x,y) =spRecruit[parameters](x,y) -spRecruit[parameters](x,y)*disturbance;
							//spSeed[parameters](x,y) =spSeed[parameters](x,y) -spSeed[parameters](x,y)*disturbanceSeeds; //disturbance to seeds milder		
							spSeed[parameters](x,y) =spSeed[parameters](x,y) -spSeed[parameters](x,y)*disturbance;
							float AreaTotalafterDistbrub =(spAdult[parameters](x,y)*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(x,y))))*allometabolConstTree[13]))) +spRecruit[parameters](x,y)*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*TempKelvin(x,y))))*allometabolConstTree[13]));
							//update K
							K(x,y)=K(x,y) + AreaTotalbeforeDistbrub - AreaTotalafterDistbrub;
							if(K(x,y)>maxK) K(x,y)=maxK;
							if(K(x,y)<0) K(x,y)=0;
						}
					}
				}
			}
		}
	}
	unsigned int InitialGammaRichnessR=spAdult.size();
	prespeciesExtinct=0;
	removeExtinctSpecies(spAdult,spRecruit, spSeed, spDisp,spHabt);
	RandomExt= InitialGammaRichnessR - (spAdult.size()+prespeciesExtinct);
	
	//Autocorrelated random variant:
	if (habSuitability=="Random")
	{
		//Island growth:   
		if ((timeStep==islandTimestep*growthStep)&&(growthStep<=maxIslandRadius))
		{
			for (unsigned int i=0; i<patchesPosition.size(); i++ )
			{
				unsigned int is = patchesPosition[i].first;
				unsigned int iz = patchesPosition[i].second;

				for (unsigned int ib=(is-1); ib<(is+2); ib++ )
				{
					for (unsigned int ia=(iz-1); ia<(iz+2); ia++ )
					{
						if(State(ib,ia)==-99)
						{
							State(ib,ia)=1;
							K(ib,ia)= maxK;
						}
					}
				}

			}
			growthStep++;

			//Updating habitat suitability island:
			for(map<vector<double>, Matrix<float> >::iterator itr = spHabt.begin(); itr != spHabt.end(); ++itr) 
			{
				const vector<double>& parameters = itr->first; // key (== parameterSet)
				Matrix<float>& habtsuitability = itr->second; // value (== abundances)
				for (unsigned int i=0; i<patchesPosition.size(); i++ )
				{
					habtsuitability(patchesPosition[i]) = habtsuitability(patchesPosition[i])- (float) 0.05;
					if(habtsuitability(patchesPosition[i])<0)habtsuitability(patchesPosition[i])=0;
					//To check the error below: also: the mutants will not be in the mainland!!!
					//spHabtMain[parameters](patchesPosition[i]) = spHabtMain[parameters](patchesPosition[i])-0.025;
					//if(spHabtMain[parameters](patchesPosition[i])<0)spHabtMain[parameters](patchesPosition[i])=0;
				}
			}
			//Updating habitat suitability mainland:
			for(map<vector<double>, Matrix<float> >::iterator itra = spHabtMain.begin(); itra != spHabtMain.end(); ++itra) 
			{
				const vector<double>& parameters = itra->first; // key (== parameterSet)
				Matrix<float>& habtsuitability = itra->second; // value (== abundances)
				bool contained = spHabt.find(parameters) != spHabt.end();
				if(contained) {
					for (unsigned int i=0; i<patchesPosition.size(); i++ )
					{
						habtsuitability(patchesPosition[i]) = habtsuitability(patchesPosition[i])- (float) 0.05;
						if(habtsuitability(patchesPosition[i])<0)habtsuitability(patchesPosition[i])=0;
					}
				}
			}

			//Updating habitatPatches:
		   patchesPosition.clear();
			for (unsigned int r=0;r<NY;++r)
			{
				for (unsigned int c=0;c<NX;++c)
				{
					if(State(c,r)>=1) 
					{
						//if the cell content is 1, add the cell to the patches
						patchesPosition.push_back(std::pair<int,int>(c,r));
					}

				}
			}
			//cout << "State: " << State <<endl;
			
			//vulcanic disturbance:
			if (vulcanicDisturbance==1)
			{
				int disturbanceRadius = rand () % growthStep;
				for (unsigned int x=centralx-disturbanceRadius;x<=centralx+disturbanceRadius;++x)
				{
					for (unsigned int y=centraly-disturbanceRadius;y<=centraly+disturbanceRadius;++y)
					{
						float disturbance = (rand() % 100)*(float)0.01;
						for(map<vector<double>, Matrix<float> >::iterator itr = spAdult.begin(); itr != spAdult.end(); ++itr) {
							const vector<double>& parameters = itr->first; // key (== parameterSet)
							
							float AreaTotalbeforeDistbrub =(spAdult[parameters](x,y)*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(x,y))))*allometabolConstTree[13]))) +spRecruit[parameters](x,y)*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*TempKelvin(x,y))))*allometabolConstTree[13]));
							spAdult[parameters](x,y) =spAdult[parameters](x,y) -spAdult[parameters](x,y)*disturbance;
							spRecruit[parameters](x,y) =spRecruit[parameters](x,y) -spRecruit[parameters](x,y)*disturbance;
							spSeed[parameters](x,y) =spSeed[parameters](x,y) -spSeed[parameters](x,y)*disturbance;
							
							float AreaTotalafterDistbrub =(spAdult[parameters](x,y)*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(x,y))))*allometabolConstTree[13]))) +spRecruit[parameters](x,y)*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*TempKelvin(x,y))))*allometabolConstTree[13]));
							//update K
							K(x,y)=K(x,y) + AreaTotalbeforeDistbrub - AreaTotalafterDistbrub;
							if(K(x,y)>maxK) K(x,y)=maxK;
							if(K(x,y)<0) K(x,y)=0;
						}
					}
				}
			}
			unsigned int InitialGammaRichnessV=spAdult.size();
			prespeciesExtinct=0;
			removeExtinctSpecies(spAdult,spRecruit, spSeed, spDisp,spHabt);
			VulcanicExt=InitialGammaRichnessV - (spAdult.size()+prespeciesExtinct); 
		}
	
		//Island erosion:
		if ((timeStep==islandTimestep*(erosionStep+erosionStep2+1))&&(erosionStep<=maxIslandRadius))
		{
			bool eroded =false;
			for (unsigned int i=0; i<patchesPosition.size(); i++ )
			{
				unsigned int is = patchesPosition[i].first;
				unsigned int iz = patchesPosition[i].second;
				if((is==centralx)&&(iz==centraly)&&(eroded ==false)){
					for (unsigned int ib=(is-(maxIslandRadius-(erosionStep-1))); ib<(is+(maxIslandRadius-(erosionStep-2))); ib++ )
					{
						if(State(ib,iz-(maxIslandRadius-(erosionStep-1)))==1)
							{
								State(ib,iz-(maxIslandRadius-(erosionStep-1)))=-99;
								K(ib,iz-(maxIslandRadius-(erosionStep-1)))=0;
							}
						if(State(ib,iz+(maxIslandRadius-(erosionStep-1)))==1)
							{
								State(ib,iz+(maxIslandRadius-(erosionStep-1)))=-99;
								K(ib,iz+(maxIslandRadius-(erosionStep-1)))=0;
							}		
					}
					for (unsigned int ia=(iz-(maxIslandRadius-(erosionStep-1))); ia<(iz+(maxIslandRadius-(erosionStep-2))); ia++ )
					{
						if(State(is-(maxIslandRadius-(erosionStep-1)),ia)==1)
							{
								State(is-(maxIslandRadius-(erosionStep-1)),ia)=-99;
								K(is-(maxIslandRadius-(erosionStep-1)),ia)=0;
							}	
						if(State(is+(maxIslandRadius-(erosionStep-1)),ia)==1)
							{
								State(is+(maxIslandRadius-(erosionStep-1)),ia)=-99;
								K(is+(maxIslandRadius-(erosionStep-1)),ia)=0;
							}	
					}
					eroded=true;
				}
			}
			erosionStep++;
			erosionStep2++;
			patchesPosition.clear();
			for (unsigned int r=0;r<NY;++r)
			{
				for (unsigned int c=0;c<NX;++c)
				{
					if(State(c,r)==1) 
					{
						//if the cell content is 1, add the cell to the patches
						patchesPosition.push_back(std::pair<int,int>(c,r));
					}
				}
			}
			
			for(map<vector<double>, Matrix<float> >::iterator itr = spAdult.begin(); itr != spAdult.end(); ++itr) {
				const vector<double>& parameters = itr->first; // key (== parameterSet)
				Matrix<float>& abundances = itr->second; // value (== abundances)
				for (unsigned int r=0;r<NY;++r)
				{
					for (unsigned int c=0;c<NY;++c)
					{
						if(State(c,r)!=1) 
						{
							abundances(c,r)=0;
							spRecruit[parameters](c,r)=0;
							spSeed[parameters](c,r)=0;
						}
						else
						{
							spHabt[parameters](c,r) = spHabt[parameters](c,r)- (float) 0.05;
							if(spHabt[parameters](c,r)<0)spHabt[parameters](c,r)=0;			
						}
					}
				}
			}
		
			for(map<vector<double>, Matrix<float> >::iterator iter = spHabtMain.begin(); iter != spHabtMain.end(); ++iter) {
				const vector<double>& parameters = iter->first; // key (== parameterSet)
				Matrix<float>& habtsuitability = iter->second; // value (== abundances)
				bool contained = spHabt.find(parameters) != spHabt.end();
				if(contained) {
					for (unsigned int r=0;r<NY;++r)
					{
						for (unsigned int c=0;c<NY;++c)
						{
							if(State(c,r)==1) 
							{
							   habtsuitability(c,r) = habtsuitability(c,r)- (float) 0.05;
							   if(habtsuitability(c,r)<0) habtsuitability(c,r)=0;
							}
						}
					}
				}
			}

		unsigned int InitialGammaRichnessB=spAdult.size();
		prespeciesExtinct=0;
		removeExtinctSpecies(spAdult,spRecruit, spSeed, spDisp,spHabt);

		ErosionExt= InitialGammaRichnessB - (spAdult.size()+prespeciesExtinct); 
		}
	}

	//Autocorrelated random variant:
	if (habSuitability=="Concentric")
	{
		//Island growth:   
		if ((timeStep==islandTimestep*growthStep)&&(growthStep<=maxIslandRadius))
		{
			for (unsigned int r=0;r<NY;++r)
			{
				for (unsigned int c=0;c<NY;++c)
				{
					if(State(c,r)>=1) 
					{
						//State(c,r)=State(c,r)+1;
						if (altitudinalTemp=="on"){
							TempKelvin(c,r)=TempKelvin(c,r)-1;
							TempKelvinBareIsland(r,c)=TempKelvinBareIsland(r,c)-1;
						}
					}
				}
			}
			for (unsigned int i=0; i<patchesPosition.size(); i++ )
			{
				unsigned int is = patchesPosition[i].first;
				unsigned int iz = patchesPosition[i].second;

				//Use a dispersal kernel to disperse the dispersal units (use a dispersal funtion)
				//TODO: substitute for (DX and DY length -1) /2
				for (unsigned int ib=(is-1); ib<(is+2); ib++ )
				{
					for (unsigned int ia=(iz-1); ia<(iz+2); ia++ )
					{
						if(State(ib,ia)==-99)
						{
							State(ib,ia)=1;
							K(ib,ia)= maxK;
							TempKelvin(ib,ia)=temperature;
							TempKelvinBareIsland(ib,ia)=temperature;
						}
					}
				}

			}
			growthStep++;			
			for(map<vector<double>, Matrix<float> >::iterator itr = spHabt.begin(); itr != spHabt.end(); ++itr) 
			{
				const vector<double>& parameters = itr->first; // key (== parameterSet)
				Matrix<float>& habtsuitability = itr->second; // value (== suitability)
				vector<double> parameterSet = parameters;
				habitatSuitability(parameterSet,spHabt);
			}
			
			//Updating habitat suitability mainland:
			for(map<vector<double>, Matrix<float> >::iterator itra = spHabtMain.begin(); itra != spHabtMain.end(); ++itra) 
			{
				const vector<double>& parameters = itra->first; // key (== parameterSet)
				Matrix<float>& habtsuitability = itra->second; // value (== abundances)
				bool contained = spHabt.find(parameters) != spHabt.end();
				if(contained) {
					habtsuitability=spHabt[parameters];
				}
				else{
					vector<double> parameterSet = parameters;
					habitatSuitability(parameterSet,spHabtMain);
				}
			}

			//Updating habitatPatches:
		    patchesPosition.clear();
			for (unsigned int r=0;r<NY;++r)
			{
				for (unsigned int c=0;c<NX;++c)
				{
					if(State(c,r)>=1) 
					{
						//if the cell content is 1, add the cell to the patches
						patchesPosition.push_back(std::pair<int,int>(c,r));
					}
				}
			}

			//vulcanic disturbance:
			if (vulcanicDisturbance==1)
			{
				int disturbanceRadius = rand () % growthStep;
				for (unsigned int x=centralx-disturbanceRadius;x<=centralx+disturbanceRadius;++x)
				{
					for (unsigned int y=centraly-disturbanceRadius;y<=centraly+disturbanceRadius;++y)
					{
						float disturbance = (rand() % 100)*(float)0.01;
						for(map<vector<double>, Matrix<float> >::iterator itr = spAdult.begin(); itr != spAdult.end(); ++itr) {
							const vector<double>& parameters = itr->first; // key (== parameterSet)
							
							float AreaTotalbeforeDistbrub =(spAdult[parameters](x,y)*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(x,y))))*allometabolConstTree[13]))) +spRecruit[parameters](x,y)*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*TempKelvin(x,y))))*allometabolConstTree[13]));
							spAdult[parameters](x,y) =spAdult[parameters](x,y) -spAdult[parameters](x,y)*disturbance;
							spRecruit[parameters](x,y) =spRecruit[parameters](x,y) -spRecruit[parameters](x,y)*disturbance;
							spSeed[parameters](x,y) =spSeed[parameters](x,y) -spSeed[parameters](x,y)*disturbance;							
							float AreaTotalafterDistbrub =(spAdult[parameters](x,y)*(maxK/((1/pow(parameters[16],0.75))*exp((E/(Boltzmann*TempKelvin(x,y))))*allometabolConstTree[13]))) +spRecruit[parameters](x,y)*(maxK/((1/pow(parameters[20],0.75))*exp((E/(Boltzmann*TempKelvin(x,y))))*allometabolConstTree[13]));
							// update K
							K(x,y)=K(x,y) + AreaTotalbeforeDistbrub - AreaTotalafterDistbrub;
							if(K(x,y)>maxK) K(x,y)=maxK;
							if(K(x,y)<0) K(x,y)=0;
						}
					}
				}
			}
			unsigned int InitialGammaRichnessV=spAdult.size();
			prespeciesExtinct=0;
			removeExtinctSpecies(spAdult,spRecruit, spSeed, spDisp,spHabt);
			VulcanicExt=InitialGammaRichnessV - (spAdult.size()+prespeciesExtinct); 

		}
	
		//Island erosion:
		if ((timeStep==(islandTimestep)*(erosionStep+erosionStep2+1))&&(erosionStep<=maxIslandRadius))
		{
			for (unsigned int r=0;r<NY;++r)
			{
				for (unsigned int c=0;c<NY;++c)
				{
					if(State(c,r)>=1) 
					{
						if (altitudinalTemp=="on"){
							TempKelvin(c,r)=TempKelvin(c,r)+1;
							TempKelvinBareIsland(r,c)= TempKelvinBareIsland(r,c)+1;
							}
					}
				}
			}
			bool eroded =false;
			for (unsigned int i=0; i<patchesPosition.size(); i++ )
			{
				unsigned int is = patchesPosition[i].first;
				unsigned int iz = patchesPosition[i].second;
				
				if((is==centralx)&&(iz==centraly)&&(eroded ==false)){
					for (unsigned int ib=(is-(maxIslandRadius-(erosionStep-1))); ib<(is+(maxIslandRadius-(erosionStep-2))); ib++ )
					{
						if(State(ib,iz-(maxIslandRadius-(erosionStep-1)))>=1)
							{
								State(ib,iz-(maxIslandRadius-(erosionStep-1)))=-99;
								K(ib,iz-(maxIslandRadius-(erosionStep-1)))=0;
							}
						if(State(ib,iz+(maxIslandRadius-(erosionStep-1)))>=1)
							{
								State(ib,iz+(maxIslandRadius-(erosionStep-1)))=-99;
								K(ib,iz+(maxIslandRadius-(erosionStep-1)))=0;
							}		
					}
					for (unsigned int ia=(iz-(maxIslandRadius-(erosionStep-1))); ia<(iz+(maxIslandRadius-(erosionStep-2))); ia++ )
					{
						if(State(is-(maxIslandRadius-(erosionStep-1)),ia)>=1)
							{
								State(is-(maxIslandRadius-(erosionStep-1)),ia)=-99;
								K(is-(maxIslandRadius-(erosionStep-1)),ia)=0;
							}	
						if(State(is+(maxIslandRadius-(erosionStep-1)),ia)>=1)
							{
								State(is+(maxIslandRadius-(erosionStep-1)),ia)=-99;
								K(is+(maxIslandRadius-(erosionStep-1)),ia)=0;
							}	
					}
					eroded=true;
				}
			}
			erosionStep++;
			erosionStep2++;
			patchesPosition.clear();
			for (unsigned int r=0;r<NY;++r)
			{
				for (unsigned int c=0;c<NX;++c)
				{
					if(State(c,r)>=1) 
					{
						//if the cell content is 1, add the cell to the patches
						patchesPosition.push_back(std::pair<int,int>(c,r));
					}
				}
			}
			
			//Updating Habitat suitability of island species
			for(map<vector<double>, Matrix<float> >::iterator itr = spAdult.begin(); itr != spAdult.end(); ++itr) {
				const vector<double>& parameters = itr->first; // key (== parameterSet)
				Matrix<float>& abundances = itr->second; // value (== abundances)
				for (unsigned int r=0;r<NY;++r)
				{
					for (unsigned int c=0;c<NY;++c)
					{
						if(State(c,r)<1) 
						{
							abundances(c,r)=0;
							spRecruit[parameters](c,r)=0;
							spSeed[parameters](c,r)=0;
						}
					}
				}				
				vector<double> parameterSet = parameters;
				habitatSuitability(parameterSet,spHabt);
			}
			
			//updating Mainland species pool:
			for(map<vector<double>, Matrix<float> >::iterator itra = spHabtMain.begin(); itra != spHabtMain.end(); ++itra) 
			{
				const vector<double>& parameters = itra->first; // key (== parameterSet)
				Matrix<float>& habtsuitability = itra->second; // value (== abundances)
				bool contained = spHabt.find(parameters) != spHabt.end();
				if(contained) {
					habtsuitability=spHabt[parameters];
				}
				else{					
					vector<double> parameterSet = parameters;
					habitatSuitability(parameterSet,spHabtMain);
				}
			}
			unsigned int InitialGammaRichnessB=spAdult.size();
			prespeciesExtinct=0;
			removeExtinctSpecies(spAdult,spRecruit, spSeed, spDisp,spHabt);
			ErosionExt= InitialGammaRichnessB - (spAdult.size()+prespeciesExtinct); 
		}
	}
}