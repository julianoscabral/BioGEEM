//---------------------------------------------------------------------------

#pragma hdrstop
#include <iostream>
#include <string>
#include<exception>
#include "time.h"
#include "Simulation.h"
#include <fstream>
#include <sstream>
#include <exception>
#include <math.h>
#include <time.h>
#include "ConfigFile.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

//---------------------------------------------------------------------------
void test(const ConfigFile& config) {

	//Create simulation
	Simulation simulation(config);
	simulation.perform();
}

int main(int argc, char* argv[])
{
	try {
		// initialise random generator
	    srand((unsigned int)time(NULL));
	    
		string configurationFileName = "Configuration_simplified.txt";

		// if a configuration file was given as parameter, use it
		// else, use default "Configuration.txt"
		if(argc == 2) {
			configurationFileName = argv[1];
		}

		ConfigFile config(configurationFileName);

		string simulationType = config.read<string>("simulationType");

		unsigned int begin = clock();	

		if (simulationType == "Test") {
			test(configurationFileName);
		}
		if (simulationType == "Default") {
			Simulation simulation(config);
	        simulation.perform();
		}
		else {
			throw std::exception((string("Unknown Simulation Type ") + simulationType).c_str());
		}

		//Simulation simulation(std::string("Configuration.txt"));

		unsigned int end = clock();
		cout << "Simulation Done (" << end-begin << " ms)" << std::endl;
	}
	catch(std::exception& e) {
		cout << e.what();
	}

	return 0;
}
//---------------------------------------------------------------------------
