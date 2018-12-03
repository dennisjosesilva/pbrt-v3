#include "plugin/multipole/multipole.hpp"
#include <vector>

int main(int argc, char const *argv[])
{
	std::vector<MultipoleLayer> layers = { MultipoleLayer{1.0f, 1.0f, 0.05, 0.35, 0.005}, 
		MultipoleLayer{1.0f, 1.0f, 0.05f, 0.35f, 0.005f}};
	
	MultipoleOptions options{0.05, 62};

	MultipoleTable table = ComputeMultipoleDiffusionProfile(layers, options);


	for (auto i = 0; i < table.NSamples(); i++)
		std::cout << table.squaredDistance(i) << "\n";

	std::cout << "table size: " << table.NSamples() << "\n";  
	std::cout <<  "DONE..." << "\n";

	return 0;
}