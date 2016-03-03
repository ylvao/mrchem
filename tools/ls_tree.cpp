/** Diff two trees on disk
 * \author Jonas Juselius
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <boost/program_options.hpp>

#include "TelePrompter.h"
#include "constants.h"
#include "Plot.h"
#include "FunctionTree.h"
#include "MWNode.h"
#include "parallel.h"

namespace po = boost::program_options;
using namespace std;


int main(int argc, char **argv) {
	mpi::environment env(argc, argv);
	mpi::communicator world;
	int rank = world.rank();
	cout << scientific << setprecision(10);

	int prank;

	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("rank,r", po::value<int>(&prank)->default_value(0), "rank for distributed trees")
		("verbose,v", "be verbose" )
		("file,f", po::value<string>(), "input file name")
	;
	po::positional_options_description p;
	p.add("file", 1);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).
			  options(desc).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("help")) {
		cout << desc << "\n";
		return 1;
	}

	int pl = vm.count("verbose") + 1;
	SET_PRINT_LEVEL(pl);

	stringstream f1;
	if (vm.count("file") == 1) {
		f1 << vm["file"].as<string>();
	} else {
		MSG_ERROR("Missing file argument!");
		return -1;
	}

	if (world.size() > 1) {
		f1 << "." << rank;
	}

	vector<vector<MWNode<3> *> > table;
	FunctionTree<3> t1;

	println(1, "Open file: " << f1.str());
	if (not t1.loadTree(f1.str())) {
		println(0, "File error: " << f1.str());
		return 1;
	}
	t1.makeNodeTable(table);

	for (unsigned int i = 0; i < table.size(); i++) {
		for (unsigned int j = 0; j < table[i].size(); j++) {
			if (world.size() > 1) {
				if (rank == prank)  {
					println(0, "$" << rank << ": " << *table[i][j]);
				}
			} else {
				println(0, *table[i][j]);
			}
		}
	}
}


