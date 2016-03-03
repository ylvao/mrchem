/** Diff two trees on disk
 * \author Jonas Juselius
 */

//TOOD (short term): splitThrs, backTransform & GenNode

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>

#include "TelePrompter.h"
#include "constants.h"
#include "Plot.h"
#include "FunctionTree.h"
#include "parallel.h"

using namespace std;
namespace po = boost::program_options;

int main(int argc, char **argv) {
	mpi::environment env(argc, argv);
	mpi::communicator world;
	cout << scientific << setprecision(10);
	int rank = world.rank();
	int prank;

	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("rank,r", po::value<int>(&prank)->default_value(0), "rank for distributed trees")
		("verbose,v", "be verbose" )
		("file1", po::value<string>(), "file 1")
		("file2", po::value<string>(), "file 2")
	;
	po::positional_options_description p;
	p.add("file1", 1);
	p.add("file2", 1);

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

	stringstream f1, f2;
	if (vm.count("file1") == 1) {
		f1 << vm["file1"].as<string>();
	} else {
		MSG_ERROR("Missing file argument!");
		return -1;
	}

	if (vm.count("file2") == 1) {
		f2 << vm["file2"].as<string>();
	} else {
		MSG_ERROR("Missing file argument!");
		return -1;
	}

	if (world.size() > 1) {
		f1 << "." << rank;
		f2 << "." << rank;
	}

	FunctionTree<3> t1;
	FunctionTree<3> t2;

	t1.loadTree(f1.str());
	t2.loadTree(f2.str());
	bool diff = t1.diffTree(t2);
	if (diff) {
		println(1, "Trees differ.");
	}
	return diff;
}


