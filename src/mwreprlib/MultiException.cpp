/*
 * 
 *
 *  \date Jul 2, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include "MultiException.h"

using namespace std;

bool MultiException::strict = false;
bool MultiException::verbose = true;

MultiException::MultiException() {
}

MultiException::MultiException(const string &err):
	msg(err) {
	if (verbose or strict) {
		cout << "Error: " << msg << endl;
	}
	if (strict) {
		cout << "Exiting..." << endl;
		exit(1);
	}
}

MultiException::MultiException(ostringstream &err) {
	this->msg = err.str();
	if (verbose or strict) {
		cout << "Error: " << msg << endl;
	}
	if (strict) {
		cout << "Exiting..." << endl;
		exit(1);
	}
}
MultiException::~MultiException() throw() {
	// TODO Auto-generated destructor stub
}

void MultiException::trigger(const string &err)
{
	msg = err;
	if (verbose or strict) {
		cout << "Error: " << msg << endl;
	}
	if (strict) {
		cout << "Exiting..." << endl;
		exit(1);
	}
	throw *this;
}

void MultiException::setVerbose(bool flag)
{
	verbose = flag;
}

void MultiException::setStrict(bool flag)
{
	strict = flag;
}
