/*
 * 
 *
 *  \date Jul 2, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#ifndef MULTIEXCEPTION_H_
#define MULTIEXCEPTION_H_

#include <cstdlib>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>

//TODO: define ABORT_EXCEPTION and redefine THROW not to throw, but abort.

#define THROW_ERROR(X) {std::ostringstream _err; \
	_err << "Error: " << __func__ << ",  line " << __LINE__ << \
			" in  " << __FILE__ << ": " << X << std::endl; \
	throw MultiException(_err); }

class MultiException: public std::exception {
public:
	MultiException();
	MultiException(const std::string &err);
	MultiException(std::ostringstream &err);
	virtual ~MultiException() throw();
	void trigger(const std::string &msg);
	static void setVerbose(bool flag);
	static void setStrict(bool flag);
	friend std::ostream& operator<<(std::ostream &o, const MultiException &e) {
			return o << e.msg;
	}
//	virtual const char *what() {
//		return err.c_str();
//	}
private:
	std::string msg;
	static bool verbose;
	static bool strict;
};

#endif /* MULTIEXCEPTION_H_ */
