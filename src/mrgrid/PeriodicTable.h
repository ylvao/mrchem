/** 
 *
 * \date Jun 7, 2009
 * \author Jonas Juselius <jonas.juselius@uit.no> \n
 *         CTCC, University of Tromsø
 *
 *
 */

#ifndef PERIODICTABLE_H_
#define PERIODICTABLE_H_

#include <string>
#include <map>

#include "AtomicElement.h"

class PeriodicTable {
public:
	PeriodicTable() { }
	virtual ~PeriodicTable() { }
	const AtomicElement &getAtomicElement(int Z) const;
	const AtomicElement &getAtomicElement(const char *id) const;
	static const int nAtomicElements=112;
protected:
	typedef std::map<std::string, const AtomicElement *> map_t;
	static const AtomicElement elements[];
	static map_t byName; // TODO should be const
	static map_t bySymbol; // TODO should be const
private:
	static map_t _init_byname();
	static map_t _init_bysymbol();
};

#endif /* PERIODICTABLE_H_ */
