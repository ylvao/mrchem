/*
 *
 *
 *  \date Oct 7, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif Simple storage class for scale and translation indexes.
 * The usefulness of the class becomes evident when examining
 * the parallel algorithms for projection & friends.
 */

#ifndef NODEINDEX_H_
#define NODEINDEX_H_

#include <boost/serialization/serialization.hpp>
#include <iostream>

#include "TelePrompter.h"

template<int D> class NodeIndexComp;

template<int D>
class NodeIndex {
public:
	NodeIndex(int n = 0, const int *l = 0) {
		N = (short int) n;
		setTranslation(l);
		rankId = 0;
	}
	virtual ~NodeIndex() {
	}
	NodeIndex(const NodeIndex<D> &idx) {
		N = idx.N;
		rankId = idx.rankId;
		setTranslation(idx.L);
	}
	NodeIndex<D> &operator=(const NodeIndex<D> &idx) {
		if (&idx == this) {
			return *this;
		}
		N = idx.N;
		rankId = idx.rankId;
		setTranslation(idx.L);
		return *this;
	}

	void setTranslation(const int *l) {
		for (int i = 0; i < D; i++) {
			if (l == 0) {
				L[i] = 0;
			} else {
				L[i] = l[i];
			}
		}
	}
	void setScale(const int n) {
		N = (short int) n;
	}
	void setRankId(const int n) {
		rankId = (unsigned short int) n;
	}

	const int *getTranslation() const {
		return L;
	}
	int getScale() const {
		return N;
	}
	int scale() const {
		return N;
	}
	int getRankId() const {
		return rankId;
	}
	int transl(int i) {
		assert(i >= 0 or i < D);
		return L[i];
	}
	const int &operator[](int i) const {
		assert(i >= 0 or i < D);
		return L[i];
	}
	int &operator[](int i) {
		assert(i >= 0 or i < D);
		return L[i];
	}
	bool operator==(const NodeIndex<D> &idx) const {
		if (N != idx.N) {
			return false;
		}
		for (int i = 0; i < D; i++) {
			if (L[i] != idx.L[i]) {
				return false;
			}
		}
		return true;
	}
	bool operator!=(const NodeIndex<D> &idx) const {
		if (*this == idx) {
			return false;
		}
		return true;
	}
	bool operator>=(const NodeIndex<D> &idx) const {
		for (int d = 0; d < D; d++) {
			if (this->L[d] < idx.getTranslation()[d]) {
				return false;
			}
		}
		return true;
	}
	bool operator<=(const NodeIndex<D> &idx) const {
		for (int d = 0; d < D; d++) {
			if (this->L[d] > idx.getTranslation()[d]) {
				return false;
			}
		}
		return true;
	}
	bool operator>(const NodeIndex<D> &idx) const {
		for (int d = 0; d < D; d++) {
			if (this->L[d] <= idx.getTranslation()[d]) {
				return false;
			}
		}
		return true;
	}
	bool operator<(const NodeIndex<D> &idx) const {
		for (int d = 0; d < D; d++) {
			if (this->L[d] >= idx.getTranslation()[d]) {
				return false;
			}
		}
		return true;
	}

	// Integer comparsion
	bool operator>=(int l) const {
		for (int d = 0; d < D; d++) {
			if (this->L[d] < l) {
				return false;
			}
		}
		return true;
	}
	bool operator<=(int l) const {
		for (int d = 0; d < D; d++) {
			if (this->L[d] > l) {
				return false;
			}
		}
		return true;
	}
	bool operator>(int l) const {
		for (int d = 0; d < D; d++) {
			if (this->L[d] <= l) {
				return false;
			}
		}
		return true;
	}
	bool operator<(int l) const {
		for (int d = 0; d < D; d++) {
			if (this->L[d] >= l) {
				return false;
			}
		}
		return true;
	}

	friend std::ostream& operator<<(std::ostream &o, const NodeIndex<D> &idx) {
		o << "[ " << idx.N << " | ";
		for (int i = 0; i < D - 1; i++) {
			o << idx.L[i] << ", ";
		}
		o << idx.L[D - 1] << "] @" << idx.rankId;
		return o;
	}
	friend class NodeIndexComp<D>;
private:
	short int N;
	int L[D];
	unsigned short int rankId;

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & N;
		ar & L;
		ar & rankId;
	}
};

template<int D>
class NodeIndexComp {
public:
	bool operator()(const NodeIndex<D> &a, const NodeIndex<D> &b) const {
		if (a->rankId < b->rankId) {
			return true;
		}
		if (a->rankId > b->rankId) {
			return false;
		}
		if (a.N < b.N) {
			return true;
		}
		if (a.N > b.N) {
			return false;
		}
		if (a.N == b.N) {
			for (int i = 0; i < D; i++) {
				if (a.L[i] == b.L[i]) {
					continue;
				}
				if (a.L[i] < b.L[i]) {
					return true;
				}
				return false;
			}
		}
		return false;
	}
	bool operator()(const NodeIndex<D> *a, const NodeIndex<D> *b) const {
		if (a->rankId < b->rankId) {
			return true;
		}
		if (a->rankId > b->rankId) {
			return false;
		}
		if (a->N < b->N) {
			return true;
		}
		if (a->N > b->N) {
			return false;
		}
		if (a->N == b->N) {
			for (int i = 0; i < D; i++) {
				if (a->L[i] == b->L[i]) {
					continue;
				}
				if (a->L[i] < b->L[i]) {
					return true;
				}
				return false;
			}
		}
		return false;
	}
};

#endif /* NODEINDEX_H_ */
