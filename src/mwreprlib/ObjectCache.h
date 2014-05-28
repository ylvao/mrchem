/*
 * 
 *
 *  \date Jul 25, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 *
 */

#ifndef OBJECTCACHE_H_
#define OBJECTCACHE_H_

#include <string>
#include <vector>
#include "TelePrompter.h"
#include "parallel.h"
#include "macros.h"

#define getObjectCache(T,X) \
	ObjectCache<T> &X = ObjectCache<T>::getInstance();

#ifdef OPENMP
#define SET_CACHE_LOCK() omp_set_lock(&this->cache_lock)
#define UNSET_CACHE_LOCK() omp_unset_lock(&this->cache_lock)
#define TEST_CACHE_LOCK() omp_test_lock(&this->cache_lock)
#else
#define SET_CACHE_LOCK()
#define UNSET_CACHE_LOCK()
#define TEST_CACHE_LOCK()
#endif

template<class T>
class ObjectCache {
public:
	static ObjectCache<T> &getInstance() {
		static ObjectCache<T> theObjectCache;
		return theObjectCache;
	}

	virtual void clear() {
		for (unsigned int i = 0; i < objs.size(); i++) {
			if (objs[i] != 0) {
				unload(i);
			}
		}
	}

	virtual void load(int id) {
		MSG_INFO("This routine does nothing in this class.")
	}

	void load(int id, T *new_o, int memory) {
		if (id >= highWaterMark) {
			for (int i = 0; i < id - highWaterMark + 1; i++) {
				objs.push_back(0);
				mem.push_back(0);
			}
			highWaterMark = id;
		}
		if (objs[id] != 0)
			return;
		mem[id] = memory;
		memLoaded += memory;
		objs[id] = new_o;

	}
	virtual void unload(int id) {
		if (id < 0 or id > highWaterMark) {
			THROW_ERROR("Id out of bounds:" << id)
		}
		if (objs[id] == 0) {
			MSG_WARN ("Object not loaded.")
			return;
		}
		memLoaded -= mem[id];
		mem[id] = 0;
		delete objs[id];
		objs[id] = 0;
	}
	virtual T &get(int id) {
		if (id < 0) {
			THROW_ERROR("Id out of bounds:" << id)
		}
		if (objs[id] == 0) {
			THROW_ERROR("Object not loaded!")
		}
		return *(objs[id]);
	}
	virtual T &operator[](int id) {
		return get(id);
	}
	bool hasId(int id) {
		if (id > highWaterMark)
			return false;
		if (objs[id] == 0)
			return false;
		return true;
	}

	int getNObjs() {
		return objs.size();
	}

	int getMem() {
		return memLoaded;
	}
	int getMem(int id) {
		return mem[id];
	}

protected:
	ObjectCache() {
		highWaterMark = 0;
		memLoaded = 0;
		objs.push_back(0);
		mem.push_back(0);
#ifdef OPENMP
		omp_init_lock(&cache_lock);
#endif
	}

	virtual ~ObjectCache() {
		SET_CACHE_LOCK();
		clear();
		UNSET_CACHE_LOCK();
#ifdef OPENMP
		omp_destroy_lock(&cache_lock);
#endif
	}
	ObjectCache(ObjectCache<T> const &oc) {
	}
	ObjectCache<T> &operator=(ObjectCache<T> const &oc) {
		return *this;
	}
#ifdef OPENMP
	omp_lock_t cache_lock;
#endif
private:
	int highWaterMark;
	int memLoaded; ///< memory occupied by loaded objects
	std::vector<T *> objs; ///< objects store
	std::vector<int> mem; ///< mem per object
};

#endif /* OBJECTCACHE_H_ */
