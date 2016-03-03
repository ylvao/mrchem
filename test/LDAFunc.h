/** 
 *
 * \date Jul 13, 2010
 * \author Stig Rune Jensen \n
 *		   CTCC, University of Tromsø
 *
 *
 */

#ifndef LDAFUNC_H
#define LDAFUNC_H

#include <boost/timer.hpp>

#include "RepresentableObject.h"
#include "MathUtils.h"
#include "FunctionTree.h"
#include "Filter.h"
#include "ProjectedNode.h"

template<int D> class NodeIndex;

template<int D>
class LDAFunc : public RepresentableObject<D> {
public:
	LDAFunc() { }
	~LDAFunc() { }

	FunctionTree<D> operator()(FunctionTree<D> &in) {
		FunctionTree<D> out = in.shallowCopy();
		this->density = &in;
		out.applyFunctional(*this);
		return out;
	}

	void apply(FunctionTree<D> &in, FunctionTree<D> &out) {
		boost::timer rolex;
		rolex.restart();
		println(0, "\nApplying LDA");
		this->density = &in;
		out.applyFunctional(*this);
		println(0, "Elapsed time:     " << rolex.elapsed());
	}

	void calcWaveletCoefs(MWNode<D> &node) {
		assert(density != 0);
		node.setCoefs(density->getNode(node.getNodeIndex()).getCoefs());
		node.mwTransform(Reconstruction);
		node.cvTransform(MWNode<D>::Forward);
		transformValues(node.getCoefs());
		node.cvTransform(MWNode<D>::Backward);
		node.mwTransform(Compression);
		node.setHasCoefs();
		node.calcNorms();
	}

	bool checkSeedNode(MWNode<D> &node) {
		assert(density != 0);
		return density->checkSeedNode(node);
	}

protected:
	FunctionTree<D> *density;

	virtual void transformValues(Eigen::VectorXd &val) = 0;
};

#endif // LDAFUNC_H
