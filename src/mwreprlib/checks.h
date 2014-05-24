/*
 * 
 *
 *  \date Jul 8, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#ifndef CHECKS_H_
#define CHECKS_H_

#define out_of_bounds(X,Y,Z)   ((X) < (Y) || (X) > (Z))
#define same_node_size(X,Y)  ( ( ((X)->dim == (Y)->dim) && ((X)->order == (Y)->order) ) ? (1):(-1) )


#endif /* CHECKS_H_ */
