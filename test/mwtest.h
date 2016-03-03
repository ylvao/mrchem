/*
 * This file has common definitions and convenience variables to simplify
 * and shorten tests.
 */
#ifndef MWTEST_H
#define MWTEST_H

#include <gtest/gtest.h>
#include <boost/timer.hpp>
#include <Getkw.h>
#include "FunctionTree.h"
#include "parallel.h"
#include "MREnv.h"
#include "TelePrompter.h"
#include "config.h"

#define ELAPSED(X,S) if (rank == 0) { \
	println(1, "  @Timing for " << S << ": " << X.elapsed()); X.restart();}

double testThrs = 1.0e-9;

bool Debug;
Getkw Input;

#endif
