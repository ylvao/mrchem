/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2019 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#include <iomanip>
#include <iostream>

#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "print_utils.h"

#include "qmfunctions/QMFunction.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

void print_utils::qmfunction(int level, const std::string &txt, const QMFunction &func, mrcpp::Timer &timer) {
    auto nodes = func.getNNodes(NUMBER::Total);
    auto memory = func.getSizeNodes(NUMBER::Total);
    auto time = timer.elapsed();
    mrcpp::print::tree(level, txt, nodes, memory, time);
}

void print_utils::coord(int level, const std::string &txt, const mrcpp::Coord<3> &val, int p, bool s) {
    std::stringstream o;
    o << " ";
    for (int i = 0; i < txt.size(); i++) o << txt[i];
    for (int i = txt.size(); i < 19; i++) o << " ";
    if (s) {
        for (auto r : val) o << std::setw(13) << std::setprecision(p) << std::scientific << r;
    } else {
        for (auto r : val) o << std::setw(13) << std::setprecision(p) << std::fixed << r;
    }
    println(level, o.str());
}

void print_utils::scalar(int level, const std::string &txt, const std::string &unit, double val, int p, bool s) {
    std::stringstream o;
    o << " ";
    for (int i = 0; i < txt.size(); i++) o << txt[i];
    for (int i = txt.size(); i < 30; i++) o << " ";
    for (int i = 0; i < unit.size(); i++) o << unit[i];
    for (int i = unit.size(); i < 8; i++) o << " ";
    if (s) {
        o << std::setw(20) << std::setprecision(p) << std::scientific << val;
    } else {
        o << std::setw(20) << std::setprecision(p) << std::fixed << val;
    }
    println(level, o.str());
}

void print_utils::vector(int level, const std::string &txt, const DoubleVector &val, int p, bool s) {
    std::stringstream o;
    o << " ";
    for (int i = 0; i < txt.size(); i++) o << txt[i];
    for (int i = txt.size(); i < 19; i++) o << " ";
    for (int i = 0; i < val.size(); i++) {
        if (s) {
            o << std::setw(13) << std::setprecision(p) << std::scientific << val(i);
        } else {
            o << std::setw(13) << std::setprecision(p) << std::fixed << val(i);
        }
    }
    println(level, o.str());
}

void print_utils::matrix(int level, const std::string &txt, const DoubleMatrix &val, int p, bool s) {
    std::stringstream o;
    for (int i = 0; i < val.rows(); i++) {
        o << " ";
        if (i == 0) {
            for (int i = 0; i < txt.size(); i++) o << txt[i];
            for (int i = txt.size(); i < 19; i++) o << " ";
        } else {
            for (int i = 0; i < 19; i++) o << " ";
        }
        for (int j = 0; j < val.cols(); j++) {
            if (s) {
                o << std::setw(13) << std::setprecision(p) << std::scientific << val(i, j);
            } else {
                o << std::setw(13) << std::setprecision(p) << std::fixed << val(i, j);
            }
        }
        o << std::endl;
    }
    println(level, o.str());
}

} // namespace mrchem
