/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#pragma once

namespace mrchem {
inline constexpr auto program_version() noexcept {
  return "1.2.0-alpha";
}

inline constexpr auto git_commit_hash() noexcept {
  return "ae0973689d691a075f6d-dirty";
}

inline constexpr auto git_commit_author() noexcept {
  return "Luca Frediani";
}

inline constexpr auto git_commit_date() noexcept {
  return "Wed Nov 27 13:21:59 2024 +0100";
}

inline constexpr auto git_branch() noexcept {
  return "master";
}
} // namespace mrchem
