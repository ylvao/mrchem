/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRCPP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRCPP.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRCPP, see:
 * <https://mrcpp.readthedocs.io/>
 */

#pragma once

namespace mrcpp {
inline constexpr auto program_version() noexcept {
  return "1.6.0-alpha";
}

inline constexpr auto git_commit_hash() noexcept {
  return "720133372c9717134c5a";
}

inline constexpr auto git_commit_author() noexcept {
  return "gitpeterwind";
}

inline constexpr auto git_commit_date() noexcept {
  return "Thu Feb 15 13:44:31 2024 +0100";
}

inline constexpr auto git_branch() noexcept {
  return "HEAD";
}
} // namespace mrcpp
