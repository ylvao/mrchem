#pragma once

#include <nlohmann/json.hpp>

namespace mrchem {
namespace mrenv {

nlohmann::json fetch_input(int argc, char **argv);
void initialize(const nlohmann::json &input);
void finalize(double wt);
void init_mpi(const nlohmann::json &json_mpi);

} // namespace mrenv
} // namespace mrchem
