#include <catch2/catch_test_macros.hpp>

#include <ipc/potentials/normal_adhesion_potential.hpp>
#include <ipc/potentials/tangential_adhesion_potential.hpp>

using namespace ipc;

TEST_CASE("Normal adhesion potential", "[potential][adhesion]")
{
    NormalAdhesionPotential potential(1e-3, 2e-3, 1e3, 0.5);
}

TEST_CASE("Tangetial adhesion potential", "[potential][adhesion]")
{
    TangentialAdhesionPotential potential(1e-3);
}