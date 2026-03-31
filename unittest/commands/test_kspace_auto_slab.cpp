/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing authors: Jiuyang Liang and Xuanzhao Gao (Flatiron Institute)
------------------------------------------------------------------------- */

#include "domain.h"
#include "force.h"
#include "input.h"
#include "kspace.h"
#include "math_const.h"
#include "utils.h"

#include "../testing/core.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <mpi.h>
#include <vector>

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

bool verbose = false;

namespace LAMMPS_NS {
using namespace MathConst;

class KSpaceAutoSlabTest : public LAMMPSTest {
 protected:
  void InitSystem() override
  {
    if (!info->has_style("atom", "full")) GTEST_SKIP();
    if (!info->has_style("pair", "coul/long")) GTEST_SKIP();

    command("boundary p p f");
    command("variable input_dir index \"" STRINGIFY(TEST_INPUT_FOLDER) "\"");
    command("include \"${input_dir}/in.fourmol\"");

    command("pair_style coul/long 8.0");
    command("pair_coeff * *");
    command("pair_modify table 0");
  }

  double expected_volfactor(double force_accuracy, double alpha) const
  {
    const double xprd = lmp->domain->xprd;
    const double yprd = lmp->domain->yprd;
    const double zprd = lmp->domain->zprd;
    const double logeps = std::log(1.0 / force_accuracy);
    const double lateral = std::max(xprd, yprd) * logeps / MY_2PI;
    const double reciprocal = std::sqrt(logeps) / alpha;
    return std::max((zprd + std::max(lateral, reciprocal)) / zprd, 1.0);
  }
};

TEST_F(KSpaceAutoSlabTest, PPPMUsesPaperFormulaWithManualGewald)
{
  if (!info->has_style("kspace", "pppm")) GTEST_SKIP();

  HIDE_OUTPUT([&] {
    command("kspace_style pppm 1.0e-6");
    command("kspace_modify gewald 0.3");
    command("kspace_modify slab auto");
    command("run 0 post no");
  });

  ASSERT_NE(lmp->force->kspace, nullptr);
  EXPECT_NEAR(lmp->force->kspace->slab_volfactor, expected_volfactor(lmp->force->kspace->accuracy, 0.3), 1.0e-12);
}

TEST_F(KSpaceAutoSlabTest, EwaldUsesPaperFormulaWithManualGewald)
{
  if (!info->has_style("kspace", "ewald")) GTEST_SKIP();

  HIDE_OUTPUT([&] {
    command("pair_modify mix arithmetic");
    command("kspace_style ewald 1.0e-6");
    command("kspace_modify gewald 0.3");
    command("kspace_modify slab auto");
    command("run 0 post no");
  });

  ASSERT_NE(lmp->force->kspace, nullptr);
  EXPECT_NEAR(lmp->force->kspace->slab_volfactor, expected_volfactor(lmp->force->kspace->accuracy, 0.3), 1.0e-12);
}

}    // namespace LAMMPS_NS

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleMock(&argc, argv);

  if (const char *var = getenv("TEST_ARGS")) {
    std::vector<std::string> env = LAMMPS_NS::utils::split_words(var);
    for (auto &arg : env) {
      if (arg == "-v") verbose = true;
    }
  }

  if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

  const int rv = RUN_ALL_TESTS();
  MPI_Finalize();
  return rv;
}
