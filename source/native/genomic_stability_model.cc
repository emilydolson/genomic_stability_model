//  This file is part of genomic_stability_model
//  Copyright (C) Emily Dolson, 2019.
//  Released under MIT license; see LICENSE

// This is the main function for the NATIVE version of this project.

#include <iostream>

#include "../genomic_stability_model.h"
#include "base/vector.h"
#include "config/command_line.h"

int main(int argc, char* argv[])
{
  InstabilityConfig config;
  auto args = emp::cl::ArgManager(argc, argv);
  if (args.ProcessConfigOptions(config, std::cout, "InstabilityConfig.cfg", "Instability-macros.h") == false) exit(0);
  if (args.TestUnknown() == false) exit(0);  // If there are leftover args, throw an error.

  // Write to screen how the experiment is configured
  std::cout << "==============================" << std::endl;
  std::cout << "|    How am I configured?    |" << std::endl;
  std::cout << "==============================" << std::endl;
  config.Write(std::cout);
  std::cout << "==============================\n" << std::endl;

  emp::Random rnd(config.SEED());

  InstabilityWorld world(rnd);
  world.Setup(config);

  world.Run();
}


