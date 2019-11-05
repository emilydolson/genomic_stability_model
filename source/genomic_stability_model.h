#ifndef _GENOMIC_STABILITY_H
#define _GENOMIC_STABILITY_H

#include "config/ArgManager.h"
#include "Evolve/World.h"


// TODO:
// - Implement non-spatial
// - Implement ability to choose mutation distribution

EMP_BUILD_CONFIG( InstabilityConfig,
  GROUP(MAIN, "Global settings"),
  VALUE(SEED, int, -1, "Random number generator seed"),
  VALUE(TIME_STEPS, int, 1000, "Number of time steps to run for"),
  VALUE(WORLD_X, int, 1000, "Length of plate in cell widths"),
  VALUE(WORLD_Y, int, 1000, "Width of plate in cell widths"),
  VALUE(INIT_POP_SIZE, int, 1, "Number of cells to seed population with"),
  VALUE(DATA_RESOLUTION, int, 10, "How many updates between printing data?"),
  VALUE(SPATIAL, bool, true, "Is the model spatial?"),
  VALUE(CELL_DEATH_PROB, double, .01, "Probability of stochastic cell death"),
  VALUE(MAX_CELLS, int, 20000, "Maximum number of cells to allow before we end model"),
  VALUE(MAX_FITNESS, double, 10.0, "Maximum fitness (for purposes of calculating reproduction probability)"),  
);

struct Cell {
    double fitness = .1;

    Cell() {;}

    bool operator< (const Cell & other) const {
      return fitness < other.fitness;
    }

};

class InstabilityWorld : public emp::World<Cell> {
  protected:
  int TIME_STEPS;
  int WORLD_X;
  int WORLD_Y;
  int INIT_POP_SIZE;
  int DATA_RESOLUTION;
  bool SPATIAL;
  double CELL_DEATH_PROB;
  int MAX_CELLS;
  int MAX_FITNESS;

  public:
  InstabilityWorld(emp::Random & r) : emp::World<Cell>(r) {;}
  InstabilityWorld() {;}

  ~InstabilityWorld() {;}

  void InitConfigs(InstabilityConfig & config) {
    TIME_STEPS = config.TIME_STEPS();
    INIT_POP_SIZE = config.INIT_POP_SIZE();
    WORLD_X = config.WORLD_X();
    WORLD_Y = config.WORLD_Y();
    DATA_RESOLUTION = config.DATA_RESOLUTION();
    SPATIAL = config.SPATIAL();
    MAX_CELLS = config.MAX_CELLS();
    MAX_FITNESS = config.MAX_FITNESS();    
  }

  size_t GetWorldX() {
    return WORLD_X;
  }

  size_t GetWorldY() {
    return WORLD_Y;
  }


  void InitPop() {
    // for (int cell_id = 0; cell_id < WORLD_X * WORLD_Y; cell_id++) {
    //   InjectAt(Cell(CELL_STATE::HEALTHY), cell_id);
    // }
    pop.resize(WORLD_X*WORLD_Y);
    size_t initial_spot = random_ptr->GetUInt(WORLD_Y*WORLD_X);
    InjectAt(Cell(), initial_spot);

    for (size_t cell_id = 1; cell_id < (size_t)INIT_POP_SIZE; cell_id++) {
      size_t spot = random_ptr->GetUInt(WORLD_Y*WORLD_X);
      AddOrgAt(emp::NewPtr<Cell>(), spot, initial_spot);
    }
  }

  void Reset(InstabilityConfig & config, bool web = false) {
    Clear();
    Setup(config, web);    
  }

  void Setup(InstabilityConfig & config, bool web = false) {
    InitConfigs(config);

    // emp::Ptr<emp::Systematics<Cell, int> > sys;
    // sys.New([](const Cell & c){return c.clade;});
    // AddSystematics(sys);

    // SetupFitnessFile().SetTimingRepeat(config.DATA_RESOLUTION());
    // SetupSystematicsFile().SetTimingRepeat(config.DATA_RESOLUTION());
    // SetupPopulationFile().SetTimingRepeat(config.DATA_RESOLUTION());

    // emp::DataFile & phylodiversity_file = SetupFile("phylodiversity.csv");
    // sys->AddEvolutionaryDistinctivenessDataNode();
    // sys->AddPairwiseDistanceDataNode();
    // sys->AddPhylogeneticDiversityDataNode();

    // phylodiversity_file.AddVar(update, "generation", "Generation");
    // phylodiversity_file.AddStats(*sys->GetDataNode("evolutionary_distinctiveness") , "evolutionary_distinctiveness", "evolutionary distinctiveness for a single update", true, true);
    // phylodiversity_file.AddStats(*sys->GetDataNode("pairwise_distance"), "pairwise_distance", "pairwise distance for a single update", true, true);
    // phylodiversity_file.AddCurrent(*sys->GetDataNode("phylogenetic_diversity"), "current_phylogenetic_diversity", "current phylogenetic diversity", true, true);
    // phylodiversity_file.PrintHeaderKeys();
    // phylodiversity_file.SetTimingRepeat(config.DATA_RESOLUTION());

    // emp::DataFile & phylodiversity_file = SetupFile("phylodiversity.csv");

    SetPopStruct_Grid(WORLD_X, WORLD_Y, true);
    InitPop();

    // SetSynchronousSystematics(true);
  }


  /// Determine if cell can divide (i.e. is space available). If yes, return
  /// id of cell that it can divide into. If not, return -1.
  int CanDivide(size_t cell_id) {
    emp::vector<int> open_spots;
    int x_coord = (int)(cell_id % WORLD_X);
    int y_coord = (int)(cell_id / WORLD_X);
    
    // Iterate over 9-cell neighborhood. Currently checks focal cell uneccesarily,
    // but that shouldn't cause problems because it will never show up as invasible.
    for (int x = std::max(0, x_coord-1); x < std::min((int)WORLD_X, x_coord + 2); x++) {
      for (int y = std::max(0, y_coord-1); y < std::min((int)WORLD_Y, y_coord + 2); y++) {
        int this_cell = y*(int)WORLD_X + x;
        // Cells can be divided into if they are empty or if they are healthy and the
        // dividing cell is cancerous
        if (!IsOccupied((size_t)this_cell)) {
          open_spots.push_back(this_cell);
        }
      }
    }
    
    // -1 is a sentinel value indicating no spots are available
    if (open_spots.size() == 0) {
      return -1;
    }

    // If there are one or more available spaces, return a random spot
    return open_spots[random_ptr->GetUInt(0, open_spots.size())];
  }

  int Mutate(emp::Ptr<Cell> c){
    c->fitness += random_ptr->GetRandNormal(0, 1);

    return 0;
  }


  void Quiesce(size_t cell_id) {
    // Quiescence - stick the cell back into the population in
    // the same spot but don't change anything else
    emp::Ptr<Cell> cell = emp::NewPtr<Cell>(*pop[cell_id]);
    AddOrgAt(cell, emp::WorldPosition(cell_id,1), cell_id);      
  }

  void RunStep() {
    std::cout << update << std::endl;

    for (int cell_id = 0; cell_id < WORLD_X * WORLD_Y; cell_id++) {
      if (!IsOccupied(cell_id)) {
        // Don't need to do anything for dead/empty cells
        continue;
      }

      //   size_t x = cell_id % WORLD_X;
      //   size_t y = cell_id / WORLD_X;

      // Die with specified probability
      if (random_ptr->P(CELL_DEATH_PROB)) {
        continue; // Don't add cell to next population (i.e. it died)
      }

      // Check for space for division
      int potential_offspring_cell = CanDivide(cell_id);

      // If space and cell gets sufficiently lucky, divide
      double divide_prob = pop[cell_id]->fitness/MAX_FITNESS;
      if (divide_prob > 1){
          divide_prob = 1;
      }
      if (potential_offspring_cell != -1 && random_ptr->P(divide_prob)) {
        // Cell divides

        // Handle daughter cell in previously empty spot
        before_repro_sig.Trigger(cell_id);
        emp::Ptr<Cell> offspring = emp::NewPtr<Cell>(*pop[cell_id]);
        Mutate(offspring);
        offspring_ready_sig.Trigger(*offspring, cell_id);
        AddOrgAt(offspring, emp::WorldPosition((size_t)potential_offspring_cell, 1), cell_id);

        // Handle daughter cell in current location
        before_repro_sig.Trigger(cell_id);
        offspring = emp::NewPtr<Cell>(*pop[cell_id]);
        Mutate(offspring);
        offspring_ready_sig.Trigger(*offspring, cell_id);
        AddOrgAt(offspring, emp::WorldPosition(cell_id,1), cell_id);

      } else {        
        Quiesce(cell_id);
      }
    }

    Update();
  }

  void Run() {
      for (int u = 0; u <= TIME_STEPS; u++) {
          RunStep();
          // Check end conditions
          if (pop.size() == 0 || pop.size() > (size_t) MAX_CELLS) {
            // We either killed everything, or population got too large
            break;
          }
      }
    //   systematics[0].DynamicCast<emp::Systematics<Cell, int>>()->Snapshot("memic_phylo.csv");  
  }

};

#endif