//  This file is part of Project Name
//  Released under the MIT Software license; see doc/LICENSE

#include "web/web.h"
#include "web/color_map.h"
#include "../genomic_stability_model.h"
#include "config/config_web_interface.h"
#include "tools/spatial_stats.h"

namespace UI = emp::web;

class InstabilityWebInterface : public UI::Animate, public InstabilityWorld{
  // friend class InstabilityWorld;
  // friend class UI::Animate;

  using color_fun_t = std::function<std::string(int)>;
  using should_draw_fun_t = std::function<bool(int)>;

  InstabilityConfig config;
  emp::ConfigWebUI config_ui;
  emp::Random r;
  // InstabilityWorld world;
  // UI::Animate anim;
  UI::Document cell_area;
  UI::Document controls;
  UI::Document stats_area;
  UI::Canvas cell_display;
  // UI::Canvas clade_display;
  const double display_cell_size = 10;

  UI::Button toggle;
  UI::Style button_style;
  bool draw_cells = true;


  color_fun_t cell_color_fun;
  UI::Selector cell_color_control;
  color_fun_t fitness_color_fun = [this](int cell_id) {
                                        double fitness = pop[cell_id]->fitness/MAX_FITNESS;
                                        if (fitness > 1) {
                                          fitness = 1;
                                        } else if (fitness < 0) {
                                          return emp::ColorHSL(0,100,100);
                                        }
                                        double hue =  fitness * 280.0;
                                        // std::cout << pop[cell_id]->fitness << " " << hue << std::endl;
                                        return emp::ColorHSL(hue,50,50);
                                     };

  color_fun_t stability_color_fun = [this](int cell_id) {
                                        double stability = pop[cell_id]->stability/10;
                                        if (stability > 1) {
                                          stability = 1;
                                        } 
                                        double hue = stability * 280.0;
                                        // std::cout << pop[cell_id]->fitness << " " << hue << std::endl;
                                        return emp::ColorHSL(hue,50,50);
                                     };

  should_draw_fun_t should_draw_cell_fun;
  should_draw_fun_t draw_if_occupied = [this](int cell_id){return IsOccupied(cell_id);};
  should_draw_fun_t always_draw = [](int cell_id){return true;};

  public:
  InstabilityWebInterface() : config_ui(config), cell_area("cell_area"), controls("control_area"), stats_area("stats_area")
    , cell_display(100, 100, "cell_display"),
    cell_color_control("cell_color_control")
    // : anim([this](){DoFrame();}, oxygen_display, cell_display) 
  {
    SetupInterface();   
  }

  void SetupInterface() {
    GetRandom().ResetSeed(config.SEED());
    Setup(config, true);

    cell_area.SetWidth(WORLD_X * display_cell_size);
    cell_display.SetSize(WORLD_X * display_cell_size, WORLD_Y * display_cell_size);
    cell_display.Clear("black");

    cell_area << "<h1 class='text-center'>Cells</h1>" << cell_display;
    controls << "<h1 class='text-center'>Controls</h1>";
    stats_area << "<h1 class='text-center'>Statistics</h1>";

    cell_color_fun = fitness_color_fun;
    should_draw_cell_fun = draw_if_occupied;


    cell_color_control.SetOption("Fitness", 
                                 [this](){
                                     cell_color_fun = fitness_color_fun;
                                     should_draw_cell_fun = draw_if_occupied;
                                     RedrawCells();
                                 }, 0);

    cell_color_control.SetOption("Stability", 
                                 [this](){
                                     cell_color_fun = stability_color_fun;
                                     should_draw_cell_fun = draw_if_occupied;
                                     RedrawCells();
                                 }, 1);


    toggle = GetToggleButton("but_toggle");
    button_style.AddClass("btn");
    button_style.AddClass("btn-primary");
    toggle.SetCSS(button_style);

    UI::Button reset_button([this](){Reset(config, true); RedrawCells();}, "Reset");
    reset_button.SetCSS(button_style);

    controls << toggle;
    controls << " " << reset_button << " " << cell_color_control << "<br>";

    RedrawCells();
    
    emp::web::OnDocumentReady([](){
        // EM_ASM(d3.select("#but_toggle").classed("btn btn-primary", true););
        EM_ASM($('select').selectpicker('setStyle', 'btn-primary'););
    });

    config_ui.SetOnChangeFun([this](const std::string & val){ std::cout << "New val: " << val<<std::endl;;InitConfigs(config);});
    config_ui.ExcludeConfig("SEED");
    config_ui.ExcludeConfig("TIME_STEPS");
    config_ui.ExcludeConfig("TREATMENT_START");
    config_ui.ExcludeConfig("TREATMENT_MUT_PROB");
    config_ui.ExcludeConfig("DATA_RESOLUTION");
    config_ui.Setup();
    controls << config_ui.GetDiv();

    stats_area << "<br>Time step: " << emp::web::Live( [this](){ return GetUpdate(); } );
    stats_area << "<br>Mean fitness: " << emp::web::Live( [this](){ return GetFitnessDataNode()->GetMean(); } );
    stats_area << "<br>Variance of fitness: " << emp::web::Live( [this](){ return GetFitnessDataNode()->GetVariance(); } );
    stats_area << "<br>Max fitness: " << emp::web::Live( [this](){ return GetFitnessDataNode()->GetMax(); } );
    stats_area << "<br>Min fitness: " << emp::web::Live( [this](){ return GetFitnessDataNode()->GetMin(); } );
    // stats_area << "<br>Extant taxa: " << emp::web::Live( [this](){ return systematics[0].DynamicCast<emp::Systematics<Cell, int>>()->GetNumActive(); } );
    // stats_area << "<br>Shannon diversity: " << emp::web::Live( [this](){ return systematics[0].DynamicCast<emp::Systematics<Cell, int>>()->CalcDiversity(); } );
  }

  void DoFrame() {
    // std::cout << frame_count << " " << GetStepTime() << std::endl;
    RunStep();
    RedrawCells();
    stats_area.Redraw();
  }

  void RedrawCells(){
    cell_display.SetSize(WORLD_X * display_cell_size, WORLD_Y * display_cell_size);
    cell_display.Freeze();
    cell_display.Clear("black");

    for (size_t x = 0; x < WORLD_X; x++) {
      for (size_t y = 0; y < WORLD_Y; y++) {
        size_t cell_id = x + y * WORLD_X;
        if (should_draw_cell_fun(cell_id)) {
          std::string color = cell_color_fun(cell_id);


          cell_display.Rect(x*display_cell_size, y*display_cell_size, display_cell_size, display_cell_size, color, color);

        }        
      }
    }

    cell_display.Activate();
  }

};


InstabilityWebInterface interface;

int main()
{
}
