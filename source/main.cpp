// Includes from PrEW
#include "GlobalVar/Chiral.h"
#include "Output/Printer.h"

// Includes from PrEWUtils
#include "Runners/ParallelRunner.h"
#include "SetupHelp/AccBoxInfo.h"
#include "SetupHelp/AccBoxPolynomialInfo.h"
#include "SetupHelp/ConstEffInfo.h"
#include "SetupHelp/CrossSectionInfo.h"
#include "SetupHelp/RunInfo.h"
#include "SetupHelp/TGCInfo.h"
#include "Setups/FitModifier.h"
#include "Setups/GeneralSetup.h"

// Standard library
#include <string>

#include "spdlog/spdlog.h"


int main (int /*argc*/, char **/*argv*/) {
  spdlog::set_level(spdlog::level::info);
  
  int energy = 250;
  int n_threads = 2;
  int n_toys = 2;
  std::string minuit_minimizers = "Combined(1000000,1000000,0.001)";
  std::string prew_minimizer = "PoissonNLL";
  std::string output_path = "../output/fit_results.out";
    
  spdlog::info("Start test.");
  
  spdlog::info("Create setup.");
  PrEWUtils::Setups::GeneralSetup setup(energy);

  spdlog::info("Add files.");
  setup.add_input_files("/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/4f_WW_sl/PrEWInput", ".*\\.csv", "CSV");
  setup.add_input_files("/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/4f_sW_sl/PrEWInput", ".*\\.csv", "CSV");
  setup.add_input_files("/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/2f_Z_h/PrEWInput", ".*\\.csv", "CSV");
  setup.add_input_files("/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/2f_Z_l/PrEWInput", ".*\\.csv", "CSV");
  

  spdlog::info("Selecting distributions.");
  setup.use_distr("SingleW_eminus");
  setup.use_distr("SingleW_eplus");
  setup.use_distr("WW_muminus");
  setup.use_distr("WW_muplus");
  setup.use_distr("2f_uds_81to101");
  setup.use_distr("2f_uds_180to275");
  setup.use_distr("2f_c_81to101");
  setup.use_distr("2f_c_180to275");
  setup.use_distr("2f_b_81to101");
  setup.use_distr("2f_b_180to275");
  setup.use_distr("2f_mu_81to101");
  setup.use_distr("2f_mu_180to275");
  
  
  spdlog::info("Setting up linking info.");
  
  spdlog::info("Creating run info.");
  PrEWUtils::SetupHelp::RunInfo run(energy);

  // Set up the systematics
  run.set_lumi(2000);
  run.add_pol("ePol-", 0.8);
  run.add_pol("ePol+", 0.8);
  run.add_pol("pPol-", 0.3);
  run.add_pol("pPol+", 0.3);
  
  // Constraints on the systematics
  run.add_lumi_constr(2000, 2000*3e-3);
  run.add_pol_constr("ePol-", 0.8, 0.8*2.5e-3);
  run.add_pol_constr("ePol+", 0.8, 0.8*2.5e-3);
  run.add_pol_constr("pPol-", 0.3, 0.3*2.5e-3);
  run.add_pol_constr("pPol+", 0.3, 0.3*2.5e-3);
  
  // Set up the desired polarisation sharings
  run.add_pol_config("e-p+", "ePol-", "pPol+", "-", "+", 0.45);
  run.add_pol_config("e+p-", "ePol+", "pPol-", "+", "-", 0.45);
  run.add_pol_config("e-p-", "ePol-", "pPol-", "-", "-", 0.05);
  run.add_pol_config("e+p+", "ePol+", "pPol+", "+", "+", 0.05);
  
  setup.set_run(run); // Add run to setup
  
  
  // Take cTGCs into account in fit
  // spdlog::info("Creating TGC info.");
  // PrEWUtils::SetupHelp::TGCInfo cTGC_info (
  //   {"SingleW_eminus", "SingleW_eplus", "WW_muplus", "WW_muminus"},
  //   "linear"
  // );
  // setup.add(cTGC_info);
  
  // Set asymmetries and total chiral cross sections scalings as free parameters
  spdlog::info("Creating chiral cross section infos.");

  PrEWUtils::SetupHelp::CrossSectionInfo xs_SingleW_eplus(
    "SingleW_eplus",
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL, PrEW::GlobalVar::Chiral::eLpL}
  );
  PrEWUtils::SetupHelp::CrossSectionInfo xs_SingleW_eminus(
    "SingleW_eminus",
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL, PrEW::GlobalVar::Chiral::eRpR}
  );
  PrEWUtils::SetupHelp::CrossSectionInfo xs_WW_muminus(
    "WW_muminus",
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL}
  );
  PrEWUtils::SetupHelp::CrossSectionInfo xs_WW_muplus(
    "WW_muplus",
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL}
  );
  
  PrEWUtils::SetupHelp::CrossSectionInfo xs_2f_uds_81to101(
    "2f_uds_81to101", 
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL}
  );
  PrEWUtils::SetupHelp::CrossSectionInfo xs_2f_c_81to101(
    "2f_c_81to101", 
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL}
  );
  PrEWUtils::SetupHelp::CrossSectionInfo xs_2f_b_81to101(
    "2f_b_81to101", 
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL}
  );
  PrEWUtils::SetupHelp::CrossSectionInfo xs_2f_mu_81to101(
    "2f_mu_81to101", 
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL}
  );
  
  PrEWUtils::SetupHelp::CrossSectionInfo xs_2f_uds_180to275(
    "2f_uds_180to275", 
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL}
  );
  PrEWUtils::SetupHelp::CrossSectionInfo xs_2f_c_180to275(
    "2f_c_180to275", 
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL}
  );
  PrEWUtils::SetupHelp::CrossSectionInfo xs_2f_b_180to275(
    "2f_b_180to275", 
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL}
  );
  PrEWUtils::SetupHelp::CrossSectionInfo xs_2f_mu_180to275(
    "2f_mu_180to275", 
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL}
  );
  
  xs_SingleW_eplus.use_chiral_asymmetries();
  xs_SingleW_eminus.use_chiral_asymmetries();
  xs_WW_muminus.use_chiral_asymmetries();
  xs_WW_muplus.use_chiral_asymmetries();
  xs_2f_uds_81to101.use_chiral_asymmetries({"Ae_Zpole"});
  xs_2f_c_81to101.use_chiral_asymmetries({"Ae_Zpole"});
  xs_2f_b_81to101.use_chiral_asymmetries({"Ae_Zpole"});
  xs_2f_mu_81to101.use_chiral_asymmetries({"Ae_Zpole"});
  xs_2f_uds_180to275.use_chiral_asymmetries({"Ae_highQ2"});
  xs_2f_c_180to275.use_chiral_asymmetries({"Ae_highQ2"});
  xs_2f_b_180to275.use_chiral_asymmetries({"Ae_highQ2"});
  xs_2f_mu_180to275.use_chiral_asymmetries({"Ae_highQ2"});
  
  xs_SingleW_eplus.use_total_chiral_cross_section();
  xs_SingleW_eminus.use_total_chiral_cross_section();
  xs_WW_muminus.use_total_chiral_cross_section();
  xs_WW_muplus.use_total_chiral_cross_section();
  xs_2f_uds_81to101.use_total_chiral_cross_section();
  xs_2f_c_81to101.use_total_chiral_cross_section();
  xs_2f_b_81to101.use_total_chiral_cross_section();
  xs_2f_mu_81to101.use_total_chiral_cross_section();
  xs_2f_uds_180to275.use_total_chiral_cross_section();
  xs_2f_c_180to275.use_total_chiral_cross_section();
  xs_2f_b_180to275.use_total_chiral_cross_section();
  xs_2f_mu_180to275.use_total_chiral_cross_section();
  
  setup.add(xs_SingleW_eplus);
  setup.add(xs_SingleW_eminus);
  setup.add(xs_WW_muminus);
  setup.add(xs_WW_muplus);
  setup.add(xs_2f_uds_81to101);
  setup.add(xs_2f_c_81to101);
  setup.add(xs_2f_b_81to101);
  setup.add(xs_2f_mu_81to101);
  setup.add(xs_2f_uds_180to275);
  setup.add(xs_2f_c_180to275);
  setup.add(xs_2f_b_180to275);
  setup.add(xs_2f_mu_180to275);
  
  
  // spdlog::info("Creating distribution efficiency infos.");
  // PrEWUtils::SetupHelp::ConstEffInfo 2f_b_eff ("2f_b_81to101", 0.3);
  // 2f_b_eff.constrain(0.3, 0.001);
  // setup.add(2f_b_eff);

  
  spdlog::info("Creating acceptance box (w/ polynomial) infos.");
  PrEWUtils::SetupHelp::AccBoxPolynomialInfo muon_box("MuonAcc");
  muon_box.add_distr("WW_muminus");
  muon_box.add_distr("WW_muplus");
  muon_box.add_distr("xs_2f_mu_81to101");
  muon_box.add_distr("xs_2f_mu_180to275");
  setup.add(muon_box);
  
  //7___________________________________________________________________________
  spdlog::info("Setting up changes of the fit wrt. the toys.");
  PrEWUtils::Setups::FitModifier fit_modifier(energy);
  
  // spdlog::info("Use Af parametrisation in fit.");
  // PrEWUtils::SetupHelp::AfInfo Af_2f_b_81to101("2f_b_81to101");
  // PrEWUtils::SetupHelp::AfInfo Af_2f_mu_81to101("2f_mu_81to101");
  // fit_modifier.add(Af_2f_b_81to101);
  // fit_modifier.add(Af_2f_mu_81to101);
  
  
  //7___________________________________________________________________________
  spdlog::info("Finalizing linking info.");
  setup.complete_setup(); // This must come last in linking setup

  spdlog::info("Create runner (incl. setting up toy generator).");
  PrEWUtils::Runners::ParallelRunner runner ( 
    setup, 
    minuit_minimizers, 
    prew_minimizer 
  );
  
  spdlog::info("Performing instructed fit modifications.");
  runner.modify_fit(fit_modifier);
  
  spdlog::info("Run toys.");
  auto results = runner.run_toy_fits( energy, n_toys, n_threads );
  
  spdlog::info("All threads done, printing first result.");
  spdlog::info(results.at(0));
  
  spdlog::info("Write results to: {}", output_path);
  PrEW::Output::Printer printer (output_path);
  printer.new_setup( energy, runner.get_data_connector() );
  printer.add_fits( results );
  printer.write();
  
  spdlog::info("Complete test done!");
}