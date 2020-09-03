// Includes from PrEW
#include "GlobalVar/Chiral.h"
#include "Output/Printer.h"

// Includes from PrEWUtils
#include "Runners/ParallelRunner.h"
#include "SetupHelp/AccBoxInfo.h"
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
  setup.add_input_file("/afs/desy.de/group/flc/pool/beyerjac/TGCAnalysis/PrEW/testdata/RK_examplefile_500_250_2018_04_03.root", "RK");
  setup.add_input_file("/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/WW_charge_separated/distributions/combined/Distribution_250GeV_WW_semilep_MuAntiNu.root", "RK");
  setup.add_input_file("/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/WW_charge_separated/distributions/combined/Distribution_250GeV_WW_semilep_AntiMuNu.root", "RK");
  setup.add_input_file("/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/Z_2f_separated/PrEW_inputs/Z2f_separated_cosTheta_distributions.root", "RK");

  spdlog::info("Selecting distributions.");
  setup.use_distr("singleWplussemileptonic");
  setup.use_distr("singleWminussemileptonic");
  setup.use_distr("WW_semilep_MuAntiNu");
  setup.use_distr("WW_semilep_AntiMuNu");
  setup.use_distr("Zuds_81to101");
  setup.use_distr("Zuds_180to275");
  setup.use_distr("Zcc_81to101");
  setup.use_distr("Zcc_180to275");
  setup.use_distr("Zbb_81to101");
  setup.use_distr("Zbb_180to275");
  setup.use_distr("Zmumu_81to101");
  setup.use_distr("Zmumu_180to275");
  
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
  spdlog::info("Creating TGC info.");
  PrEWUtils::SetupHelp::TGCInfo cTGC_info (
    {"singleWplussemileptonic", "singleWminussemileptonic", "WW_semilep_MuAntiNu", "WW_semilep_AntiMuNu"},
    "linear"
  );
  setup.add(cTGC_info);
  
  // Set asymmetries and total chiral cross sections scalings as free parameters
  spdlog::info("Creating chiral cross section infos.");

  PrEWUtils::SetupHelp::CrossSectionInfo xs_singleWminus(
    "singleWminussemileptonic",
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL, PrEW::GlobalVar::Chiral::eLpL}
  );
  PrEWUtils::SetupHelp::CrossSectionInfo xs_singleWplus(
    "singleWplussemileptonic",
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL, PrEW::GlobalVar::Chiral::eRpR}
  );
  PrEWUtils::SetupHelp::CrossSectionInfo xs_WW_MuAntiNu(
    "WW_semilep_MuAntiNu",
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL}
  );
  PrEWUtils::SetupHelp::CrossSectionInfo xs_WW_AntiMuNu(
    "WW_semilep_AntiMuNu",
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL}
  );
  
  PrEWUtils::SetupHelp::CrossSectionInfo xs_Zuds_81to101(
    "Zuds_81to101", 
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL}
  );
  PrEWUtils::SetupHelp::CrossSectionInfo xs_Zcc_81to101(
    "Zcc_81to101", 
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL}
  );
  PrEWUtils::SetupHelp::CrossSectionInfo xs_Zbb_81to101(
    "Zbb_81to101", 
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL}
  );
  PrEWUtils::SetupHelp::CrossSectionInfo xs_Zmumu_81to101(
    "Zmumu_81to101", 
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL}
  );
  
  PrEWUtils::SetupHelp::CrossSectionInfo xs_Zuds_180to275(
    "Zuds_180to275", 
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL}
  );
  PrEWUtils::SetupHelp::CrossSectionInfo xs_Zcc_180to275(
    "Zcc_180to275", 
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL}
  );
  PrEWUtils::SetupHelp::CrossSectionInfo xs_Zbb_180to275(
    "Zbb_180to275", 
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL}
  );
  PrEWUtils::SetupHelp::CrossSectionInfo xs_Zmumu_180to275(
    "Zmumu_180to275", 
    {PrEW::GlobalVar::Chiral::eLpR, PrEW::GlobalVar::Chiral::eRpL}
  );
  
  xs_singleWminus.use_chiral_asymmetries();
  xs_singleWplus.use_chiral_asymmetries();
  xs_WW_MuAntiNu.use_chiral_asymmetries();
  xs_WW_AntiMuNu.use_chiral_asymmetries();
  xs_Zuds_81to101.use_chiral_asymmetries({"Ae_Zpole"});
  xs_Zcc_81to101.use_chiral_asymmetries({"Ae_Zpole"});
  xs_Zbb_81to101.use_chiral_asymmetries({"Ae_Zpole"});
  xs_Zmumu_81to101.use_chiral_asymmetries({"Ae_Zpole"});
  xs_Zuds_180to275.use_chiral_asymmetries({"Ae_highQ2"});
  xs_Zcc_180to275.use_chiral_asymmetries({"Ae_highQ2"});
  xs_Zbb_180to275.use_chiral_asymmetries({"Ae_highQ2"});
  xs_Zmumu_180to275.use_chiral_asymmetries({"Ae_highQ2"});
  
  xs_singleWminus.use_total_chiral_cross_section();
  xs_singleWplus.use_total_chiral_cross_section();
  xs_WW_MuAntiNu.use_total_chiral_cross_section();
  xs_WW_AntiMuNu.use_total_chiral_cross_section();
  xs_Zuds_81to101.use_total_chiral_cross_section();
  xs_Zcc_81to101.use_total_chiral_cross_section();
  xs_Zbb_81to101.use_total_chiral_cross_section();
  xs_Zmumu_81to101.use_total_chiral_cross_section();
  xs_Zuds_180to275.use_total_chiral_cross_section();
  xs_Zcc_180to275.use_total_chiral_cross_section();
  xs_Zbb_180to275.use_total_chiral_cross_section();
  xs_Zmumu_180to275.use_total_chiral_cross_section();
  
  setup.add(xs_singleWminus);
  setup.add(xs_singleWplus);
  setup.add(xs_WW_MuAntiNu);
  setup.add(xs_WW_AntiMuNu);
  setup.add(xs_Zuds_81to101);
  setup.add(xs_Zcc_81to101);
  setup.add(xs_Zbb_81to101);
  setup.add(xs_Zmumu_81to101);
  setup.add(xs_Zuds_180to275);
  setup.add(xs_Zcc_180to275);
  setup.add(xs_Zbb_180to275);
  setup.add(xs_Zmumu_180to275);
  
  
  spdlog::info("Creating distribution efficiency infos.");
  PrEWUtils::SetupHelp::ConstEffInfo Zbb_eff ("Zbb_81to101", 0.3);
  Zbb_eff.constrain(0.3, 0.001);
  setup.add(Zbb_eff);
  
  
  spdlog::info("Creating acceptance box infos.");
  PrEWUtils::SetupHelp::AccBoxInfo muon_box("MuonAcceptance", "costheta", 0, 1.9);
  muon_box.add_distr("Zmumu_81to101");
  muon_box.add_distr("Zmumu_180to275");
  setup.add(muon_box);
  
  
  //7___________________________________________________________________________
  spdlog::info("Setting up changes of the fit wrt. the toys.");
  PrEWUtils::Setups::FitModifier fit_modifier(energy);
  
  spdlog::info("Use Af parametrisation in fit.");
  PrEWUtils::SetupHelp::AfInfo Af_Zbb_81to101("Zbb_81to101");
  PrEWUtils::SetupHelp::AfInfo Af_Zmumu_81to101("Zmumu_81to101");
  fit_modifier.add(Af_Zbb_81to101);
  fit_modifier.add(Af_Zmumu_81to101);
  
  
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