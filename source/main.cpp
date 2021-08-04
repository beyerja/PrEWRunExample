// Includes from PrEW
#include "GlobalVar/Chiral.h"
#include "Output/Printer.h"

// Includes from PrEWUtils
#include "Runners/ParallelRunner.h"
#include "SetupHelp/SetupInfos.h"
#include "Setups/FitModifier.h"
#include "Setups/GeneralSetup.h"

// Standard library
#include <string>

#include "spdlog/spdlog.h"

int main(int /*argc*/, char ** /*argv*/) {
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
  setup.add_input_files("/nfs/dust/ilc/group/ild/beyerjac/TGCAnalysis/SampleProduction/NewMCProduction/2f_Z_l/PrEWInput/MuAcc_costheta_0.9925", ".*\\.csv", "CSV");

  spdlog::info("Selecting distributions.");
  setup.use_distr("2f_mu_81to101_FZ");
  setup.use_distr("2f_mu_81to101_BZ");
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
  run.add_lumi_constr(2000, 2000 * 3e-3);
  run.add_pol_constr("ePol-", 0.8, 0.8 * 2.5e-3);
  run.add_pol_constr("ePol+", 0.8, 0.8 * 2.5e-3);
  run.add_pol_constr("pPol-", 0.3, 0.3 * 2.5e-3);
  run.add_pol_constr("pPol+", 0.3, 0.3 * 2.5e-3);

  // Set up the desired polarisation sharings
  run.add_pol_config("e-p+", "ePol-", "pPol+", "-", "+", 0.45);
  run.add_pol_config("e+p-", "ePol+", "pPol-", "+", "-", 0.45);
  run.add_pol_config("e-p-", "ePol-", "pPol-", "-", "-", 0.05);
  run.add_pol_config("e+p+", "ePol+", "pPol+", "+", "+", 0.05);

  setup.set_run(run); // Add run to setup

  // spdlog::info("Creating distribution efficiency infos.");
  // PrEWUtils::SetupHelp::ConstEffInfo 2f_mu_eff ("2f_mu_180to275", 0.3);
  // 2f_mu_eff.constrain(0.3, 0.001);
  // setup.add(2f_mu_eff);

  spdlog::info("Creating acceptance box (w/ polynomial) infos.");
  PrEWUtils::SetupHelp::AccBoxPolynomialInfo muon_box("MuonAcc");
  muon_box.add_distr("2f_mu_81to101_FZ");
  muon_box.add_distr("2f_mu_81to101_BZ");
  muon_box.add_distr("2f_mu_180to275");
  setup.add(muon_box);

  // 7___________________________________________________________________________
  spdlog::info("Setting up changes of the fit wrt. the toys.");
  PrEWUtils::Setups::FitModifier fit_modifier(energy);

  spdlog::info("Use generalised 2f parametrisation in fit.");
  auto pars_2f_mu_81to101 = PrEWUtils::SetupHelp::DifermionPars()
                                .s0("s0_2f_mu_81to101")
                                .Ae("Ae_2f_mu_81to101", 0.2)
                                .Af("Af_2f_mu_81to101", 0.2)
                                .ef("ef_2f_mu_81to101", 0.02)
                                .k0("k0_2f_mu_81to101", 0.07)
                                .dk("dk_2f_mu_81to101", 0.0006);
  fit_modifier.add(PrEWUtils::SetupHelp::DifermionParamInfo(
      "2f_mu_81to101_FZ", pars_2f_mu_81to101));
  fit_modifier.add(PrEWUtils::SetupHelp::DifermionParamInfo(
      "2f_mu_81to101_BZ", pars_2f_mu_81to101));

  auto pars_2f_mu_180to275 =
      PrEWUtils::SetupHelp::DifermionPars().Ae(0.11).Af(0.03).ef(1.4);
  fit_modifier.add(PrEWUtils::SetupHelp::DifermionParamInfo(
      "2f_mu_180to275", pars_2f_mu_180to275));

  // 7___________________________________________________________________________
  spdlog::info("Finalizing linking info.");
  setup.complete_setup(); // This must come last in linking setup

  spdlog::info("Create runner (incl. setting up toy generator).");
  PrEWUtils::Runners::ParallelRunner runner(setup, minuit_minimizers,
                                            prew_minimizer);

  spdlog::info("Performing instructed fit modifications.");
  runner.modify_fit(fit_modifier);

  spdlog::info("Run toys.");
  auto results = runner.run_toy_fits(energy, n_toys, n_threads);

  spdlog::info("All threads done, printing first result.");
  spdlog::info(results.at(0));

  spdlog::info("Write results to: {}", output_path);
  PrEW::Output::Printer printer(output_path);
  printer.new_setup(energy, runner.get_data_connector());
  printer.add_fits(results);
  printer.write();

  spdlog::info("Complete test done!");
}