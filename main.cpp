#include "main.hpp"
#include "tree.hpp"
#include <argparse/argparse.hpp>
#include <utilities/ame.hpp>
#include <progressbar/progressbar.hpp>

argparse::ArgumentParser get_arguments(int argc, char *argv[]);
void filter(TChain *&chain, TTree *&tree, argparse::ArgumentParser &program);
void debug(argparse::ArgumentParser &program);
void save_file(const std::string &output_file, TTree *tree);

int main(int argc, char *argv[])
{
    argparse::ArgumentParser program = get_arguments(argc, argv);
    debug(program);
    TChain *chain = new TChain(program.get<std::string>("--chain-name").c_str());
    Initialize_TChain(chain, program.get<std::vector<std::string>>("--input-files"));
    TTree *tree = new TTree(program.get<std::string>("--chain-name").c_str(), "");
    Initialize_TTree(tree);

    filter(chain, tree, program);
    save_file(program.get<std::string>("--output-file"), tree);
    return 0;
}

argparse::ArgumentParser get_arguments(int argc, char *argv[])
{
    argparse::ArgumentParser program("e15190-filter", "0.1.0");
    program.add_argument("-r", "--reaction", "reaction name")
        .default_value(std::string{"Ca40Ni58E56"})
        .action([](const std::string &value)
                {
            static const std::vector<std::string> choices = { 
                "Ca40Ni58E56", "Ca40Ni58E140", "Ca48Ni64E56",  "Ca48Ni64E140",
                "Ca40Sn112E56", "Ca40Sn112E140", "Ca48Sn124E56", "Ca48Sn124E140"
            };
            if (std::find(choices.begin(), choices.end(), value) != choices.end()) {
                return value;
            }
            return std::string{ "Ca40Ni58E56" }; })
        .required();

    program.add_argument("-i", "--input-files", "list of input paths, separated by comma")
        .nargs(argparse::nargs_pattern::at_least_one)
        .required();

    program.add_argument("-o", "--output-file", "output file path")
        .default_value(std::string("output.root"))
        .required();

    program.add_argument("-n", "--chain-name", "name of the TChain")
        .default_value(std::string("AMD"))
        .required();

    program.add_argument("-f", "--reference-frame", "reference frame of the input particles")
        .default_value(std::string("cms"))
        .action([](const std::string &value)
                {
            static const std::vector<std::string> choices = { "cms", "lab" };
            if (std::find(choices.begin(), choices.end(), value) != choices.end()) {
                return value;
            }
            return std::string{ "cms" }; })
        .required();

    program.add_argument("-d", "--debug", "debug mode")
        .default_value(false)
        .implicit_value(true);

    program.add_argument("--progress-bar", "show progress bar")
        .default_value(false)
        .implicit_value(true);

    try
    {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error &err)
    {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    return program;
}

void filter(TChain *&chain, TTree *&tree, argparse::ArgumentParser &program)
{
    auto ame_mass_table = get_ame_mass_table();
    int beamA, beamZ, targetA, targetZ, beam_energy;
    read_reaction(program.get<std::string>("--reaction"), beamA, beamZ, targetA, targetZ, beam_energy);
    double beam_mass = ame_mass_table[{beamZ, beamA}];
    double target_mass = ame_mass_table[{targetZ, targetA}];
    double betacms = get_reaction_beta(beam_mass, target_mass, beam_energy, beamA);

    microball *uball = new microball();
    Initialize_MicroBall(uball, program.get<std::string>("--reaction"));
    hira *hira10 = new hira();
    nwb *nwallb = new nwb();

    ProgressBar pbar(chain->GetEntries(), "Progress");
    for (int ievt = 0; ievt < chain->GetEntries(); ievt++)
    {
        chain->GetEntry(ievt);
        uball->ResetCsIHitMap();
        hira10->ResetCounter();
        nwallb->ResetCounter();

        for (unsigned int i = 0; i < amd.multi; i++)
        {
            double mass = ame_mass_table[{amd.Z[i], amd.N[i] + amd.Z[i]}];
            particle ptcl = {amd.N[i], amd.Z[i], amd.px[i], amd.py[i], amd.pz[i], mass, "cms"};
            ptcl.initialize(betacms);

            // phi is calculated according to microball detector, if the particle is not covered by microball, phi is not correct and should be in the range of [-pi, pi].
            correct_phi_value(ptcl, uball);

            int uball_multi = uball->GetCsIHits();
            int hira_multi = hira10->GetCountPass();
            int nwb_multi = nwallb->GetCountPass();

            if (ReadMicroballParticle(uball, ptcl))
            {
                e15190.uball_N[uball_multi] = ptcl.n;
                e15190.uball_Z[uball_multi] = ptcl.z;
                e15190.uball_px[uball_multi] = ptcl.px;
                e15190.uball_py[uball_multi] = ptcl.py;
                e15190.uball_pz[uball_multi] = ptcl.pz_lab;
                uball->AddCsIHit(ptcl.theta_lab_deg, ptcl.phi_deg);
            }

            if (ReadHiRAParticle(hira10, ptcl))
            {
                e15190.hira_N[hira_multi] = ptcl.n;
                e15190.hira_Z[hira_multi] = ptcl.z;
                e15190.hira_px[hira_multi] = ptcl.px;
                e15190.hira_py[hira_multi] = ptcl.py;
                e15190.hira_pz[hira_multi] = ptcl.pz_lab;
                hira10->CountPass();
            }

            if (ReadNWBParticle(nwallb, ptcl))
            {
                e15190.nwb_N[nwb_multi] = ptcl.n;
                e15190.nwb_Z[nwb_multi] = ptcl.z;
                e15190.nwb_px[nwb_multi] = ptcl.px;
                e15190.nwb_py[nwb_multi] = ptcl.py;
                e15190.nwb_pz[nwb_multi] = ptcl.pz_lab;
                nwallb->CountPass();
            }
        }

        e15190.hira_multi = hira10->GetCountPass();
        e15190.nwb_multi = nwallb->GetCountPass();
        e15190.uball_multi = uball->GetCsIHits();
        e15190.b = amd.b;
        // if microball multi is 0, in experiment we don't see the event. We still keep the event here as this data can be easily removed in the analysis.
        tree->Fill();

        if (program.get<bool>("--progress-bar") == true)
            pbar.Update();
    }
    if (program.get<bool>("--progress-bar") == true)
        pbar.Finish();
}

void save_file(const std::string &output_file, TTree *tree)
{
    TFile *outputfile = new TFile(output_file.c_str(), "recreate");
    outputfile->cd();
    tree->Write();
    outputfile->Write();
    outputfile->Close();
    return;
}

void debug(argparse::ArgumentParser &program)
{
    if (program.get<bool>("--debug"))
    {
        std::cout << "reaction: " << program.get<std::string>("--reaction") << std::endl;
        std::cout << "input files: ";
        for (const auto &file : program.get<std::vector<std::string>>("--input-files"))
        {
            std::cout << file << " ";
        }
        std::cout << std::endl;
        std::cout << "output file: " << program.get<std::string>("--output-file") << std::endl;
        std::cout << "chain name: " << program.get<std::string>("--chain-name") << std::endl;
        std::cout << "reference frame: " << program.get<std::string>("--reference-frame") << std::endl;
        std::exit(0);
    }
    return;
}
