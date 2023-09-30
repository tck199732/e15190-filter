#include <map>
#include <vector>
#include <regex>
#include <string>
#include <stdlib.h>
#include <getopt.h>
#include <filesystem>
namespace fs = std::filesystem;

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include "microball.hpp"
#include "hira.hpp"
#include "nwb.hpp"
#include <utilities/particle.hpp>
#include <utilities/physics.hpp>

struct E15190
{
    const static int MAX_MULTI = 128;
    // impact parameter
    double b;

    // microball
    int uball_multi;
    std::array<int, MAX_MULTI> uball_N;
    std::array<int, MAX_MULTI> uball_Z;
    std::array<double, MAX_MULTI> uball_px;
    std::array<double, MAX_MULTI> uball_py;
    std::array<double, MAX_MULTI> uball_pz;

    // hira
    int hira_multi;
    std::array<int, MAX_MULTI> hira_N;
    std::array<int, MAX_MULTI> hira_Z;
    std::array<double, MAX_MULTI> hira_px;
    std::array<double, MAX_MULTI> hira_py;
    std::array<double, MAX_MULTI> hira_pz;

    // veto wall

    // neutron wall
    int nwb_multi;
    std::array<int, MAX_MULTI> nwb_N;
    std::array<int, MAX_MULTI> nwb_Z;
    std::array<double, MAX_MULTI> nwb_px;
    std::array<double, MAX_MULTI> nwb_py;
    std::array<double, MAX_MULTI> nwb_pz;
};

E15190 e15190;

void Initialize_TTree(TTree *&tree)
{
    // impact parameter
    tree->Branch("b", &e15190.b, "b/D");

    // microball
    tree->Branch("uball_multi", &e15190.uball_multi, "uball_multi/I");
    tree->Branch("uball_N", &e15190.uball_N[0], "uball_N[uball_multi]/I");
    tree->Branch("uball_Z", &e15190.uball_Z[0], "uball_Z[uball_multi]/I");
    tree->Branch("uball_px", &e15190.uball_px[0], "uball_px[uball_multi]/D");
    tree->Branch("uball_py", &e15190.uball_py[0], "uball_py[uball_multi]/D");
    tree->Branch("uball_pz", &e15190.uball_pz[0], "uball_pz[uball_multi]/D");

    // hira
    tree->Branch("hira_multi", &e15190.hira_multi, "hira_multi/I");
    tree->Branch("hira_N", &e15190.hira_N[0], "hira_N[hira_multi]/I");
    tree->Branch("hira_Z", &e15190.hira_Z[0], "hira_Z[hira_multi]/I");
    tree->Branch("hira_px", &e15190.hira_px[0], "hira_px[hira_multi]/D");
    tree->Branch("hira_py", &e15190.hira_py[0], "hira_py[hira_multi]/D");
    tree->Branch("hira_pz", &e15190.hira_pz[0], "hira_pz[hira_multi]/D");

    // nwb
    tree->Branch("nwb_multi", &e15190.nwb_multi, "nwb_multi/I");
    tree->Branch("nwb_N", &e15190.nwb_N[0], "nwb_N[nwb_multi]/I");
    tree->Branch("nwb_Z", &e15190.nwb_Z[0], "nwb_Z[nwb_multi]/I");
    tree->Branch("nwb_px", &e15190.nwb_px[0], "nwb_px[nwb_multi]/D");
    tree->Branch("nwb_py", &e15190.nwb_py[0], "nwb_py[nwb_multi]/D");
    tree->Branch("nwb_pz", &e15190.nwb_pz[0], "nwb_pz[nwb_multi]/D");
}

void Initialize_MicroBall(microball *&uball, const std::string &reaction)
{
    fs::path project_dir = std::getenv("PROJECT_DIR");
    fs::path database_dir = project_dir / "database/e15190/microball/acceptance";
    fs::path path_config = database_dir / "config.dat";
    fs::path path_geometry = database_dir / "geometry.dat";
    fs::path path_threshold = database_dir / "fitted_threshold.dat";

    uball->ConfigurateSetup(reaction, path_config.string());
    uball->ReadGeometryMap(path_geometry.string());
    uball->ReadThresholdKinergyMap(path_threshold.string());
}

bool ReadMicroballParticle(microball *&mb, const particle &ptcl)
{
    bool pass_charge = mb->IsChargedParticle(ptcl.z);
    bool pass_coverage = mb->IsCovered(ptcl.theta_lab_deg, ptcl.phi_deg);
    bool pass_threshold = mb->IsAccepted(ptcl.kinergy_lab, ptcl.theta_lab_deg, ptcl.a, ptcl.z);
    return pass_charge && pass_coverage && pass_threshold;
}

bool ReadHiRAParticle(hira *&hira, const particle &ptcl)
{
    return hira->PassAngularCut(ptcl.theta_lab_deg, ptcl.phi_deg) &&
           hira->PassCharged(ptcl.z) &&
           hira->PassKinergyCut(ptcl.a, ptcl.z, ptcl.kinergy_lab);
}

bool ReadNWBParticle(nwb *&nwb, const particle &ptcl)
{
    return nwb->PassAngularCut(ptcl.theta_lab_deg, ptcl.phi_deg) &&
           nwb->PassKinergyCut(ptcl.a, ptcl.z, ptcl.kinergy_lab);
}

void correct_phi_value(particle &ptcl, microball *&microball)
{
    int ring = microball->GetRingID(ptcl.theta_lab_deg);
    if (ring == -1)
    {
        return;
    }
    double phi_min_in_ring = microball->GetPhiMinInRing(ring);
    double phi_max_in_ring = microball->GetPhiMaxInRing(ring);

    if (ptcl.phi_deg < phi_min_in_ring)
    {
        ptcl.phi_rad += 2. * TMath::Pi();
        ptcl.phi_deg += 360.;
    }
    if (ptcl.phi_deg > phi_max_in_ring)
    {
        ptcl.phi_rad -= 2. * TMath::Pi();
        ptcl.phi_deg -= 360.;
    }
}

void read_reaction(const std::string &reaction, int& beamA, int& beamZ, int& targetA, int& targetZ, int& beam_energy)
{
    {
        std::regex pattern("[A-Z][a-z]");
        std::vector<std::string> tokens;
        std::sregex_iterator iter(reaction.begin(), reaction.end(), pattern);
        std::sregex_iterator end;

        while (iter != end)
        {
            std::smatch match = *iter;
            tokens.push_back(match.str());
            ++iter;
        }

        std::string beam = tokens[0];
        std::string target = tokens[1];

        auto GetZ = [](const std::string &nuclei) -> double
        {
            if (nuclei == "Ca" || nuclei == "ca")
                return 20;
            else if (nuclei == "Ni" || nuclei == "ni")
                return 28;
            else if (nuclei == "Sn" || nuclei == "sn")
                return 50;
            else
                return 0;
        };
        beamZ = GetZ(tokens[0]);
        targetZ = GetZ(tokens[1]);
    }

    {
        std::regex pattern("[0-9]+");
        std::vector<int> tokens;

        std::sregex_iterator iter(reaction.begin(), reaction.end(), pattern);
        std::sregex_iterator end;

        while (iter != end)
        {
            std::smatch match = *iter;
            tokens.push_back(std::stoi(match.str()));
            ++iter;
        }
        beamA = tokens[0];
        targetA = tokens[1];
        beam_energy = tokens[2];
    }
}
