#ifndef nwb_hh
#define nwb_hh

#include <iostream>
#include <fstream>
#include <array>
#include <map>
#include <filesystem>
namespace fs = std::filesystem;

#include "TMath.h"
#include "TString.h"

class nwb
{
protected:

    // kinergy cut for experiment e15190
    std::map<std::string, std::array<double, 2>> KINERGYCUT = {
        {"n", {10.0, 308.0}},
    };

    std::map<std::pair<int, int>, std::string> MASS_CHARGE_TO_NAME = {
        {{1, 0}, "n"},
    };

    std::array<double, 2> THETACUT = {28., 51.5};
    std::array<double, 2> PHICUT = {-360, 360};

public:
    nwb();
    ~nwb(){};
    bool PassKinergyCut(const std::string &particle_name, const double &kinergy_per_nucleon);
    bool PassKinergyCut(const int &A, const int &Z, const double &kinergy);
    bool PassAngularCut(const double &theta_deg, const double &phi_deg);
    void ResetCounter() { mCounterPass = mCounterFail = 0; }
    long GetCountPass() { return mCounterPass; }
    long GetCountFail() { return mCounterFail; }
    void CountPass() { mCounterPass++; }
    void CountFail() { mCounterFail++; }

private:
    long mCounterPass, mCounterFail;
    std::map<std::string, std::array<double, 2>> mKinergyCut;
    std::array<double, 2> mThetaCut;
    std::array<double, 2> mPhiCut;
    std::map<std::pair<int, int>, std::string> mMassChargeToName;
};

#endif
