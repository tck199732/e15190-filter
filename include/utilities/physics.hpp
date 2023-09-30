#ifndef physics_hh
#define physics_hh

#include "TMath.h"
double get_reaction_beta(const double &beam_mass, const double &targ_mass, const double &beam_ea, const int &beam_a) 
{
    double beam_ke = beam_ea * beam_a;
    double beam_energy_tot = beam_ke + beam_mass;
    double mom_beam = std::sqrt(pow(beam_ke, 2.) + 2. * beam_ke * beam_mass);
    double gamma = beam_energy_tot / beam_mass;
    return mom_beam / (gamma * beam_mass + targ_mass);
}

double get_beam_rapidity(const double &beam_mass, const double &beam_ea, const int &beam_a)
{
    double beam_ke = beam_ea * beam_a;
    double beam_energy_tot = beam_ke + beam_mass;
    double mom_beam = std::sqrt(pow(beam_ke, 2.) + 2. * beam_ke * beam_mass);
    return 0.5 * TMath::Log((beam_energy_tot + mom_beam) / (beam_energy_tot - mom_beam));
}

double get_pt(const double &px, const double &py) {return std::sqrt(pow(px, 2.) + pow(py, 2.));}

double get_pmag(const double &pt, const double &pz) { return std::sqrt(pow(pt, 2.) + pow(pz, 2.));}

double get_pmag(const double &px, const double &py, const double &pz) { return get_pmag(get_pt(px, py), pz); }

double get_kinergy(const double &mass, const double &p) { return std::sqrt(pow(p, 2.) + pow(mass, 2.)) - mass; }

double get_kinergy(const double &mass, const double &px, const double &py, const double &pz) { return get_kinergy(mass, get_pmag(px, py, pz));}

double get_phi_rad(const double &px, const double &py) { return TMath::ATan2(py, px); }
double get_theta_rad(const double &pt, const double &pz) { return TMath::ATan2(pt, pz); }
double get_phi_deg(const double &px, const double &py) { return get_phi_rad(px, py) * TMath::RadToDeg(); }
double get_theta_deg(const double &pt, const double &pz) { return get_theta_rad(pt, pz) * TMath::RadToDeg(); }

double boostz(const double &mass, const double &pz, const double &ekin, const double &betacms) { return 1. / std::sqrt(1 - pow(betacms, 2.)) * (pz - betacms * (ekin + mass));}

double get_rapidity(const double &ekin, const double &pz, const double &mass){ return 0.5 * TMath::Log((ekin + pz + mass) / (ekin - pz + mass));}

#endif