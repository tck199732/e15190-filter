#ifndef particle_hh
#define particle_hh

#include <string>
#include "physics.hpp"

struct particle
{
    int n, z;
    double px_a, py_a, pz_a, mass;
    std::string frame;

    int a;
    // same in lab and cms
    double px, py, phi_rad, phi_deg, pmag_trans;

    // cms quantities
    double pz_cms;
    double theta_cms_rad, theta_cms_deg, kinergy_cms, pmag_cms, rapidity_cms;

    // lab quantities
    double pz_lab;
    double theta_lab_rad, theta_lab_deg, kinergy_lab, pmag_lab, rapidity_lab;

    void initialize(const double &beta)
    {
        // initialize frame-independent quantities
        this->a = this->n + this->z;
        this->px = this->px_a * this->a;
        this->py = this->py_a * this->a;
        this->phi_rad = get_phi_rad(this->px, this->py);
        this->phi_deg = this->phi_rad * TMath::RadToDeg();
        this->pmag_trans = get_pt(this->px, this->py);
        double gamma = 1. / std::sqrt(1 - pow(beta, 2.));

        if (this->frame == "cms")
        {
            // initialize cms quantities
            this->pz_cms = this->pz_a * this->a;
            this->pmag_cms = get_pmag(this->pmag_trans, this->pz_cms);
            this->kinergy_cms = get_kinergy(this->mass, this->pmag_cms);
            this->theta_cms_rad = get_theta_rad(this->pmag_trans, this->pz_cms);
            this->theta_cms_deg = this->theta_cms_rad * TMath::RadToDeg();
            this->rapidity_cms = get_rapidity(this->kinergy_cms, this->pz_cms, this->mass);
            // construct lab quantities
            this->pz_lab = boostz(this->mass, this->pz_cms, this->kinergy_cms, -beta);
            this->pmag_lab = get_pmag(this->pmag_trans, this->pz_lab);
            this->kinergy_lab = get_kinergy(this->mass, this->pmag_lab);
            this->theta_lab_rad = get_theta_rad(this->pmag_trans, this->pz_lab);
            this->theta_lab_deg = this->theta_lab_rad * TMath::RadToDeg();
            this->rapidity_lab = get_rapidity(this->kinergy_lab, this->pz_lab, this->mass);
        }

        else if (this->frame == "lab")
        {
            // initialize lab quantities
            this->pz_lab = this->pz_a * this->a;
            this->pmag_lab = get_pmag(this->pmag_trans, this->pz_lab);
            this->kinergy_lab = get_kinergy(this->mass, this->pmag_lab);
            this->theta_lab_rad = get_theta_rad(this->pmag_trans, this->pz_lab);
            this->theta_lab_deg = this->theta_lab_rad * TMath::RadToDeg();
            this->rapidity_lab = get_rapidity(this->kinergy_lab, this->pz_lab, this->mass);
            // construct cms quantities
            this->pz_cms = boostz(this->mass, this->pz_lab, this->kinergy_lab, beta);
            this->pmag_cms = get_pmag(this->pmag_trans, this->pz_cms);
            this->kinergy_cms = get_kinergy(this->mass, this->pmag_cms);
            this->theta_cms_rad = get_theta_rad(this->pmag_trans, this->pz_cms);
            this->theta_cms_deg = this->theta_cms_rad * TMath::RadToDeg();
            this->rapidity_cms = get_rapidity(this->kinergy_cms, this->pz_cms, this->mass);
        }
    }
};

#endif