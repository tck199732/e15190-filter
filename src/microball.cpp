#include "microball.hpp"

microball::microball()
{
    this->Is_apply_cut_charged_particle = 1;
    this->Is_apply_cut_multiple_hit = 1;
    this->Is_apply_cut_kinergy = 1;
    this->Is_apply_cut_coverage = 1;
}

void microball::ReadGeometryMap(const std::string &filename)
{
    if (!fs::exists(filename))
    {
        std::string msg = Form("file does not exist : %s", filename.c_str());
        throw std::invalid_argument(msg.c_str());
    }
    std::ifstream infile(filename.c_str());
    infile.ignore(99, '\n');

    int ring, det;
    double theta1, theta2, phi1, phi2;
    while (infile >> ring)
    {
        infile >> det >> theta1 >> theta2 >> phi1 >> phi2;
        if (this->DetectorSetupMap.count(ring) == 1)
        {
            this->ThetaMap[ring] = {theta1, theta2};
            if (std::find(this->DetectorSetupMap[ring].begin(), this->DetectorSetupMap[ring].end(), det) == this->DetectorSetupMap[ring].end())
            {
                continue;
            }
            this->PhiMap[{ring, det}] = {phi1, phi2};
        }
    }
    infile.close();
}

void microball::ConfigurateSetup(const std::string &reaction, const std::string &filename)
{
    if (!fs::exists(filename))
    {
        std::string msg = Form("file does not exist : %s", filename.c_str());
        throw std::invalid_argument(msg.c_str());
    }
    std::ifstream stream(filename.c_str());
    stream.ignore(99, '\n');
    std::string line;
    while (std::getline(stream, line))
    {
        std::istringstream iss(line);
        std::string sys;
        int ring_id, det_id;
        iss >> sys >> ring_id;
        if (sys != reaction)
        {
            continue;
        }

        this->DetectorSetupMap[ring_id] = std::vector<int>();
        while (iss >> det_id)
        {
            this->DetectorSetupMap[ring_id].push_back(det_id);
        }
    }
    return;
}

void microball::ReadThresholdKinergyMap(const std::string &filename)
{
    if (!fs::exists(filename))
    {
        std::string msg = Form("file does not exist : %s", filename.c_str());
        throw std::invalid_argument(msg.c_str());
    }
    std::ifstream infile(filename.c_str());
    infile.ignore(99, '\n');

    int ring, A, Z;
    double kinergy_MeV;
    while (infile >> ring)
    {
        if (this->DetectorSetupMap.count(ring) == 0)
        {
            continue;
        }
        infile >> A >> Z >> kinergy_MeV;
        if (A > this->MaxA || Z > this->MaxZ)
        {
            continue;
        }
        this->KinergyThresholdMap[{ring, A, Z}] = kinergy_MeV;
    }
}

std::pair<int, int> microball::GetRingDetID(const double &thetalab, const double &phi)
{
    int ring_id = this->GetRingID(thetalab);
    if (ring_id == -1)
    {
        return {-1, -1};
    }

    for (const auto &[ring_det, phi_range] : this->PhiMap)
    {
        if (ring_det[0] != ring_id)
        {
            continue;
        }
        if (phi >= phi_range[0] && phi < phi_range[1])
        {
            return {ring_det[0], ring_det[1]};
        }
    }
    return {-1, -1};
}

int microball::GetRingID(const double &thetalab)
{
    for (auto &[ring, theta_range] : this->ThetaMap)
    {
        if (thetalab >= theta_range[0] && thetalab < theta_range[1])
        {
            return ring;
        }
    }
    return -1;
}

int microball::GetDetID(const double &thetalab, const double &phi)
{
    return GetRingDetID(thetalab, phi).second;
}

double microball::GetPhiMinInRing(const int &ring)
{
    double phi_min = DBL_MAX;
    for (const auto &[ring_det, phi_range] : this->PhiMap)
    {
        if (ring_det[0] == ring && phi_range[0] < phi_min)
        {
            phi_min = phi_range[0];
        }
    }
    return phi_min;
}

double microball::GetPhiMaxInRing(const int &ring)
{
    double phi_max = DBL_MIN;
    for (const auto &[ring_det, phi_range] : this->PhiMap)
    {
        if (ring_det[0] == ring && phi_range[1] > phi_max)
        {
            phi_max = phi_range[1];
        }
    }
    return phi_max;
}

double microball::GetThresholdKinergy(const int &ring_id, const int &aid, const int &zid)
{
    if (this->DetectorSetupMap.count(ring_id) == 0)
    {
        return DBL_MAX;
    }
    return this->KinergyThresholdMap[{ring_id, aid, zid}];
}
double microball::GetThresholdKinergy(const double &thetalab, const int &aid, const int &zid)
{
    int ring_id = this->GetRingID(thetalab);
    return this->GetThresholdKinergy(ring_id, aid, zid);
}

bool microball::IsChargedParticle(const int &Z)
{
    if (!this->Is_apply_cut_charged_particle)
    {
        return true;
    }
    return Z > 0;
}

bool microball::IsCovered(const double &thetalab, const double &phi)
{
    if (!this->Is_apply_cut_coverage)
    {
        return true;
    }
    int ring_id = this->GetRingID(thetalab);
    int det_id = this->GetDetID(thetalab, phi);
    return (ring_id != -1 && det_id != -1);
}

bool microball::IsAccepted(const double &ekinlab, const double &thetalab, const int &aid, const int &zid)
{
    if (!this->Is_apply_cut_kinergy)
    {
        return true;
    }
    int ring_id = this->GetRingID(thetalab);
    return (ring_id != -1 && aid <= this->MaxA && zid <= this->MaxZ) ? ekinlab >= this->KinergyThresholdMap[{ring_id, aid, zid}] : false;
}

void microball::ResetCsIHitMap()
{
    for (const auto &[ring_det, _] : this->CsIHitMap)
    {
        this->CsIHitMap[ring_det] = 0;
    }
    return;
}

int microball::GetCsIHits()
{
    int Mch = 0;
    for (auto &[ring_det, hits] : this->CsIHitMap)
    {
        int ring = ring_det[0];
        int det = ring_det[1];

        if (!Is_apply_cut_multiple_hit)
        {
            if (ring == -1 && det == -1)
            {
                if (!Is_apply_cut_coverage)
                    Mch += hits;
            }
            else
            {
                Mch += hits;
            }
        }
        else
        {
            if (hits > 0)
            {
                if (ring == -1 && det == -1)
                {
                    if (!Is_apply_cut_coverage)
                        Mch += 1;
                }
                else
                {
                    Mch += 1;
                }
            }
        }
    }
    return Mch;
}

int microball::GetCsIHits(const int &ring, const int &det)
{
    if (this->CsIHitMap.count({ring, det}) == 0)
    {
        return 0;
    }
    return this->CsIHitMap[{ring, det}];
}

void microball::AddCsIHit(const double &thetalab, const double &phi)
{
    int ring_id = this->GetRingID(thetalab);
    int det_id = this->GetDetID(thetalab, phi);

    if (this->CsIHitMap.count({ring_id, det_id}) == 0)
    {
        this->CsIHitMap[{ring_id, det_id}] = 1;
    }
    else if (this->CsIHitMap.count({ring_id, det_id}) == 1)
    {
        this->CsIHitMap[{ring_id, det_id}]++;
    }
    return;
}

void microball::ViewDetectorSetupMap()
{
    for (auto &[ring, det_array] : this->DetectorSetupMap)
    {
        std::cout << ring << " ";
        for (auto &det : det_array)
        {
            std::cout << det << " ";
        }
        std::cout << std::endl;
    }
}

void microball::ViewGeometryMap()
{
    for (auto &[ring_det, _] : this->PhiMap)
    {
        int ring = ring_det[0];
        int det = ring_det[1];

        auto print = [](auto... args)
        {
            ((std::cout << args << " "), ...) << std::endl;
        };

        print(ring, det, this->ThetaMap[ring][0], this->ThetaMap[ring][1], this->PhiMap[{ring, det}][0], this->PhiMap[{ring, det}][1]);
    }
}

void microball::ViewThresholdKinergyMap()
{
    for (auto &[ring_A_Z, kinergy] : this->KinergyThresholdMap)
    {
        int ring = ring_A_Z[0];
        int A = ring_A_Z[1];
        int Z = ring_A_Z[2];

        auto print = [](auto... args)
        {
            ((std::cout << args << " "), ...) << std::endl;
        };

        print(ring, A, Z, kinergy);
    }
}