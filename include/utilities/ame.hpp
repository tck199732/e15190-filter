#ifndef ame_hh
#define ame_hh
#include <map>
#include <fstream>
#include <filesystem>
namespace fs = std::filesystem;

std::map<std::pair<int, int>, double> get_ame_mass_table()
{
    std::map<std::pair<int, int>, double> ame_mass_table;
    fs::path PROJECT_DIR = std::getenv("PROJECT_DIR");
    fs::path path = PROJECT_DIR / "database/ame/ame20.csv";
    std::ifstream infile(path.c_str());
    infile.ignore(99, '\n');

    std::string line;
    while (std::getline(infile, line)) {
        std::vector<std::string> row;
        std::stringstream lineStream(line);
        std::string cell;
        // Split each line into individual cells using a comma as the delimiter
        while (std::getline(lineStream, cell, ',')) {
            row.push_back(cell);
        }
        int Z = std::stoi(row[0]);
        int A = std::stoi(row[1]);
        double mass = std::stod(row[2]);
        std::pair pair = std::make_pair(Z, A);
        ame_mass_table[pair] = mass;
    }
    return ame_mass_table;
}

#endif