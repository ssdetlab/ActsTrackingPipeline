//
// Created by Alon Levi on 09/04/2024.
//

// 3D LOOKUP - takes too long to load, is useful since it resolves 2D histograms which are not 1-to-1
/*
std::unordered_map <Scalar, std::vector<Scalar>> read3DLookup(std::string file) {
    std::ifstream lookupFile(&file[0]);
    std::unordered_map <Scalar, std::vector<Scalar>> lookupTable;
    Scalar x, y, z;
    while (lookupFile >> x >> y >> z) {
        std::vector<Scalar> v{y,z};
        lookupTable[x] = v;
    }
    lookupFile.close();
    return lookupTable;
}
*/