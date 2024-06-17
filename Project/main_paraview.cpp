//Librerie
#include "src_paraview/paraview.hpp"
#include <iostream>
#include "src/FracturesNetworkLibrary.hpp"
#include "src/Utils.hpp"
#include "src/Geometry.hpp"

using namespace std;
using namespace FractureNetwork;

int main()
{
    string filepath = "DFN_files/FR10_data.txt";

    vector<Fracture> fratture;
    vector<Trace> tracce;

    bool done = ImportFractures(filepath, fratture);
    if(!done){
        cout << "Errore nell'apertura del file" << endl;
        return 1;
    }

    CalculateTraces(fratture,tracce);
    SortingFractureTraces(fratture,tracce);

    Mesh mesh;
    list<MatrixXd> cutted = cutting(fratture,tracce);
    vector<MatrixXd> c;
    c.reserve(cutted.size());
    move(begin(cutted), end(cutted), back_inserter(c));

    Export_paraview(fratture,tracce);
    Export_paraview2(c);

    return 0;
}

