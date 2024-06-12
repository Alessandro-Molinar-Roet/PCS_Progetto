//Librerie
#include "FracturesNetworkLibrary.hpp"
#include "Utils.hpp"
#include "Geometry.hpp"
#include <iostream>

using namespace std;
using namespace FractureNetwork;

int main()
{
    Define_tol();

    string filepath = "DFN_files/FR3_data.txt";

    vector<Fracture> fratture;
    vector<Trace> tracce;

    bool done = ImportFractures(filepath, fratture);
    if(!done){
        cout << "Errore nell'apertura del file" << endl;
        return 1;
    }


    // //CHECK:
    // // Stampa il contenuto del vettore di fratture
    // cout << fratture.size() << endl;
    // cout << "Contenuto del vettore di fratture:\n";
    // for (const auto& f : fratture) {

    CalculateTraces(fratture,tracce);

    // //CHECK:
    // // Stampa il contenuto del vettore di tracce
    // cout << "Contenuto del vettore di tracce:\n";
    // for (const auto& traccia : tracce) {
    //     cout << fixed << traccia.vertices << endl << endl;
    // }

    SortingFractureTraces(fratture,tracce);

    // //CHECK:
    // // Stampa di fake per verificare lo spostamento
    // for (unsigned int i = 0; i<fratture.size();i++) {
    //     if(!fratture[i].passing.empty() || !fratture[i].not_passing.empty()){
    //         cout << i << " ";
    //         if(!fratture[i].passing.empty()){
    //             for(auto l: fratture[i].passing){
    //                 cout << l << " ";
    //             }
    //         }
    //         if(!fratture[i].not_passing.empty()){
    //             for(auto l: fratture[i].not_passing){
    //                 cout << l << " ";
    //             }
    //         }
    //         cout << endl;
    //     }
    // }


    //Print
    string traceFile = "tracce.txt";
    bool printed = PrintTrace(traceFile, tracce);
    if(!printed){
        cout << "Errore nell'apertura/creazione del file tracce.txt" << endl;
        return 1;
    }
    string tracceFile2 = "fratture_tracce.txt";
    bool printed2 = PrintFractureTraces(tracceFile2, fratture, tracce);
    if(!printed2){
         cout << "Errore nell'apertura/creazione del file fratture_tracce.txt " << endl;
        return 1;
    }

    //     //PUNTO 2:
    //     //
    //     //
    Mesh mesh;
    list<MatrixXd> cutted = cutting(fratture,tracce);
    extractinfo(cutted,mesh);

    cout << "okay" << endl;
    return 0;
}

