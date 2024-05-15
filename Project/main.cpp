//Librerie
#include "FracturesNetworkLibrary.hpp"
#include "Utils.hpp"
#include "Geometry.hpp"
#include <iostream>
#include <iomanip>

using namespace std;
using namespace FractureNetwork;

int main()
{
    string filepath = "DFN_files/FR3_data.txt";
    Fractures fratture;
    Traces tracce;

    bool done = ImportFractures(filepath, fratture);
    if(!done){
        cout << "Errore nell'importazione del file" << endl;
        return 1;
    }

 /*
    //CHECK:
    // Stampa il contenuto del vettore di fratture
    cout << fratture.num << endl;
    cout << "Contenuto del vettore di fratture:\n";
    for (const auto& matrice : fratture.f_Vertices) {
        cout << setprecision(10) << fixed <<matrice << endl << endl;
    }
*/
    CalculateFracture(fratture,tracce);
    SortingFractureTraces(tracce);

/*
    //CHECK:
    // Stampa il contenuto del vettore di tracce
    cout << "Contenuto del vettore di tracce:\n";
    for (const auto& matrice : tracce.t_Vertices) {
        cout << setprecision(10) << fixed << matrice << endl << endl;
    }

 */

    string traceFile = "tracce.txt";
    bool printed = PrintTrace(traceFile, tracce);
    if(!printed){
        return 1;
    }
    string tracceFile2 = "fratture_tracce.txt";
    bool printed2 = PrintFractureTraces(tracceFile2, tracce);
    if(!printed2){
        return 1;
    }

//     //PUNTO 2:
//     //
//     //

cout << "okay" << endl;
  return 0;
}
