//Librerie
#include "FracturesNetworkLibrary.hpp"
#include "Utils.hpp"
#include "Geometry.hpp"
#include <iostream>

using namespace std;
using namespace FractureNetwork;

int main()
{
    //PARTE 1:
    Define_tol();

    string filepath = "DFN_files/FR3_data.txt";

    //Inizializzazione delle strutture:
    vector<Fracture> fratture = {};
    vector<Trace> tracce = {};

    //Letture da file:
    bool done = ImportFractures(filepath, fratture);
    if(!done){
        cerr << "Errore nell'apertura del file: " << filepath << endl;
        return 1;
    }

    //Calcolo delle tracce:
    CalculateTraces(fratture,tracce);

    //Ordino tracce di ogni frattura:
    SortingFractureTraces(fratture,tracce);

    //Print su file:
    string traceFile = "tracce.txt";
    bool printed = PrintTrace(traceFile, tracce);
    if(!printed){
        cerr << "Errore nell'apertura/creazione del file tracce.txt" << endl;
        return 1;
    }

    string tracceFile2 = "fratture_tracce.txt";
    bool printed2 = PrintFractureTraces(tracceFile2, fratture, tracce);
    if(!printed2){
         cout << "Errore nell'apertura/creazione del file fratture_tracce.txt " << endl;
        return 1;
    }

    //PARTE 2:
    //Inizializzo struttura:
    Mesh mesh = {};

    //Taglio tutte le fratture:
    list<MatrixXd> cutted = cutting(fratture, tracce);

    //Converto e salvo lista di matrici in una mesh
    extractinfo(cutted, mesh);

    return 0;
}

