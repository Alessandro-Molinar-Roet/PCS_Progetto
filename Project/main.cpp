//Librerie
#include"FracturesNetworkLibrary.hpp"
#include"Utils.hpp"
#include<iomanip>

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

    //CHECK:
    // Stampa il contenuto del vettore di fratture
    cout << fratture.num << endl;
    cout << "Contenuto del vettore di fratture:\n";
    for (const auto& matrice : fratture.f_Vertices) {
        cout << setprecision(10) << fixed <<matrice << endl << endl;
    }

    CalculateFracture(fratture,tracce);

    //CHECK:
    // Stampa il contenuto del vettore di traccie
    cout << "Contenuto del vettore di traccie:\n";
    for (const auto& matrice : tracce.t_Vertices) {
        cout << setprecision(10) << fixed << matrice << endl << endl;
    }

    //PUNTO 1:
    //prossimi passi:
    //stampare traccie PrintTrace
    //ordinare per stampare per fratture le sue traccie Sort e PrintDFN

    //PUNTO 2:
    //
    //

    cout << "okay" << endl;
  return 0;
}
