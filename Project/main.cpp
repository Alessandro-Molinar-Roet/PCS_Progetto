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
    MatrixXd frattura1(3, 4);
    frattura1 << 0.0, 1.0, 1.0, 0.0,
        0.0, 0.0, 1.0, 1.0,
        0.0, 0.0, 0.0, 0.0;

    MatrixXd frattura2(3, 4);
    frattura2 << 0.8, 0.8, 0.8, 0.8,
        0.0, 0.0, 1.0, 1.0,
        -0.1, 0.3, 0.3, -0.1;

    MatrixXd frattura3(3, 4);
    frattura3 << -0.237778, 0.3161837, 0.3161837, -0.237778,
        0.5, 0.5, 0.5, 0.5,
        -0.34444, -0.34444, 0.4528389, 0.4528389;

    Vector3d E1, E2;
    bool tips1, tips2;

    if (TracciaTraPoligoni(frattura1, frattura2, E1, E2, tips1, tips2)) {
        // Se la funzione ritorna true, stampa i vettori E1 ed E2
        cout << "E1: " << E1.transpose() << endl;
        cout << "E2: " << E2.transpose() << endl;
        cout << "tips1: "<<tips1<<endl;
    } else {
        cout << "Nessuna traccia tra i poligoni trovata." << endl;
    }



  return 0;
}
