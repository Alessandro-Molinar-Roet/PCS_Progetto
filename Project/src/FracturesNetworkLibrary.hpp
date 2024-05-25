// File contente strutture
#pragma once

#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

namespace FractureNetwork{

//struttura dati per salvare le fratture
//num = numero di fratture totali
//f_ID = vettore di ID
//f_Vertices = vettore di matrici contenenti coordinate {[x1 x2 x3 x4...]; [y1 y2 y3 y4...]; [z1 z2 z3 z3...]}
struct Fractures{
    unsigned int num = 0;
    vector<unsigned int> f_ID = {};
    vector<MatrixXd> f_Vertices {}; // ATTENZIONE arrotonda numeri Errore o va bene? //ATTENZIONE 3xn
    //ATTENZIONE serve ridefinire operatore = >= e <= per includere tolleranza ovunque? per farlo devo ridefiniri sturuttura vettori?
};

//struttura dati per salvare le tracce
//t_ID = vettore di ID
//t_length = vettore di lunghezza
//t_Vertices = vettore di matrici con vertici su colonne
//t_generator = vettore di Id di poligoni generatori (2i,2i-1)
struct Traces{
    unsigned int num = 0;
    vector<unsigned int> t_ID = {}; // ATTENZIONE meglio lista non so quante intersezioni ho? ma devo accedervi spesso creo il doppio?
    vector<double> t_length = {};
    vector<Matrix<double, 3, 2>> t_Vertices = {};
    vector<unsigned int> t_generator = {}; // ATTENZIONE come salvo gli ID dei piani che generano uno frattura? posizione = 2i+1 e 2i+2?

    //ATTENZIONE struttura dati adatta? logicamente devo metterlo qui o a parte?
    map<unsigned int,  vector<vector<unsigned int>>> Dfn = {};
};


}
