// File contente strutture
#pragma once

#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

namespace FractureNetwork{


struct Fracture{
    // struttura dati per salvare una frattura
    // ID = Id poligono supposto uguale alla poszione in cui viene letto
    // vertices =  matrice contenente coordinate {[x1 x2 x3 x4...]; [y1 y2 y3 y4...]; [z1 z2 z3 z3...]}
    // passing = vettore di id tracce passanti
    // not_passing = vettore di id tracce non passanti

    unsigned int num = 0;
    unsigned int ID = {};
    MatrixXd vertices {}; // ATTENZIONE arrotonda numeri //ATTENZIONE 3xn

    vector<unsigned int> passing;
    vector<unsigned int> not_passing;
};


struct Trace{
    // struttura dati per salvare una traccia
    // ID = id traccia corrispondente a posizione in cui viene trovata
    // length = lunghezza della traccia
    // vertices = matrice contenente coordinate estremi {[x1 x2]; [y1 y2]; [z1 z2]}
    // first_generator =  Id generatore 1
    // second_generator = Id generatore 2

    unsigned int ID = {};
    double length = {};
    Matrix<double, 3, 2> vertices = {};

    unsigned int first_generator = {};
    unsigned int second_generator = {};
};


struct Mesh
{
    // Struttura per salvare mesh generata dai tagli delle fratture
    // Celle0D = punti (numero, ID, coordinate)
    // Cell01D = lati (numero, ID, ID estremi lato)
    // Celle2D = aree (numero, ID, ID vertici, ID lati)

    unsigned int NumberCell0D = 0;
    vector<unsigned int> Cell0DId = {};
    vector<Vector2d> Cell0DCoordinates = {};

    unsigned int NumberCell1D = 0;
    vector<unsigned int> Cell1DId = {};
    vector<Vector2i> Cell1DVertices = {};

    unsigned int NumberCell2D = 0;
    vector<unsigned int> Cell2DId = {};
    vector<vector<unsigned int>> Cell2DVertices = {};
    vector<vector<unsigned int>> Cell2DEdges = {};
};


}
