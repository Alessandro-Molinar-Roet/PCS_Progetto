// File contente strutture
#pragma once

#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;


namespace FractureNetwork{


// Struttura dati per salvare una frattura
//
// ID = Id poligono supposto uguale alla poszione in cui viene letto
// vertices = matrice contenente coordinate {[x1 x2 x3 x4...]; [y1 y2 y3 y4...]; [z1 z2 z3 z3...]}
// passing = vettore di Id delle tracce passanti
// not_passing = vettore di id tracce non passanti
struct Fracture{
    unsigned int ID = {};
    MatrixXd vertices {};

    vector<unsigned int> passing;
    vector<unsigned int> not_passing;
};


// Struttura dati per salvare una traccia
//
// ID = id traccia corrispondente a posizione in cui viene trovata
// length = lunghezza della traccia
// vertices = matrice contenente coordinate estremi {[x1 x2]; [y1 y2]; [z1 z2]}
// first_generator =  Id generatore 1
// second_generator = Id generatore 2
struct Trace{
    unsigned int ID = {};
    double length = {};
    Matrix<double, 3, 2> vertices = {};

    unsigned int first_generator = {};
    unsigned int second_generator = {};
};


// Struttura per salvare mesh generata dai tagli delle fratture lungo le tracce
//
// NumberCell0D = numero di punti della mesh
// Cell0DId = vettore di unsigned int con gli Id dei punti della mesh
// Cell0DCoordinates = vettore di Vector3d con le coordinate dei punti della mesh
// NumberCell1D = numero di lati della mesh
// Cell1DId = vector di unsigned int con gli Id dei lati della mesh
// Cell1DVertices = vector di Vector2i con gli Id dei due estremi di ongi lato
// NumberCell2D = numero di aree della mesh
//
// Cell2DId = vector di unsinged int con gli ID delle aree della mesh
// Cell1DVertices = vector di vector di unsigned int con gli Id degli estremi di ongi area
// Cell2DEdges = vector di vector di unsigned int con gli Id dei lati di ongi area
struct Mesh
{
    unsigned int NumberCell0D = 0;
    vector<unsigned int> Cell0DId = {};
    vector<Vector3d> Cell0DCoordinates = {};

    unsigned int NumberCell1D = 0;
    vector<unsigned int> Cell1DId = {};
    vector<Vector2i> Cell1DVertices = {};

    unsigned int NumberCell2D = 0;
    vector<unsigned int> Cell2DId = {};
    vector<vector<unsigned int>> Cell2DVertices = {};
    vector<vector<unsigned int>> Cell2DEdges = {};
};


// Funzione di hash per Vector3d.
//
//vec = Vector3d
// return -> il valore di hash associato al vec
struct Vector3dHash {
    size_t operator()(const Vector3d& vec) const {
        size_t hx = hash<double>()(vec.x());
        size_t hy = hash<double>()(vec.y());
        size_t hz = hash<double>()(vec.z());
        return hx ^ (hy << 1) ^ (hz << 2) ;  // << sposta il numero di bit a sinsitra
    }
};


// Funzione di hash per pair di unsigned int
//
// p = pair di unsigned int
// return -> il valore di hash associato alla pair
struct pair_hash {
    size_t operator()(const pair<unsigned int, unsigned int>& p) const {
        size_t p1 = hash<unsigned int>()(p.first);
        size_t p2 = hash<unsigned int>()(p.second);
        return p1^p2 + p2^p1;
    }
};




}
