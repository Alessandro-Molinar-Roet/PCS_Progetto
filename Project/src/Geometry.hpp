#pragma once

#include"FracturesNetworkLibrary.hpp"
#include"tol.hpp"

namespace FractureNetwork{

// prende vettore di fratture e ne calcola tutte le intersezioni salvandole in vettore di tracce,
// inoltre per ogni fratuura salva il vetotre con le sue tracce
void CalculateTraces(vector<Fracture>& fractures, vector<Trace>& traces);

// Prende due frattura e calcola l se i loro bounding box si intersecano
bool near1(const MatrixXd& fracture1, const MatrixXd& fracture2);

// Funzione alternativa. Prende due fratture e calcola se le palle in cui sono inscritte si intersecano
bool near2(const MatrixXd& fracture1, const MatrixXd& fracture2);

// Prende una frattura e calcola la normale al piano che la contiene
Vector3d normalP(const MatrixXd& fracture);

// Prende due fratutre e restituisce true se si intersecano, inoltre se esiste calcola la loro intersezione e se è passante per le rispettive fratture
bool fracturesIntersection(const MatrixXd& fracture1,const MatrixXd& fracture2, Vector3d& E1, Vector3d& E2, bool& tips1, bool& tips2);

// Prende una fratuttra e una retta e ne calcola l'intersezione
bool lineFractIntersect(const MatrixXd& fracture, const Vector3d& pointOnLine,const Vector3d& direction, vector<Vector3d>& intersections);

Vector3d applyThreshold(const Vector3d& vec);

//******************************************************************************************************************************
//Parte due

void cutting(vector<Fracture>& fratture, vector<Trace>& tracce);

void split(const MatrixXd& polygon, const vector<unsigned int>& passing, const vector<Trace>& tracce, unsigned int counter, list<MatrixXd>& cutted);

bool IsInside(const MatrixXd& frattura, const Vector3d& tr1, const Vector3d& direzione, vector<Vector3d>& intersezioni, vector<unsigned int>& lato);

void extractinfo(const MatrixXd& polygon, Mesh& mesh, const unsigned int& counter);

}
