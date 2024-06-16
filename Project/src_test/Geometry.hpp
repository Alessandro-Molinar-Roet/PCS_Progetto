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
bool lineFractIntersect(const MatrixXd& fracture, const Vector3d& pointOnLine,const Vector3d& direction,vector<Vector3d>& intersections);

Vector3d applyThreshold(const Vector3d& vec);

















// Serve Per punto due ?
// trova la posizione relativa di un punto e un poligono che giacciono sullo stesso piano ( 0 se all'esterno, 1 all'intenro, 2 sulla frontiera)
unsigned int IsInside(const MatrixXd& polygon, Vector3d normal, Vector3d point);


}
