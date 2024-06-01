#pragma once

#include"FracturesNetworkLibrary.hpp"
#include"tol.hpp"

namespace FractureNetwork{

// prende vettore di fratture e ne calcola tutte le intersezioni salvandole in vettore di tracce,
// inoltre per ogni fratuura salva il vetotre con le sue tracce
void CalculateTraces(vector<Fracture>& fratture, vector<Trace>& tracce);

// Prende due frattura e calcola l se i loro bounding box si intersecano
bool near1(const MatrixXd& first_polygon, const MatrixXd& secodn_polygon);

// Funzione alternativa. Prende due fratture e calcola se le palle in cui sono inscritte si intersecano
bool near2(const MatrixXd& first_polygon, const MatrixXd& secodn_polygon);

// Prende una frattura e calcola la normale al piano che la contiene
Vector3d normaleP(const MatrixXd& frattura);

// Prende due fratutre e restituisce true se si intersecano, inoltre se esiste calcola la loro intersezione e se è passante per le rispettive fratture
bool TracciaTraPoligoni(const MatrixXd& frattura1,const MatrixXd& frattura2, Vector3d& E1, Vector3d& E2, bool& tips1, bool& tips2);

// Prende una fratuttra e una retta e ne calcola l'intersezione
bool intersRettaPoly(const MatrixXd& frattura, const Vector3d& puntoRetta,const Vector3d& direzione,vector<Vector3d>& intersezioni);

Vector3d applyThreshold(const Vector3d& vec);

















// Serve Per punto due ?
// trova la posizione relativa di un punto e un poligono che giacciono sullo stesso piano ( 0 se all'esterno, 1 all'intenro, 2 sulla frontiera)
unsigned int IsInside(const MatrixXd& polygon, Vector3d normal, Vector3d point);


}
