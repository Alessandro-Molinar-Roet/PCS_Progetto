#pragma once

#include"FracturesNetworkLibrary.hpp"

namespace FractureNetwork{

//prende fratture lette e calcola le intersezioni salvandole nelle tracce
void CalculateTraces(const Fractures& fratture, Traces& tracce);// non voglio modificare fratture o si? const & vs const & const

//prende due matrici contenenti vertici poligono e controlla se potrebbero avere intersezioni
bool near(const MatrixXd& first_polygon, const MatrixXd& secodn_polygon);

//Calcola la normale al piano che contiene la frattura
Vector3d normaleP(const MatrixXd& frattura);

//Trova se c'è interesezione tra due fratture, nel caso calcola gli estremi del segmento e dice se la traccia è passante o meno
bool TracciaTraPoligoni(const MatrixXd& frattura1,const MatrixXd& frattura2, Vector3d& E1, Vector3d& E2, bool& tips1, bool& tips2);




















// Serve Per punto due ?
// trova la posizione relativa di un punto e un poligono che giacciono sullo stesso piano ( 0 se all'esterno, 1 all'intenro, 2 sulla frontiera)
unsigned int IsInside(const MatrixXd& polygon, Vector3d normal, Vector3d point);


}
