#pragma once

#include"FracturesNetworkLibrary.hpp"

namespace FractureNetwork{

//prende fratture lette e calcola le intersezioni salvandole nelle tracce
void CalculateFracture(const Fractures& fratture, Traces& tracce);// non voglio modificare fratture o si? const & vs const & const

//prende due matrici contenenti vertici poligono e controlla se potrebbero avere intersezioni
bool near(const MatrixXd& first_polygon, const MatrixXd& secodn_polygon);

//trova l'intersezione tra i lati di first_polygon e secodn_polygon
bool FindIntersection(const MatrixXd& first_polygon, const MatrixXd& second_polyogn, Matrix<double, 3, 2>& vertices, bool& find, bool& complete, bool& tips2);

//trova la posizione relativa di un punto e un poligono che giacciono sullo stesso piano ( 0 se all'esterno, 1 all'intenro, 2 sulla frontiera)
unsigned int IsInside(const MatrixXd& polygon, Vector3d normal, Vector3d point);

//ordina le traccie della fratture
void SortingFractureTraces(const Fractures& fratture, Traces& tracce);

Vector3d normaleP(const MatrixXd& frattura);

bool TracciaTraPoligoni(const MatrixXd& frattura1,const MatrixXd& frattura2, Vector3d& E1, Vector3d& E2, bool& tips1, bool& tips2);

}
