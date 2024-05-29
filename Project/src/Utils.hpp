#pragma once

#include<iostream>
#include"FracturesNetworkLibrary.hpp"

using namespace std;

namespace FractureNetwork {

//legge contenuto del file e salva in struttura fratture
bool ImportFractures(const string& filepath, Fractures& fratture);

//***************************************************************************************************************

//prende fratture lette e calcola le intersezioni salvandole nelle tracce
void CalculateFracture(const Fractures& fratture, Traces& tracce);// non voglio modificare fratture o si? const & vs const & const

//prende due matrici contenenti vertici poligono e controlla se potrebbero avere intersezioni
bool near1(const MatrixXd& first_polygon, const MatrixXd& secodn_polygon);

bool near2(const MatrixXd& first_polygon, const MatrixXd& secodn_polygon);

//trova l'intersezione tra i lati di first_polygon e secodn_polygon
void FindIntersection(const MatrixXd& first_polygon, const MatrixXd& second_polyogn, Matrix<double, 3, 2>& vertices, bool& find, bool& complete);

//trova la posizione relativa di un punto e un poligono che giacciono sullo stesso piano ( 0 se all'esterno, 1 all'intenro, 2 sulla frontiera)
unsigned int IsInside(const MatrixXd& polygon, Vector3d normal, Vector3d point);


}
