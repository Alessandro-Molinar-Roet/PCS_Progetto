#pragma once
#include<vector>
#include<iostream>
#include"FracturesNetworkLibrary.hpp"

using namespace std;

namespace FractureNetwork {

//legge contenuto del file e salva in struttura fratture
bool ImportFractures(const string& filepath, Fractures& fratture);

//ordina le traccie della fratture
void SortingFractureTraces(Traces& tracce);

//stampa su file traccie
bool PrintTrace(const string& filepath, const Traces& tracce);


//********************************************************************************************************************
Vector3d normaleP(const MatrixXd& frattura);

Vector4d piano(const MatrixXd& frattura);

bool TracciaTraPoligoni(const MatrixXd& frattura1,const MatrixXd& frattura2, Vector3d& E1, Vector3d& E2, bool& tips1, bool& tips2);


//stampa su file frattura con relative traccie ordinate (tips-lunghezza)
bool PrintFractureTraces(const string& filepath, const Traces& tracce);

}


