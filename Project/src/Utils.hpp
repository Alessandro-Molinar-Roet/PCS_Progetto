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

//stampa su file frattura con relative traccie ordinate (tips-lunghezza)
bool PrintFractureTraces(const string& filepath, const Traces& tracce);

}


