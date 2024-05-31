#pragma once

#include"FracturesNetworkLibrary.hpp"

using namespace std;

namespace FractureNetwork {

// legge contenuto del file e salva in vettore di fratture
bool ImportFractures(const string& filepath, vector<Fracture>& fratture);

// ordina per lunghezza le traccie passanti e nont per ogni fratture
void SortingFractureTraces(vector<Fracture>& fratture, vector<Trace>& tracce);

// stampa su file traccee
bool PrintTrace(const string& filepath, const vector<Trace>& tracce);

// stampa su file frattura con relative tracce
bool PrintFractureTraces(const string& filepath, const vector<Fracture>& fratture, const vector<Trace>& tracce);

}


