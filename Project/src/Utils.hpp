#pragma once

#include"FracturesNetworkLibrary.hpp"

using namespace std;

namespace FractureNetwork {

// legge contenuto del file e salva in vettore di fratture
bool ImportFractures(const string& filepath, vector<Fracture>& fractures);

// ordina per lunghezza le traccie passanti e nont per ogni fratture
void SortingFractureTraces(vector<Fracture>& fractures, vector<Trace>& traces);

// stampa su file traccee
bool PrintTrace(const string& filepath, const vector<Trace>& traces);

// stampa su file frattura con relative tracce
bool PrintFractureTraces(const string& filepath, const vector<Fracture>& fractures, const vector<Trace>& traces);

}


