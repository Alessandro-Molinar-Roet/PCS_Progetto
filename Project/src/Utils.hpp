#pragma once

#include"FracturesNetworkLibrary.hpp"

using namespace std;

namespace FractureNetwork {


// Funzione per definire tolleranze
//
// Dopo aver chiesto all'utente le tolleranze mono e bi-dimensionali controlla che queste non siano inferiori ad una tolleranza standard
// definita a priori e le salva in due extern const
void Define_tol();


// Funzione per leggere il file e salvare i dati in un vettore di fratture
//
// filepath = percorso dove si trova il file da leggere
// fractures = vettore passato per referenza in cui vengono salvate le fratture
// return -> un booleano che indica se il file è stato aperto correttamente
bool ImportFractures(const string& filepath, vector<Fracture>& fractures);


// Funzione per odinare le tracce passanti e non di ogni frattura secondo la loro lunghezza
//
// filepath = percorso dove si trova il file da leggere
// fractures = vettore passato per referenza in cui vengono salvate le fratture
// return -> un booleano che indica se il file è stato aperto correttamente
void SortingFractureTraces(vector<Fracture>& fractures, vector<Trace>& traces);


// Funzione per stampare il vettore di tracce su un file
//
// filepath = nome del file su cui stampare
// traces = vettore di tracce da stampare
// return -> booleano che indica se la creazione/ apertura del file è andata a buon fine
bool PrintTrace(const string& file_name, const vector<Trace>& traces);


// Funzione per stampare le tracce associate ad ogni frattura
//
// file_name = nome del file su cui stampare
// fractures = vettore di fratture
// traces = vettore di tracce
// return -> booleano che indica se la creazione/ apertura del file è andata a buon fine
bool PrintFractureTraces(const string& file_name, const vector<Fracture>& fractures, const vector<Trace>& traces);


}


