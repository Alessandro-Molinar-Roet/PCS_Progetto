#pragma once

#include"FracturesNetworkLibrary.hpp"

namespace FractureNetwork{

//PARTE 1

// Funzione per calcolare tutte le tracce date dall'intersezioni tra fratture
//
//fractures = vettore di fratture, passato per referenza per andare a inserie gli ID delle tracce che contiene
//traces = vettore di tracce, passato per refernza per riempirlo con le tracce trovate
void CalculateTraces(vector<Fracture>& fractures, vector<Trace>& traces);


// Funzione per escludere alcuni casi in cui due fratture molto lontane non possono intersecarsi,
// sfrutta metodo dei bounding box
//
// fracture1 = MatrixXd conenente i vertici della frattura numero 1
// fracture2 = MatrixXd conenente i vertici della frattura frattura numero 2
//
// return -> un boleano che indica se l'intersezione è possbile o meno
bool near1(const MatrixXd& fracture1, const MatrixXd& fracture2);


// Funzione alternativa a near 1 (NON utilizzata, solo per confronto empirico)
// sfrutta metodo delle circoferenze inscirtte
//
// fracture1 = MatrixXd conenente i vertici della frattura numero 1
// fracture2 = MatrixXd conenente i vertici della frattura numero 2
//
// return -> un boleano che indica se l'intersezione è possbile o meno
bool near2(const MatrixXd& fracture1, const MatrixXd& fracture2);


// Funzione che calcola la normale al piano contenente una frattura
//
// fracture = MatrixXd conenente i vertici di una frattura di cui voglio calcolare la normale
//
// return -> Vector3d contenente la normale
Vector3d normalP(const MatrixXd& fracture);


// Funzione che calcola se esistono i punti di intersezione tra due fratture
//
// fracture1 = MatrixXd conenente i vertici di una frattura numero 1
// fracture2 = MatrixXd conenente i vertici di una frattura numero 2
// E1 = Vector3d contenete, se esiste, un estremo del segmento di intersezione delle due fratture
// E2 = Vector3d contenete, se esiste, un estremo del segmento di intersezione delle due fratture
// tips1 = booleano che indica se la traccia trovata è non passante per la frattura numero 1
// tips2 = booleano che indica se la traccia trovata è non passante per la frattura numero 2
//
// return -> boolenao che indica se esiste la traccia
bool fracturesIntersection(const MatrixXd& fracture1,const MatrixXd& fracture2, Vector3d& E1, Vector3d& E2, bool& tips1, bool& tips2);


// Funzione che calcola l'intersezione tra una frattura e una retta
//
// fracture = MatrixXd conenente i vertici di una frattura
// pointOnLine = Vector3d contenente un punto appartenente alla retta
// direction = Vector3d contenente la direzione della retta
// intersections = Vector3d contenente i punti di intersezione tra i lati del poligono e la retta
//
// return -> booleano che indica se la retta si interseca con la frattura
bool lineFractIntersect(const MatrixXd& fracture, const Vector3d& pointOnLine,const Vector3d& direction, vector<Vector3d>& intersections);


//******************************************************************************************************************************
// PARTE 2

// Funzione che taglia tutte le fratture lungo i prolungamenti di tutte le tracce che contiene
//
// fratture = vettore di fratture da tagliare
// tracce = vettore di tracce con le inforamzione necessarie per tagliare le fratture
list<MatrixXd> cutting(vector<Fracture>& fratture, vector<Trace>& tracce);


// Funzione ricorsvia che divide taglia una singola frattura lungo una traccia (ricorsivamente taglia lungo tutte le tracce)
//
// polygon = MatrixXd contenente i vertici della frattura
// all = vettore di unsigned int con gli ID di tutte le tracce contenute dalla frattura
// tracce = vettore di tracce
// counter = contatore che segna per quale traccia tagliare alla chiamta ricorsiva attuale
// cutted = lista di MatrixXd in cui salvare tutti i vertici dei sotto poligoni generati dal taglio
void split(const MatrixXd& polygon, const vector<unsigned int>& all, const vector<Trace>& tracce, unsigned int counter, list<MatrixXd>& cutted);


// Funzione controlla se una traccia e contenuta in un dato poligono e allo stesso tempo calcola l'intersezioni e i lati su cui essa avviene
//
// frattura = MatrixXd con vertici di un poligono
// tr1 = Vectro3d con uno dei due estremi della traccia
// direzione = Vector3d con la direzione della traccia
// intersezioni = vector<Vector3d>> passato per refenza per salvare punti di intersezione
// lati = vector di unsigned int che passato per refenza per salvare i lati del poligono su cui avvengono le intersezioni
//
// return -> booleano che indica se almeno un estremo della traccia è contenuto nel poligono
bool IsInside(const MatrixXd& frattura, const Vector3d& tr1, const Vector3d& direzione, vector<Vector3d>& intersezioni, vector<unsigned int>& lati);


// Funzione che salva i sotto poligoni generati dal taglio in una mesh
//
// cutted = lista di MatrixXd contenenti i vertici di poligoni
// mesh = struttura in cui salvare le informazioni
void extractinfo(const list<MatrixXd>& cutted, Mesh& mesh);

}
