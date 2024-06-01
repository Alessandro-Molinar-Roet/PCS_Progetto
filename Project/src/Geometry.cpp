#include "Geometry.hpp"
#include "math.h"
#include "tol.hpp"
#include <iostream>

namespace FractureNetwork {

//calcola tutte le traccie e se sono passatnti o non passatni e le salva DOVE E COME??????
void CalculateTraces(vector<Fracture>& fratture, vector<Trace>& tracce){
    unsigned int counter = 0; //tiene il conto tracce trovate (per ID)
    unsigned int dim = fratture.size(); //numero di fratture

    list<Trace> temp_tracce = {}; // lista temporanea per salvare tracce
    vector<list<unsigned int>>  temp_passing(dim); // posizione i = lista id tracce passanti nel poligono i
    vector<list<unsigned int>>  temp_not_passing(dim); // posizione i = lista id tracce non passanti nel poligono i

    for(unsigned int i = 0; i<dim; i++){ //ciclo su tutte le matrici (poligoni)
        MatrixXd first_polygon = fratture[i].vertices;
        for(unsigned int j = i+1; j<dim; j++ ){ // ciclo su tutti i poligoni >i (altre coppie gia cocntrollate)
            MatrixXd second_polygon = fratture[j].vertices;
            if(near1(first_polygon, second_polygon)){ //se molto lontani non controllo nemmeno intersezione
                //inizializzo variabili utili
                Vector3d E1 = {};
                Vector3d E2 = {};
                bool tips1 = true;
                bool tips2 = true;
                bool find = TracciaTraPoligoni(first_polygon, second_polygon, E1, E2,tips1, tips2); //calcolo intersezioni
                if(find){
                    //Salvo variabili in struttura adeguata
                    Trace traccia = {};

                    Matrix<double,3,2> vertices = {};
                    vertices.col(0) = E1;
                    vertices.col(1) = E2;
                    traccia.vertices = vertices;

                    traccia.ID = counter; //ID della traccia = numero in cui è stata trovata
                    traccia.first_generator = fratture[i].ID;
                    traccia.second_generator = fratture[j].ID;
                    double length = sqrt( pow((vertices(0,0)-vertices(0,1)),2) +  pow((vertices(1,0)-vertices(1,1)),2) +  pow((vertices(2,0)-vertices(2,1)),2)); //ATTENZIONE berrone no
                    traccia.length = length;

                    temp_tracce.push_back(traccia);

                    //ATTENZIONE puschback diretto su vettore ridimensionando vettore???
                    if(!tips1){
                        temp_passing[i].push_back(counter);
                    }
                    else{
                        temp_not_passing[i].push_back(counter);
                    }
                    if(!tips2){
                        temp_passing[j].push_back(counter);
                    }
                    else{
                        temp_not_passing[j].push_back(counter);
                    }

                    counter ++;  
                }
            }
        }
    }

    //ATTENZIONE Vicini resize più sicuro ma si rompe con resize
    //ATTENZIONE serve sapere cosa fa move?
    tracce.reserve(size(temp_tracce));
    move(begin(temp_tracce), end(temp_tracce), back_inserter(tracce));

    for(unsigned int i = 0; i<dim; i++){
        fratture[i].passing.reserve(temp_passing[i].size());
        move(begin(temp_passing[i]), end(temp_passing[i]), back_inserter(fratture[i].passing));
        fratture[i].not_passing.reserve(temp_not_passing[i].size());
        move(begin(temp_not_passing[i]), end(temp_not_passing[i]), back_inserter(fratture[i].not_passing));

    }
}

//***************************************************************************************************************
Vector3d normaleP(const MatrixXd& frattura)
{
    Vector3d v1 = frattura.col(1)-frattura.col(0);
    Vector3d v2 = frattura.col(2)-frattura.col(1);
    return v1.cross(v2);
}

//***************************************************************************************************************
bool intersRettaPoly(const MatrixXd& frattura, const Vector3d& puntoRetta,const Vector3d& direzione,vector<Vector3d>& intersezioni)
{

    for(unsigned int i=0;i<frattura.cols();i++){
        //intersezione lato-retta
        Vector3d puntoipoly = frattura.col(i);
        Vector3d latoipoly =frattura.col((i+1)%frattura.cols())-frattura.col(i);
        double t = ((puntoRetta.cross(direzione)).dot((latoipoly.cross(direzione)))-(puntoipoly.cross(direzione)).dot((latoipoly.cross(direzione))))/pow((latoipoly.cross(direzione)).norm(), 2);
        if (t >= 0.0 && t <= 1.0)
        {
            intersezioni.push_back(puntoipoly + t * latoipoly);
        }
        if(intersezioni.size()==2){
            return true;
        }
    }
    return false;
}

//***************************************************************************************************************
Vector3d applyThreshold(const Vector3d& vec) {
    Vector3d result = vec; // Create a copy of the input vector
    for (int i = 0; i < result.size(); i++) {
        if (std::abs(result[i]) < tol) {
            result[i] = 0.0;
        }
    }
    return result;
}

//***************************************************************************************************************
bool TracciaTraPoligoni(const MatrixXd& frattura1,const MatrixXd& frattura2, Vector3d& E1, Vector3d& E2, bool& tips1, bool& tips2)
{
    //Considero la retta di intersezione tra i due piani che contengono i poligoni
    //Calcolo la direzione della retta come prodotto vettoriale tra le normali ai piani
    Vector3d n1 = normaleP(frattura1);
    Vector3d n2 = normaleP(frattura2);
    Vector3d direzione = n1.cross(n2).normalized();
    //Se il prodotto vettoriale è nullo, i piani sono paralleli
    if(abs(direzione.norm()) < tol){ //abs<tol
        return false;
    }
    //Interseco il piano del poligono 1 con la retta che contiene il lato 1 del poligono 2 per trovare un punto della retta
    Vector3d lato = frattura2.col(1)-frattura2.col(0);
    Vector3d punto = frattura2.col(0);
    //Controllo che quel lato non sia parallelo al piano, e se lo è, prendo il lato successivo
    if(abs(lato.dot(n1))<tol)
    {
        lato = frattura2.col(2)-frattura2.col(1);
        punto = frattura2.col(1);
    }
    //Soluzione dell'intersezione piano-retta
    double t = -(n1.dot(punto-frattura1.col(0)))/(n1.dot(lato));
    Vector3d puntoRetta = punto+t*lato;
    //Ora controlliamo se il poligono interseca la retta

    vector<Vector3d> intersezioni;
    intersezioni.reserve(4);
    vector<Vector3d> inter1;
    if(intersRettaPoly(frattura1,puntoRetta,direzione,inter1)){
        intersezioni.push_back(applyThreshold(inter1[0]));
        intersezioni.push_back(applyThreshold(inter1[1]));
    }
    else{
        return false;
    }
    vector<Vector3d> inter2;
    if(intersRettaPoly(frattura2,puntoRetta,direzione,inter2)){
        intersezioni.push_back(applyThreshold(inter2[0]));
        intersezioni.push_back(applyThreshold(inter2[1]));
    }
    else{
        return false;
    }
    //Identifico i 4 punti trovati sulla retta
    double a = ((intersezioni[0]-puntoRetta).dot(direzione))/(direzione.norm()*direzione.norm());
    double b = ((intersezioni[1]-puntoRetta).dot(direzione))/(direzione.norm()*direzione.norm());
    //a e b sono delle posizioni relative a puntoRetta, ma non so come sono ordinati perché questo dipende
    //dal senso di lettura dei lati dei poligoni e dalla posizione relativa alla retta
    if(a>b)
    {
        double atemp = a;
        Vector3d temp = intersezioni[0];
        a=b;
        b=atemp;
        intersezioni[0] = intersezioni[1];
        intersezioni[1] = temp;

    }
    double c = ((intersezioni[2]-puntoRetta).dot(direzione))/(direzione.norm()*direzione.norm());
    double d = ((intersezioni[3]-puntoRetta).dot(direzione))/(direzione.norm()*direzione.norm());

    if(c>d)
    {
        double ctemp =c;
        Vector3d temp = intersezioni[2];
        c=d;
        d=ctemp;
        intersezioni[2] = intersezioni[3];
        intersezioni[3] = temp;
    }

    if(b<=c || d<=a)// i segmenti [a,b] e [c,d] non si intersecano
    {
        return false;
    }
    //c'è un'intersezione "propria", le tracce sono non-passanti per entrambe le fratture
    if(a<c && c<b && b<d)
    {
        E1 = intersezioni[2];
        E2 = intersezioni[1];
        tips1 = true;
        tips2 = true;
        return true;
    }
    else if(c<a && a<d && d<b)
    {
        E1 = intersezioni[0];
        E2 = intersezioni[3];
        tips1 = true;
        tips2 = true;
        return true;
    }
    //un segmento è contenuto nell'altro, possibilmente con uno o entrambi gli estremi che si sovrappongono
    else if(a<=c && d<=b)
    {
        E1 = intersezioni[2];
        E2 = intersezioni[3];
        tips2 = false;
        tips1 = true;
        if(abs(a-c)<tol && abs(d-b)<tol)
        {
            tips1 = false;
        }
        return true;
    }
    else if(c<=a && b<=d)
    {
        E1 = intersezioni[0];
        E2 = intersezioni[1];
        tips1 = false;
        tips2 = true;
        if(abs(a-c)<tol && abs(d-b)<tol)
        {
            tips2 = false;
        }
        return true;
    }
    return false;
}


//***************************************************************************************************************
// Bounding box:
bool near1(const MatrixXd& first_polygon, const MatrixXd& second_polygon)
{
    bool controllo = true;

    Vector3d max1 = first_polygon.rowwise().maxCoeff();
    Vector3d min1 = first_polygon.rowwise().minCoeff();
    Vector3d max2 = second_polygon.rowwise().maxCoeff();
    Vector3d min2 = second_polygon.rowwise().minCoeff();

    if (min1(0) > max2(0) || min1(1) > max2(1) || min1(2) > max2(2))
    {
        controllo = false;
    }
    if (min2(0) > max1(0) || min2(1) > max1(1) || min2(2) > max1(2))
    {
        controllo = false;
    }

    return controllo;
}

//***************************************************************************************************************
// Calcolo punto medio del primo e del secondo poligono
// Calcolo la distanza tra essi e controllo che sia minore di una certa tolleranza
// 2 casi:
// se è <, allora c'è la possibilità che si intersechino
// se è >, allora non si intersecano
bool near2(const MatrixXd& first_polygon, const MatrixXd& second_polygon){
    bool controllo = true;
    VectorXd bar_polygon1 = first_polygon.colwise().mean();
    VectorXd bar_polygon2 = second_polygon.colwise().mean();

    double max_1 = 0.0;
    for (unsigned int j = 0; j < first_polygon.col(0).size(); j++){
        double length1 = bar_polygon1.norm() - first_polygon.col(j).norm();
        if (length1 > max_1){
            max_1 = length1;
        }
    }

    double max_2 = 0.0;
    for (unsigned int j = 0; j < second_polygon.col(0).size(); j++){
        double length2 = bar_polygon2.norm() - second_polygon.col(j).norm();
        if (length2 > max_2){
            max_2 = length2;
        }
    }

    double control_length = max_1 + max_2;
    double bar_distance = (bar_polygon1 - bar_polygon2).norm();

    if (bar_distance < control_length){
        cout << "Intersection is possible" << endl;
    }
    else{
        controllo = false;
        cout << "Intersection is impossible" << endl;
    }
    return controllo;
}




//***************************************************************************************************************
//Serve per punto 2?
unsigned int IsInside(const MatrixXd& polygon, Vector3d normal, Vector3d point){
    unsigned int counter = 1; //ATTENZIONE essite bolean a 3?
    for (int i = 0; i < polygon.row(0).size(); ++i) { // cicla su tutti i lati di un poligono (colonne matrice)
        Vector3d lato = polygon.col((i + 1) % polygon.row(0).size()) - polygon.col(i); // vertice[i]-vertice[i+1]
        Vector3d toPoint = point - polygon.col(i); // punto-vertice
        double prodotto_misto = normal.dot(lato.cross(toPoint)); //prodotto misto normale*(lato x toPoint)
        if (lato.cross(toPoint).norm()<tol){ // ATTENZIONE Tolleranza se no si rompe dove la definsico?????
            counter = 2; // sul bordo
        }
        if (prodotto_misto < 0) {
            counter = 0; // fuori
            return counter;
        }
    }
    return counter; //interno
}

}

