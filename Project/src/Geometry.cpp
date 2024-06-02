#include "Geometry.hpp"
#include "math.h"
#include "tol.hpp"
#include <iostream>

namespace FractureNetwork {

//calcola tutte le traccie e se sono passatnti o non passatni e le salva DOVE E COME??????
void CalculateTraces(vector<Fracture>& fractures, vector<Trace>& traces){
    unsigned int counter = 0; //tiene il conto tracce trovate (per ID)
    unsigned int dim = fractures.size(); //numero di fratture

    list<Trace> temp_traces = {}; // lista temporanea per salvare tracce
    vector<list<unsigned int>>  temp_passing(dim); // posizione i = lista id tracce passanti nel poligono i
    vector<list<unsigned int>>  temp_not_passing(dim); // posizione i = lista id tracce non passanti nel poligono i

    for(unsigned int i = 0; i<dim; i++){ //ciclo su tutte le matrici (poligoni)
        MatrixXd fracture1 = fractures[i].vertices;
        for(unsigned int j = i+1; j<dim; j++ ){ // ciclo su tutti i poligoni >i (altre coppie gia cocntrollate)
            MatrixXd fracture2 = fractures[j].vertices;
            if(near1(fracture1, fracture2)){ //se molto lontani non check nemmeno intersezione
                //inizializzo variabili utili
                Vector3d E1 = {};
                Vector3d E2 = {};
                bool tips1 = true;
                bool tips2 = true;
                bool find = fracturesIntersection(fracture1, fracture2, E1, E2,tips1, tips2); //calcolo intersezioni
                if(find){
                    //Salvo variabili in struttura adeguata
                    Trace traccia = {};

                    Matrix<double,3,2> vertices = {};
                    vertices.col(0) = E1;
                    vertices.col(1) = E2;
                    traccia.vertices = vertices;

                    traccia.ID = counter; //ID della traccia = numero in cui è stata trovata
                    traccia.first_generator = fractures[i].ID;
                    traccia.second_generator = fractures[j].ID;
                    double length = sqrt( pow((vertices(0,0)-vertices(0,1)),2) +  pow((vertices(1,0)-vertices(1,1)),2) +  pow((vertices(2,0)-vertices(2,1)),2)); //ATTENZIONE berrone no
                    traccia.length = length;

                    temp_traces.push_back(traccia);

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
    traces.reserve(size(temp_traces));
    move(begin(temp_traces), end(temp_traces), back_inserter(traces));

    for(unsigned int i = 0; i<dim; i++){
        fractures[i].passing.reserve(temp_passing[i].size());
        move(begin(temp_passing[i]), end(temp_passing[i]), back_inserter(fractures[i].passing));
        fractures[i].not_passing.reserve(temp_not_passing[i].size());
        move(begin(temp_not_passing[i]), end(temp_not_passing[i]), back_inserter(fractures[i].not_passing));

    }
}

//***************************************************************************************************************
Vector3d normalP(const MatrixXd& fracture)
{
    Vector3d v1 = fracture.col(1)-fracture.col(0);
    Vector3d v2 = fracture.col(2)-fracture.col(1);
    return v1.cross(v2);
}

//***************************************************************************************************************
bool lineFractIntersect(const MatrixXd& fracture, const Vector3d& pointOnLine,const Vector3d& direction,vector<Vector3d>& intersections)
{

    for(unsigned int i=0;i<fracture.cols();i++){
        //intersezione lato-retta
        Vector3d ithFractureVertex = fracture.col(i);
        Vector3d ithFractureEdge =fracture.col((i+1)%fracture.cols())-fracture.col(i);
        double t = ((pointOnLine.cross(direction)).dot((ithFractureEdge.cross(direction)))-(ithFractureVertex.cross(direction)).dot((ithFractureEdge.cross(direction))))/pow((ithFractureEdge.cross(direction)).norm(), 2);
        if (t >= -tol && t <= 1.0+tol)
        {
            intersections.push_back(ithFractureVertex + t * ithFractureEdge);
        }
        if(intersections.size()==2){
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
bool fracturesIntersection(const MatrixXd& fracture1,const MatrixXd& fracture2, Vector3d& E1, Vector3d& E2, bool& tips1, bool& tips2)
{
    //Considero la retta di intersezione tra i due piani che contengono i poligoni
    //Calcolo la direction della retta come prodotto vettoriale tra le normali ai piani
    Vector3d n1 = normalP(fracture1);
    Vector3d n2 = normalP(fracture2);
    Vector3d direction = n1.cross(n2).normalized();
    //Se il prodotto vettoriale è nullo, i piani sono paralleli
    if(abs(direction.norm()) < tol){
        return false;
    }
    //Interseco il piano del poligono 1 con la retta che contiene il lato 1 del poligono 2 per trovare un punto della retta
    Vector3d edge = fracture2.col(1)-fracture2.col(0);
    Vector3d vertex = fracture2.col(0);
    //check che quel lato non sia parallelo al piano, e se lo è, prendo il lato successivo
    if(abs(edge.dot(n1))<tol)
    {
        edge = fracture2.col(2)-fracture2.col(1);
        vertex = fracture2.col(1);
    }
    //Soluzione dell'intersezione piano-retta
    double t = -(n1.dot(vertex-fracture1.col(0)))/(n1.dot(edge));
    Vector3d pointOnLine = vertex+t*edge;
    //Ora controlliamo se il poligono interseca la retta

    vector<Vector3d> intersections;
    intersections.reserve(4);
    vector<Vector3d> inter1;
    if(lineFractIntersect(fracture1,pointOnLine,direction,inter1)){
        intersections.push_back(applyThreshold(inter1[0]));
        intersections.push_back(applyThreshold(inter1[1]));
    }
    else{
        return false;
    }
    vector<Vector3d> inter2;
    if(lineFractIntersect(fracture2,pointOnLine,direction,inter2)){
        intersections.push_back(applyThreshold(inter2[0]));
        intersections.push_back(applyThreshold(inter2[1]));
    }
    else{
        return false;
    }
    //Identifico i 4 punti trovati sulla retta
    double a = ((intersections[0]-pointOnLine).dot(direction))/(direction.norm()*direction.norm());
    double b = ((intersections[1]-pointOnLine).dot(direction))/(direction.norm()*direction.norm());
    //a e b sono delle posizioni relative a pointOnLine, ma non so come sono ordinati perché questo dipende
    //dal senso di lettura dei lati dei poligoni e dalla posizione relativa alla retta
    if(abs(a-b)<tol){ //i due punti di intersezione sono sovrapposti
        return false;
    }
    if(a>b)
    {
        swap(a,b);
        swap(intersections[0],intersections[1]);
    }
    double c = ((intersections[2]-pointOnLine).dot(direction))/(direction.norm()*direction.norm());
    double d = ((intersections[3]-pointOnLine).dot(direction))/(direction.norm()*direction.norm());
    if(abs(b-c)<tol){
        return false;
    }

    if(c>d)
    {
        swap(c,d);
        swap(intersections[2],intersections[3]);
    }

    if(b-tol<=c || d-tol<=a)// i segmenti [a,b] e [c,d] non si intersecano
    {
        return false;
    }
    //c'è un'intersezione "propria", le tracce sono non-passanti per entrambe le fratture
    if(tol<c-a && tol<b-c && tol<d-b)
    {
        E1 = intersections[2];
        E2 = intersections[1];
        tips1 = true;
        tips2 = true;
        return true;
    }
    else if(tol<a-c && tol<d-a && tol<b-d)
    {
        E1 = intersections[0];
        E2 = intersections[3];
        tips1 = true;
        tips2 = true;
        return true;
    }
    //un segmento è contenuto nell'altro, possibilmente con uno o entrambi gli estremi che si sovrappongono
    else if(a-tol<=c && d<=b+tol)
    {
        E1 = intersections[2];
        E2 = intersections[3];
        tips2 = false;
        tips1 = true;
        if(abs(a-c)<tol && abs(d-b)<tol)
        {
            tips1 = false;
        }
        return true;
    }
    else if(c-tol<=a && b<=d+tol)
    {
        E1 = intersections[0];
        E2 = intersections[1];
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
bool near1(const MatrixXd& fracture1, const MatrixXd& fracture2)
{
    bool check = true;

    Vector3d max1 = fracture1.rowwise().maxCoeff();
    Vector3d min1 = fracture1.rowwise().minCoeff();
    Vector3d max2 = fracture2.rowwise().maxCoeff();
    Vector3d min2 = fracture2.rowwise().minCoeff();

    if (min1(0) > max2(0)+tol || min1(1) > max2(1)+tol || min1(2) > max2(2)+tol)
    {
        check = false;
    }
    if (min2(0) > max1(0)+tol|| min2(1) > max1(1)+tol || min2(2) > max1(2)+tol)
    {
        check = false;
    }

    return check;
}

//***************************************************************************************************************
// Calcolo punto medio del primo e del secondo poligono
// Calcolo la distanza tra essi e check che sia minore di una certa tolleranza
// 2 casi:
// se è <, allora c'è la possibilità che si intersechino
// se è >, allora non si intersecano
bool near2(const MatrixXd& fracture1, const MatrixXd& fracture2){
    bool check = true;
    VectorXd bar_fracture1 = fracture1.colwise().mean();
    VectorXd bar_fracture2 = fracture2.colwise().mean();

    double max_1 = 0.0;
    for (unsigned int j = 0; j < fracture1.col(0).size(); j++){
        double length1 = bar_fracture1.norm() - fracture1.col(j).norm();
        if (length1 > max_1){
            max_1 = length1;
        }
    }

    double max_2 = 0.0;
    for (unsigned int j = 0; j < fracture2.col(0).size(); j++){
        double length2 = bar_fracture2.norm() - fracture2.col(j).norm();
        if (length2 > max_2){
            max_2 = length2;
        }
    }

    double control_length = max_1 + max_2;
    double bar_distance = (bar_fracture1 - bar_fracture2).norm();

    if (bar_distance + tol < control_length){
        cout << "Intersection is possible" << endl;
    }
    else{
        check = false;
        cout << "Intersection is impossible" << endl;
    }
    return check;
}




//***************************************************************************************************************
//Serve per punto 2?
unsigned int IsInside(const MatrixXd& fracture, Vector3d normal, Vector3d point){
    unsigned int counter = 1; //ATTENZIONE essite bolean a 3?
    for (int i = 0; i < fracture.row(0).size(); ++i) { // cicla su tutti i lati di un poligono (colonne matrice)
        Vector3d edge = fracture.col((i + 1) % fracture.row(0).size()) - fracture.col(i); // vertice[i]-vertice[i+1]
        Vector3d toPoint = point - fracture.col(i); // punto-vertice
        double tripleProduct = normal.dot(edge.cross(toPoint)); //prodotto misto normale*(lato x toPoint)
        if (edge.cross(toPoint).norm()<tol){ // ATTENZIONE Tolleranza se no si rompe dove la definsico?????
            counter = 2; // sul bordo
        }
        if (tripleProduct < 0) {
            counter = 0; // fuori
            return counter;
        }
    }
    return counter; //interno
}

}

