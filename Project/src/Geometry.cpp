#include "Geometry.hpp"
#include "math.h"
#include <iostream>
#include<iomanip>
#include<tol.hpp>

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
bool TracciaTraPoligoni(const MatrixXd& frattura1,const MatrixXd& frattura2, Vector3d& E1, Vector3d& E2, bool& tips1, bool& tips2)
{
    //Considero la retta di intersezione tra i due piani che contengono i poligoni
    //Calcolo la direzione della retta come prodotto vettoriale tra le normali ai piani
    Vector3d n1 = normaleP(frattura1);
    Vector3d n2 = normaleP(frattura2);
    Vector3d direzione = n1.cross(n2).normalized();
    //Se il prodotto vettoriale è nullo, i piani sono paralleli
    if(abs(direzione.norm()) < 0.00001){ //abs<tol
        return false;
    }
    //Interseco il piano del poligono 1 con la retta che contiene il lato 1 del poligono 2 per trovare un punto della retta
    Vector3d lato = frattura2.col(1)-frattura2.col(0);
    Vector3d punto = frattura2.col(0);
    //Controllo che quel lato non sia parallelo al piano, e se lo è, prendo il lato successivo
    if(lato.dot(n1)==0)
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
        intersezioni.push_back(inter1[0]);
        intersezioni.push_back(inter1[1]);
    }
    else{
        return false;
    }
    vector<Vector3d> inter2;
    if(intersRettaPoly(frattura2,puntoRetta,direzione,inter2)){
        intersezioni.push_back(inter2[0]);
        intersezioni.push_back(inter2[1]);
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
        if(a==c && d == b)
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
        if(a==c && d == b)
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
bool IsInside(const MatrixXd& polygon, Vector3d point){
    Vector3d point1 = polygon.col(0);
    Vector3d point2 = polygon.col(1);
    Vector3d N = (point-point1).cross(point-point2);

    for (unsigned int i = 0; i < polygon.cols(); i++) { // cicla su tutti i lati di un poligono (colonne matrice)
        Vector3d pi = polygon.col(i);
        Vector3d ps = polygon.col((i + 1) % polygon.cols());
        double prodotto_misto = (point-pi).cross(point-ps).dot(N); //prodotto misto normale*(lato x toPoint) //ATTENZIONE megliodot product o sequenza if?
        if(prodotto_misto < 0){
            return false; //esterno
        }
    }
    return true; //interno
}

//***************************************************************************************************************
vector<unsigned int> intersect(const MatrixXd& frattura, const Vector3d& puntoRetta,const Vector3d& direzione,vector<Vector3d>& intersezioni){
    vector<unsigned int> lato;
    unsigned int counter = 1;

    for(unsigned int i= 0; i<frattura.cols();i++){
        Vector3d puntoipoly = frattura.col(i);
        Vector3d latoipoly =frattura.col((i+1) % frattura.cols())-frattura.col(i);

        if(!(latoipoly.cross(direzione)).isZero(0)){
            double t = ((puntoRetta.cross(direzione)).dot((latoipoly.cross(direzione)))-(puntoipoly.cross(direzione)).dot((latoipoly.cross(direzione))))/(latoipoly.cross(direzione).norm()*latoipoly.cross(direzione).norm());
            if ( (t > (0)) && (t < (1.0)) ){
                intersezioni.push_back(puntoipoly + t * latoipoly);
                lato.push_back(counter);
            }
        }
        else{
            /*
            cout << "paralelli" << endl;
            cout << "controllo libro" << endl;
            */
        }
        if(intersezioni.size() == 2 ){
            return lato;
        }
        counter++; //spostare counter
    }

    cout << "error: " << intersezioni.size() << "\n" << frattura.transpose() << "\n" << direzione <<"\n" << puntoRetta <<  endl;
    return lato;
}


void split(const MatrixXd& polygon, const vector<unsigned int>& all, const vector<Trace>& tracce, unsigned int counter, list<MatrixXd>& cutted ){
    //inizializzazione variabili
    unsigned int current;
    Vector3d point1 = {};
    Vector3d point2 = {};
    Vector3d direzione = {};
    bool full1 = false;
    bool full2 = false;

    if(counter < all.size()){
        current = all[counter]; // ID traccia tagliante
        point1 = tracce[current].vertices.col(0); // estremo 1
        point2 = tracce[current].vertices.col(1); // estremo 2
        direzione = point2-point1; // direzione dell traccia

        full1 = IsInside(polygon,point1); // true se la frattura contiene l'estremo uno della traccia
        full2 = IsInside(polygon,point2); // true se la frattura contiene l'estremo due della traccia
    }

    if( !full1 && !full2){
        cutted.push_back(polygon);
        // extractinfo(polygon, mesh);
    }
    else{
        counter++;
        unsigned int dim = polygon.cols();
        MatrixXd temp_sx(3,dim+2);
        MatrixXd temp_dx(3,dim+2);
        vector<Vector3d> intersezioni;
        vector<unsigned int> lati = intersect(polygon, point1, direzione, intersezioni);

        bool change = false;
        unsigned int counter1 = 0;
        unsigned int counter2 = 0;
        if(lati[0] != lati[1]){
            for(unsigned int i = 0; i< dim ; i++ ){
                if(!change){
                    temp_sx.col(counter1) = polygon.col(i);
                    counter1 ++;
                }
                if(change){
                    temp_dx.col(counter2) = polygon.col(i);
                    counter2 ++;
                }
                if( (lati[0]-1) == i){
                    if(intersezioni[0]!=polygon.col(i)){
                        temp_sx.col(counter1) = intersezioni[0];
                        counter1++;
                    }
                    change = true;
                }
                if((lati[1]-1) == i){
                    if(intersezioni[1]!=polygon.col(i)){
                        temp_sx.col(counter1) = intersezioni[1];
                        counter1 ++;

                        temp_dx.col(counter2) = intersezioni[1];
                        counter2 ++;
                        temp_dx.col(counter2) = intersezioni[0];
                        counter2 ++;
                    }
                    change = false;
                }
            }

            // estraggo sotto matrice riempita
            MatrixXd sx = temp_sx.leftCols(counter1);
            MatrixXd dx = temp_dx.leftCols(counter2);

            // chiamate ricorsive
            split(sx, all, tracce, counter, cutted);
            split(dx, all, tracce, counter, cutted);
        }
        else{ // caso a libro non effettuo tagli passo a frattura successiva
            // chiamata ricorsive
            cout << "intersezione a libro" << endl;
            split(polygon,all,tracce,counter,cutted);
        }
    }
}

void cutting(vector<Fracture>& fratture, vector<Trace>& tracce){
    list<MatrixXd> cutted = {};  // lista con le fratture tagliate (ogni matrice ha i vertici di un taglio)
    Mesh mesh = {};
    // Ciclo su tutte le fratture
    for(unsigned int i = 0; i<fratture.size(); i++){
        MatrixXd polygon = fratture[i].vertices; // vertici frattura
        // Se ce ne sono faccio il cut per le passanti
        vector<unsigned int> passing = fratture[i].passing; // vettore id traccie passanti
        vector<unsigned int> not_passing = fratture[i].not_passing;

        vector<unsigned int> all;
        all.reserve(passing.size() + not_passing.size());
        all.insert(all.end(), passing.begin(), passing.end());
        all.insert(all.end(), not_passing.begin(), not_passing.end());

        if(!all.empty()){
            unsigned int counter = 0;
            split(polygon, all , tracce, counter, cutted); // funzione ricorsiva taglio finchè posso
        }
    }

    for (auto const& i : cutted){
        cout << i << endl << endl;
    }

    cout << cutted.size() << endl;
    cout << tracce.size() << endl;
}





}

