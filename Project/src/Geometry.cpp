#include "Geometry.hpp"
#include "math.h"
#include <iostream>

namespace FractureNetwork {

//calcola tutte le traccie e se sono passatnti o non passatni e le salva DOVE E COME??????
void CalculateTraces(const Fractures& fratture, Traces& tracce){
    unsigned int counter = 0;
    for(unsigned int i = 0; i<fratture.num; i++){ //ciclo su tutte le matrici (piani)
        MatrixXd first_polygon = fratture.f_Vertices[i];
        for(unsigned int j = i+1; j<fratture.num; j++ ){ // controllo tutte le coppie da >i (altre gia cocntrollate precendemente)
            MatrixXd second_polygon = fratture.f_Vertices[j];

            if(near(first_polygon, second_polygon)){
                Vector3d E1 = {};
                Vector3d E2 = {};
                bool tips1 = true;
                bool tips2 = true;
                bool find = TracciaTraPoligoni(first_polygon, second_polygon, E1, E2,tips1, tips2);

                if(find){
                    Matrix<double,3,2> vertices = {};
                    vertices.col(0) = E1;
                    vertices.col(1) = E2;
                    tracce.t_Vertices.push_back(vertices);
                    double length = sqrt( pow((vertices(0,0)-vertices(0,1)),2) +  pow((vertices(1,0)-vertices(1,1)),2) +  pow((vertices(2,0)-vertices(2,1)),2)); //ATTENZIONE berrone no
                    tracce.t_length.push_back(length);
                    tracce.t_ID.push_back(counter); //ID della traccia = numero in cui è stata trovata

                    tracce.t_generator.push_back(fratture.f_ID[i]); // ATTENZIONE come salvo gli ID dei piani che generano uno frattura? posizione = 2i e 2i-1?
                    tracce.t_generator.push_back(fratture.f_ID[j]);

                    tracce.Dfn.insert({ i, {{},{}} });
                    if(!tips1){
                        tracce.Dfn[i][0].push_back(counter);
                    }
                    else{
                        tracce.Dfn[i][1].push_back(counter);
                    }
                    tracce.Dfn.insert({j, {{},{}} });
                    if(!tips2){
                        tracce.Dfn[j][0].push_back(counter);
                    }
                    else{
                        tracce.Dfn[j][1].push_back(counter);
                    }

                    counter ++;
                    tracce.num = counter;
                    }
                }
            }
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
bool TracciaTraPoligoni(const MatrixXd& frattura1,const MatrixXd& frattura2, Vector3d& E1, Vector3d& E2, bool& tips1, bool& tips2)
{
    //Considero la retta di intersezione tra i due piani che contengono i poligoni
    //Calcolo la direzione della retta come prodotto vettoriale tra le normali ai piani
    Vector3d n1 = normaleP(frattura1);
    Vector3d n2 = normaleP(frattura2);
    Vector3d direzione = n1.cross(n2).normalized();
    //Se il prodotto vettoriale è nullo, i piani parallele
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
    for(unsigned int i=0;i<frattura1.cols();i++){
        //Nel momento in cui trovo due intersezioni, non mi serve controllare i lati rimanenti
        if(intersezioni.size()==2)
        {
            break;
        }
        //intersezione lato-retta
        Vector3d puntoipoly1 = frattura1.col(i);
        Vector3d latoipoly1 =frattura1.col((i+1)%frattura1.cols())-frattura1.col(i);
        double t = ((puntoRetta.cross(direzione)).dot((latoipoly1.cross(direzione)))-(puntoipoly1.cross(direzione)).dot((latoipoly1.cross(direzione))))/pow((latoipoly1.cross(direzione)).norm(), 2);
        if (t >= 0.0 && t <= 1.0)
        {
            intersezioni.push_back(puntoipoly1 + t * latoipoly1);
        }
    }
    //Se il primo poligono non ha intersezioni con la retta, non serve controllare il secondo poligono
    if (intersezioni.size() < 2)
    {
        return false;
    }

    for(unsigned int i=0;i<frattura2.cols();i++){
        if(intersezioni.size()==4)
        {
            break;
        }

        Vector3d puntoipoly2 = frattura2.col(i);
        Vector3d latoipoly2 =frattura2.col((i+1)%frattura2.cols())-frattura2.col(i);
        double t = ((puntoRetta.cross(direzione)).dot((latoipoly2.cross(direzione)))-(puntoipoly2.cross(direzione)).dot((latoipoly2.cross(direzione))))/pow((latoipoly2.cross(direzione)).norm(), 2);
        if (t >= 0.0 && t <= 1.0)
        {
            intersezioni.push_back(puntoipoly2 + t * latoipoly2);
        }

    }
    if (intersezioni.size() < 4)
    {
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

    cout << "problem" << endl;
    return false;
}


//***************************************************************************************************************
bool near(const MatrixXd& first_polygon, const MatrixXd& second_polygon){
    return true;
}


















//***************************************************************************************************************
//Serve per punto 2?
unsigned int IsInside(const MatrixXd& polygon, Vector3d normal, Vector3d point){
    unsigned int counter = 1; //ATTENZIONE essite bolean a 3?
    for (int i = 0; i < polygon.row(0).size(); ++i) { // cicla su tutti i lati di un poligono (colonne matrice)
        Vector3d lato = polygon.col((i + 1) % polygon.row(0).size()) - polygon.col(i); // vertice[i]-vertice[i+1]
        Vector3d toPoint = point - polygon.col(i); // punto-vertice
        double prodotto_misto = normal.dot(lato.cross(toPoint)); //prodotto misto normale*(lato x toPoint)
        if (lato.cross(toPoint).isZero(0.000001)){ // ATTENZIONE Tolleranza se no si rompe dove la definsico?????
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
