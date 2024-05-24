#include "Geometry.hpp"
#include "math.h"
#include <iostream>

namespace FractureNetwork {

//calcola tutte le traccie e se sono passatnti o non passatni e le salva DOVE E COME??????
void CalculateFracture(const Fractures& fratture, Traces& tracce){
    unsigned int counter = 0;
    for(unsigned int i = 0; i<fratture.num; i++){ //ciclo su tutte le matrici (piani)
        MatrixXd first_polygon = fratture.f_Vertices[i];
        for(unsigned int j = i+1; j<fratture.num; j++ ){ // controllo tutte le coppie da >i (altre gia cocntrollate precendemente)
            MatrixXd second_polygon = fratture.f_Vertices[j];

            if(near(first_polygon, second_polygon)){
                //ATTENZIONE tutto vettori ma non consoco dimensione e ho bisogno di accedervi spesso LISTA o RESERVE?
                Matrix<double, 3, 2> vertices = {};
                bool find = false;
                bool complete = false;
                bool tips1 = true;
                bool tips2 = true;
                tips1 = FindIntersection(first_polygon, second_polygon,vertices,find,complete,tips2);
                if(!complete){
                    tips2 = FindIntersection(second_polygon,first_polygon, vertices,find,complete,tips1); //ATTENZIONE funziona? come miglioro leggibilità codice?
                }
                //ATTENZIONE:
                //DEVO salvare tips?
                if(complete){
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


 bool FindIntersection(const MatrixXd& first_polygon, const MatrixXd& second_polygon, Matrix<double, 3, 2>& vertices, bool& find, bool& complete,bool& tips2){
    bool tips = true;
    unsigned int edges_counter = 0;
    unsigned int edges_on = 0;
    //equazione del piano di second_polygon ax + by + cz + d = 0
    //normal = vettore normale passante per polygon (a,b,c)
    Vector3d vector1 = second_polygon.col(0)-second_polygon.col(1);
    Vector3d vector2 = second_polygon.col(0)-second_polygon.col(2);
    Vector3d normal = vector1.cross(vector2).transpose();
    double d = -second_polygon.col(0).dot(normal);

    //oggetti first_polygon necessari per ciclo ATTENZIONE posso dichiaralri dentro o fuori
    int num = first_polygon.row(0).size();
    Vector3d point1;
    Vector3d point2;
    Vector3d distance;
    Vector3d pl_intersection;
    //ciclo sui lati di first-polygon, estendo a retta formula parametrica e faccio intersezione retta piano
    // se non in [0,1] vuol dire che sono andato o oltre o indietro rsipetto al mio segmento no intersezione
    for(int i = 0; i<num; i++){
        point1 = first_polygon.col(i);
        point2 = first_polygon.col((i + 1) % num);

        distance = point2-point1; // point2-point1 no il contrario si rompe tutto
        double t = -1;
        if(normal.dot(distance)!= 0){ //evito di divedere per 0
            t = -(normal.dot(point1) + d)/ (normal.dot(distance)); //metto (x,y,z)=f(t) in equaizone del piano e trovo t
            pl_intersection = point1 + t * distance; // espressione parametrica della retta:  point1+t*disctance
        }
        else{ // it lies on entirly
            if((point1.dot(normal) + d == 0) ||(point2.dot(normal) + d == 0) ){
                cout << "line lais on" << endl; //MANCA ATTENZIONE (controllare se entrambi i punti risolvono equazione del piano
            }
        }

        if(t>= 0.0 && t<= 1.00){
            // counter = 0 sta fuori / counter = 1 sta all'inerno / counter = 2 sta sul bordo
            unsigned int counter = IsInside(second_polygon, normal, pl_intersection);
            if(counter > 0 ){
                edges_counter ++;
                if(!find){
                    vertices.col(0) = pl_intersection;
                    find = true;
                }
                else{
                    vertices.col(1) = pl_intersection;
                    complete = true;
                }
                if(counter == 2){ //se sul bordo dell'altro poligono
                    edges_on ++;
                }
            }
        }
    }
    if(edges_counter == 2){
        tips = false; // è passante per first polygon salvare in struddura giusta tips
    }
    if(edges_on == 2){
        cout << "shortcut" << endl; //ATTENZIONE funziona? posso migliorare leggibilita?
        tips2 = false;
    }
    return tips;
}

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

bool near(const MatrixXd& first_polygon, const MatrixXd& second_polygon){
    return true;
}


}
