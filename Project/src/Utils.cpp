#include<fstream>
#include"Utils.hpp"
#include"math.h"

namespace FractureNetwork {

bool ImportFractures(const string& filepath, Fractures& fratture){
    ifstream file;
    file.open(filepath);

    if(file.fail()){
        return false;
    }

    string line = "";
    getline(file,line); //header file
    getline(file,line); // fratture.num leggo, converto e salvo
    fratture.num = stoi(line);

    // alloco dimensione vettore di matrici
    fratture.f_Vertices.resize(fratture.num);
    fratture.f_ID.resize(fratture.num);

    //ciclo su linee file
    //matrix_counter conta quante matrici ho gia letto
    //counter conta per ogni matrice a che linea di info sono arrivato (ogni matrice occupa 6 righe [0,5])
    unsigned int counter = 0;
    unsigned int matrix_counter = 0;
    unsigned int num_vertices = 0;
    while (getline(file, line)){
        istringstream converter(line); //ATTENZIONE efficente crearlo anche qundo non lo uso?
        char temp = ' ';
        if(counter == 1){ //linea ID; num_vertici
            converter >> fratture.f_ID[matrix_counter] >> temp >> num_vertices;
            fratture.f_Vertices[matrix_counter].resize(3,num_vertices);
        }
        if((counter>=3) && (counter<=5)){ //linee di coordinate vertici (e1; e1; e3; e4; ...)
            for(unsigned int i = 0; i<num_vertices; i++){
                converter >> fratture.f_Vertices[matrix_counter](counter-3,i) >> temp;
            }
        }
        if(counter == 6){ //linea intestazione (sommando subito dopo salto header)
            matrix_counter ++;
            counter = 0;
        }
        counter ++;
    }
    return true;
}

//********************************************************************************

//calcola tutte le traccie e se sono passatnti o non passatni e le salva DOVE E COME??????
void CalculateFracture(const Fractures& fratture, Traces& tracce){
    unsigned int counter = 0;
    for(unsigned int i = 0; i<fratture.num; i++){ //ciclo su tutte le matrici (piani)
        MatrixXd first_polygon = fratture.f_Vertices[i];
        for(unsigned int j = i+1; j<fratture.num; j++ ){ // controllo tutte le coppie da >i (altre gia cocntrollate precendemente)
            MatrixXd second_polygon = fratture.f_Vertices[j];

            if(near(first_polygon, second_polygon)){
                Matrix<double, 3, 2> vertices = {};
                bool find = false;
                bool complete = false;
                FindIntersection(first_polygon, second_polygon,vertices,find,complete);
                FindIntersection(second_polygon,first_polygon, vertices,find,complete);
                //ATTENZIONE:
                //ripeto due volte funzione posso evitarlo? servono tante flags
                //DEVO salvare passanti non passatni: più facile falro all'intenro di FindIntersection per passanti ma per non passanti no (bool su findIntersecition)
                if(complete){
                    tracce.t_Vertices.push_back(vertices);
                    double length = sqrt( pow((vertices(0,0)-vertices(0,1)),2) +  pow((vertices(1,0)-vertices(1,1)),2) +  pow((vertices(2,0)-vertices(2,1)),2)); //ATTENZIONE berrone no
                    tracce.t_length.push_back(length);
                    tracce.t_ID.push_back(counter); //ID della traccia = numero in cui è stata trovata

                    tracce.t_generator.push_back(fratture.f_ID[i]); // ATTENZIONE come salvo gli ID dei piani che generano uno frattura? posizione = 2i e 2i-1?
                    tracce.t_generator.push_back(fratture.f_ID[j]);

                    counter ++;
                    tracce.num = counter;
                }
            }
        }
    }
}

void FindIntersection(const MatrixXd& first_polygon, const MatrixXd& second_polygon, Matrix<double, 3, 2>& vertices, bool& find, bool& complete){
    unsigned int edges_counter = 0;
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
            if(counter == 2){
                //ATTENZIONE serve?? un punto è sul vertice del poligono se lo sono entrambi potrei non fare il viceversa ma ne vale la pena??
            }
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
            }
        }
    }
    if(edges_counter == 2){
        // è passante per first polygon salvare in struddura giusta tips
    }
}

unsigned int IsInside(const MatrixXd& polygon, Vector3d normal, Vector3d point){
    unsigned int counter = 1; //ATTENZIONE essite bolean a 3?
    for (int i = 0; i < polygon.row(0).size(); ++i) { // cicla su tutti i lati di un poligono (colonne matrice)
        Vector3d lato = polygon.col((i + 1) % polygon.row(0).size()) - polygon.col(i); // vertice[i]-vertice[i+1]
        Vector3d toPoint = point - polygon.col(i); // punto-vertice
        double prodotto_misto = normal.dot(lato.cross(toPoint)); //prodotto misto normale*(lato x toPoint)
        if (lato.cross(toPoint).isZero(0.000001)){ // ATTENZIONE Tolleranza se no si rompe dove la definsico?????
            counter = 2; // sul bordo
            return counter;
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
