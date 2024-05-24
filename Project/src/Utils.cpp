#include<fstream>
#include<vector>
#include<iostream>
#include"Utils.hpp"
#include"math.h"

using namespace std;
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
        istringstream converter(line); //ATTENZIONE efficente crearlo anche qundo non lo uso? SOLUZIONE If line[0] = #
        char temp = ' ';

        if(counter == 1){ //linea ID; num_vertici
            converter >> fratture.f_ID[matrix_counter] >> temp >> num_vertices;
            fratture.f_Vertices[matrix_counter].resize(3,num_vertices);
        }

        if((counter>=3) && (counter<=5)){ //linee di coordinate vertici (e1; e2; e3; e4; ...)
            for(unsigned int i = 0; i<num_vertices; i++){
                converter >> fratture.f_Vertices[matrix_counter](counter-3,i) >> temp;
            }
        }

        if(counter == 6){ //linea intestazione (sommo subito dopo quindi counter = 1)
            matrix_counter ++;
            counter = 0;
        }
        counter ++;
    }
    return true;
}

//********************************************************************************
void SortingFractureTraces(Traces& tracce){
    // ATTENZIONE casino
    for (auto const& [key, val] : tracce.Dfn)
    {
        vector passing  = val[0];
        unsigned int num_p = size(passing);
        vector<pair<unsigned int,double>> vp(num_p);
        for(unsigned int i = 0; i<num_p; i++){
            unsigned int id = passing[i];
            vp[i] = make_pair(id,tracce.t_length[id]);
        }
        sort(vp.begin(), vp.end(), [](auto &left, auto &right) { return left.second < right.second;});
        for(unsigned int i = 0; i<num_p; i++){
            tracce.Dfn[key][0][i] = vp[i].first;
        }

        vector not_passing  = val[1];
        unsigned int num_n = size(not_passing);
        vector<pair<unsigned int,double>> vn(num_n);
        for(unsigned int i = 0; i<num_n; i++){
            unsigned int id = not_passing[i];
            vn[i] = make_pair(id,tracce.t_length[id]);
        }
        sort(vn.begin(), vn.end(), [](auto &left, auto &right) { return left.second < right.second;});
        for(unsigned int i = 0; i<num_n; i++){
            tracce.Dfn[key][1][i] = vn[i].first;
        }
    }
}

//***************************************************************************************************************
bool PrintTrace(const string& filepath, const Traces& tracce){
    ofstream myfile;
    myfile.open(filepath);
    if(myfile.fail()){
        return false;
    }
    myfile << "# Number of Traces" << "\n";
    myfile << tracce.num << endl;
    for(unsigned int i = 0; i<tracce.t_ID.size(); i++){
        myfile << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;
        myfile << tracce.t_ID[i] << "; "
               << tracce.t_generator[2*i +1] << "; " << tracce.t_generator[2*i + 2] << "; "
               << tracce.t_Vertices[i](0,0) << "; " << tracce.t_Vertices[i](1,0) << "; " << tracce.t_Vertices[i](2,0) << "; "
               << tracce.t_Vertices[i](0,1) << "; " << tracce.t_Vertices[i](1,1) << "; " << tracce.t_Vertices[i](2,1) << endl;
    }
    myfile.close();
    return true;
}


//***************************************************************************************************************
bool PrintFractureTraces(const string& filepath, const Traces& tracce){
    ofstream myfile;
    myfile.open(filepath);

    if(myfile.fail()){
        return false;
    }

    for (auto const& [key, val] : tracce.Dfn){ //ATTENZIONE va bene auto per leggibilità o meglio tipo vero?
        myfile << "# FractureId; NumTraces" << "\n";
        myfile << key << "; " << size(val[0]) + size(val[1])  << "\n";

        for(unsigned int i = 0; i<size(val[0]); i++){
            unsigned int id = val[0][i];
            myfile << "# TraceId; Tips; Length" << "\n";
            myfile << id << "; false; " << tracce.t_length[id] << "\n";
        }

        for(unsigned int i = 0; i<size(val[1]); i++){
            unsigned int id = val[1][i];
            myfile << "# TraceId; Tips; Length" << "\n";
            myfile << id << "; true; " << tracce.t_length[id] << "\n";
        }
        myfile << "\n\n";
    }

    return true;
}
//**********************************************************************************
//Per ogni poligono, trovo la normale del piano su cui giace
Vector3d normaleP(const MatrixXd& frattura)
{
    Vector3d v1 = frattura.col(1)-frattura.col(0);
    Vector3d v2 = frattura.col(2)-frattura.col(1);
    return v1.cross(v2);
}

bool TracciaTraPoligoni(const MatrixXd& frattura1,const MatrixXd& frattura2, Vector3d& E1, Vector3d& E2, bool& tips1, bool& tips2)
{
    //Calcolo la direzione della retta
    Vector3d n1 = normaleP(frattura1);
    Vector3d n2 = normaleP(frattura2);
    Vector3d direzione = n1.cross(n2).normalized();
    //Se il prodotto vettoriale è nullo, le normali sono parallele e i piani non si intersecano
    if(abs(direzione.norm()) < 0.00001){ //abs<tol
        return false;
    }
    //Interseco il piano del poligono 1 con la retta che contiene il lato 1 del poligono 2 per trovare un punto della retta
    Vector3d lato = frattura2.col(1)-frattura2.col(0);
    Vector3d punto = frattura2.col(0);
    //Controllo che quel lato non sia parallelo al piano, e se lo è, prendo il lato successivo
    if(lato.dot(n1)==0){
        lato = frattura2.col(2)-frattura2.col(1);
        punto = frattura2.col(1);
    }
    //Soluzione dell'intersezione piano retta
    double t = -(n1.dot(punto-frattura1.col(0)))/(n1.dot(lato));
    Vector3d puntoRetta = punto+t*lato;
    //Ora controlliamo se il poligono interseca la retta


    vector<Vector3d> intersezioni;
    for(unsigned int i=0;i<frattura1.cols();i++){
        if(intersezioni.size()==2)
        {
            break;
        }
        Vector3d puntoipoly1 = frattura1.col(i);
        Vector3d latoipoly1 =frattura1.col((i+1)%frattura1.cols())-frattura1.col(i);
        double t = ((puntoRetta.cross(direzione)).dot((latoipoly1.cross(direzione)))-(puntoipoly1.cross(direzione)).dot((latoipoly1.cross(direzione))))/pow((latoipoly1.cross(direzione)).norm(), 2);
                if (t >= 0.0 && t <= 1.0)
        {
                intersezioni.push_back(puntoipoly1 + t * latoipoly1);
                }
    }
    if (intersezioni.size() < 2) {
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

    double a = ((intersezioni[0]-puntoRetta).dot(direzione))/(direzione.norm()*direzione.norm());
    double b = ((intersezioni[1]-puntoRetta).dot(direzione))/(direzione.norm()*direzione.norm());
    //a e b sono delle posizioni relative a puntoRetta, ma non so come sono ordinati perché questo dipende
    //dal senso di lettura dei lati dei poligoni e dalla posizione relativa alla retta
    if(a>b)
    {
        double atemp = a;
        a=b;
        b=atemp;
    }
    double c = ((intersezioni[2]-puntoRetta).dot(direzione))/(direzione.norm()*direzione.norm());
    double d = ((intersezioni[3]-puntoRetta).dot(direzione))/(direzione.norm()*direzione.norm());

    if(c>d)
    {
        double ctemp =c;
        c=d;
        d=ctemp;
    }


    if(b<=c || d<=a)// i segmenti [a,b] e [c,d] non si intersecano
    {
        return false;
    }
    //c'è un'intersezione "propria"
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
        if(a==c && d == b){
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
        if(a==c && d == b){
            tips2 = false;
        }
        return true;

    }




}



