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

}


