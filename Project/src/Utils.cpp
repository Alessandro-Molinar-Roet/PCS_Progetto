#include<fstream>
#include<vector>
#include<iostream>
#include "tol.hpp"
#include "Utils.hpp"
#include"math.h"
#include<iomanip>
#include <algorithm>
#include <iomanip>

using namespace std;
namespace FractureNetwork {

void Define_tol(){
    //definizione delle tolleranze:
    double tol1D_user = 0.0;
    double tol2D_user = 0.0;

    cout << "Inserire tolleranza 1D: \n";
    while (!(cin >> tol1D_user)){
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
        cout << "Errore: formato inserito non supportato.\n";
        cout << "Inserire tolleranza 1D: \n";
    }

    cout << "Inserire tolleranza 2D: \n";
    while (!(cin >> tol2D_user)) {
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
        cout << "Errore: formato inserito non supportato.\n";
        cout << "Inserire tolleranza 2D: \n";
    }

    tol2D_user = max(standard, tol1D_user*tol1D_user);

    tol1 = max(standard, tol1D_user);
    tol2 = max(standard, tol2D_user);
}



bool ImportFractures(const string& filepath, vector<Fracture>& fractures){
    ifstream file;
    file.open(filepath);

    if(file.fail()){
        return false;
    }
    string line = "";
    getline(file,line); //header file
    getline(file,line); // fratture.num leggo, converto e salvo
    unsigned int num = stoi(line);
    fractures.resize(num);

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
            converter >> fractures[matrix_counter].ID >> temp >> num_vertices;
            fractures[matrix_counter].vertices.resize(3,num_vertices);
        }

        if((counter>=3) && (counter<=5)){ //linee di coordinate vertici (e1; e2; e3; e4; ...)
            for(unsigned int i = 0; i<num_vertices; i++){
                converter >> setprecision(12)>>fractures[matrix_counter].vertices(counter-3,i) >> temp;
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
void SortingFractureTraces(vector<Fracture>& fractures, vector<Trace>& traces){
    for (unsigned int i = 0; i<fractures.size(); i++){
        if(!fractures[i].passing.empty()){
            unsigned int num = size(fractures[i].passing);
            vector<pair<unsigned int,double>> temp(num);
            for(unsigned int j = 0; j<num; j++){
                unsigned int id = fractures[i].passing[j];
                temp[j] = make_pair(id,traces[id].length);
            }
            sort(temp.begin(), temp.end(), [](auto &left, auto &right) { return left.second > right.second;});
            //ATTENZIONE ALTERNATIVA #include<algorihm> e funzione transform
            for(unsigned int k = 0; k<num; k++){
                fractures[i].passing[k] = temp[k].first;
            }
        }

        if(!fractures[i].not_passing.empty()){
            unsigned int num = size(fractures[i].not_passing);
            vector<pair<unsigned int,double>> temp(num);
            for(unsigned int j = 0; j<num; j++){
                unsigned int id = fractures[i].not_passing[j];
                temp[j] = make_pair(id,traces[id].length);
            }
            sort(temp.begin(), temp.end(), [](auto &left, auto &right) { return left.second > right.second;});
            //ATTENZIONE ALTERNATIVA #include<algorihm> e funzione transform
            for(unsigned int k = 0; k<num; k++){
                fractures[i].not_passing[k] = temp[k].first;
            }
        }
    }
}

//***************************************************************************************************************
bool PrintTrace(const string& filepath, const vector<Trace>& traces){
    ofstream myfile;
    myfile.open(filepath);
    if(myfile.fail()){
        return false;
    }
    unsigned int num = traces.size();
    myfile << "# Number of Traces" << "\n";
    myfile << num << endl;
    for(unsigned int i = 0; i< num; i++){
        myfile << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;
        myfile << traces[i].ID << "; "
               << traces[i].first_generator << "; " << traces[i].second_generator << "; "
               << fixed<<setprecision(12) << traces[i].vertices(0,0) << "; " << traces[i].vertices(1,0) << "; " << traces[i].vertices(2,0) << "; "
               << traces[i].vertices(0,1) << "; " << traces[i].vertices(1,1) << "; " << traces[i].vertices(2,1) << endl;
    }
    myfile.close();
    return true;
}


//***************************************************************************************************************
bool PrintFractureTraces(const string& filepath, const  vector<Fracture>& fractures, const vector<Trace>& traces){
    ofstream myfile;
    myfile.open(filepath);

    if(myfile.fail()){        
        return false;
    }

    for (unsigned int i = 0; i<fractures.size(); i++){

        vector<unsigned int> passing = fractures[i].passing;
        vector<unsigned int> not_passing = fractures[i].not_passing;
        if(!passing.empty() || !not_passing.empty()){
            myfile << "# FractureId; NumTraces" << "\n";
            myfile << i << "; " << passing.size() + not_passing.size() << "\n";

            for(unsigned int i = 0; i< passing.size(); i++){
                unsigned int id = passing[i];
                myfile << "# TraceId; Tips; Length" << "\n";
                myfile << id << "; false; " << traces[id].length << "\n";
            }

            for(unsigned int i = 0; i< not_passing.size(); i++){
                unsigned int id = not_passing[i];
                myfile << "# TraceId; Tips; Length" << "\n";
                myfile << id << "; true; " << traces[id].length << "\n";
            }
            myfile << "\n\n";
        }
    }
    return true;
}

}


