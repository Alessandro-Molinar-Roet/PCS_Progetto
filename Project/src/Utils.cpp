#include<fstream>
#include<vector>
#include<iostream>
#include"Utils.hpp"
#include"math.h"
#include "Eigen/Dense"

using namespace Eigen;
using namespace std;

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
    cout << counter1 << endl;
    cout << counter2 << endl;
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


bool near1(const MatrixXd& first_polygon, const MatrixXd& second_polygon)
{

    bool controllo = true;
    // Bounding box:

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


    // Vector3d max1;
    // Vector3d max2;
    // Vector3d min1;
    // Vector3d min2;


    // for (unsigned int i = 0; i < first_polygon.row(0).size(); i++)
    // {
    //     double max = first_polygon(i)(0);
    //     for (unsigned int j = 0; j < first_polygon.col(0).size(); j++)
    //     {
    //         if (first_polygon(i)(j) >= max)
    //         {
    //             max = first_polygon(i)(j);
    //             max1(i) = max;
    //         }
    //     }
    // }

    // for (unsigned int i = 0; i < second_polygon.row().size(); i++)
    // {
    //     double max = second_polygon(i)(0);
    //     for (unsigned int j = 0; j < second_polygon.col(0).size(); j++)
    //     {
    //         if (second_polygon(i)(j) >= max)
    //         {
    //             max = second_polygon(i)(j);
    //             max2(i) = max;
    //         }
    //     }
    // }

    // for (unsigned int i = 0; i < first_polygon.row().size(); i++)
    // {
    //     double min = first_polygon(i)(0);
    //     for (unsigned int j = 0; j < first_polygon.col(0).size(); j++)
    //     {
    //         if (first_polygon(i)(j) <= min)
    //         {
    //             min = first_polygon(i)(j);
    //             min1(i) = min;
    //         }
    //     }
    // }

    // for (unsigned int i = 0; i < second_polygon.row().size(); i++)
    // {
    //     double min = second_polygon(i)(0);
    //     for (unsigned int j = 0; j < second_polygon.col(0).size(); j++)
    //     {
    //         if (second_polygon(i)(j) >= min)
    //         {
    //             min = second_polygon(i)(j);
    //             min2(i) = min;
    //         }
    //     }
    // }


}
bool near2(const MatrixXd& first_polygon, const MatrixXd& second_polygon)

{
    // Calcolo punto medio del primo e del secondo poligono
    // Calcolo la distanza tra essi e controllo che sia minore di una certa tolleranza
    // 2 casi:
    // se è <, allora c'è la possibilità che si intersechino
    // se è >, allora non si intersecano

    // Vecchio metodo:

    bool controllo = true;

    VectorXd bar_polygon1 = first_polygon.colwise().mean();
    VectorXd bar_polygon2 = second_polygon.colwise().mean();

    double max_1 = 0.0;
    for (unsigned int j = 0; j < first_polygon.col(0).size(); j++)
    {
        double length1 = bar_polygon1.norm() - first_polygon.col(j).norm();
        if (length1 > max_1)
        {
            max_1 = length1;
        }
    }

    double max_2 = 0.0;
    for (unsigned int j = 0; j < second_polygon.col(0).size(); j++)
    {
        double length2 = bar_polygon2.norm() - second_polygon.col(j).norm();
        if (length2 > max_2)
        {
            max_2 = length2;
        }
    }
    double control_length = max_1 + max_2;

    double bar_distance = (bar_polygon1 - bar_polygon2).norm();

    if (bar_distance < control_length)
    {
        cout << "Intersection is possible" << endl;
    }
    else
    {
        controllo = false;
        cout << "Intersection is impossible" << endl;
    }
    return controllo;
}


