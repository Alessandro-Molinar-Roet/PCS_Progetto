#pragma once

#include "../src/FracturesNetworkLibrary.hpp"
#include "UCDUtilities.hpp"
#include <iostream>

using namespace std;
using namespace FractureNetwork;
using namespace Gedim;

void Export_paraview( vector<Fracture> &f,  vector<Trace>&t){

    // numero totale di punti
    int N = 0;
    for (const auto&  polygon: f){
        MatrixXd vertices = polygon.vertices;
        N+= (vertices.cols()-2)*3;
    }

    // metto tutti i punti in una grande matrice
    MatrixXd points_p(3,N);
    unsigned int col = 0;
    for (const auto& polygon: f) {
        MatrixXd vertices = polygon.vertices;
        Vector3d fixed = vertices.col(0);
        for (unsigned int i = 1; i < polygon.vertices.cols()-1; i++){
            points_p.col(col) = fixed;
            col++;
            points_p.col(col) = polygon.vertices.col(i);;
            col++;
            points_p.col(col) = polygon.vertices.col(i+1);
            col++;
        }
    }

    // creo vettore con id raggruppati per poligono
    vector<vector<unsigned int>> polygon_vertices;
    unsigned int start_index = 0;
    for (const auto&  polygon: f){
        MatrixXd vertices = polygon.vertices;

        //triangolazione del poligono
        vector<Matrix3d> triangolation;
        Vector3d fixed = vertices.col(0);
        for (unsigned int i = 1; i < polygon.vertices.cols()-1; i++){

            Matrix3d temp(3,3);
            temp.col(0) = fixed;
            temp.col(1) = polygon.vertices.col(i);
            temp.col(2) = polygon.vertices.col(i+1);
            triangolation.push_back(temp);
        }

        //stesso id per stesso
        for(unsigned int i = 0; i < triangolation.size(); i++){
            vector<unsigned int> polygon_ids;
            Matrix3d temp = triangolation[i];
            for(unsigned j = 0; j< 3; j++){
                polygon_ids.push_back(start_index);
                start_index ++;
            }
            polygon_vertices.push_back(polygon_ids);
        }
    }

    //inizializzo vettori propietà
    vector<UCDProperty<double>> points_properties;
    vector<UCDProperty<double>> polygons_properties;

    VectorXi material_p(polygon_vertices.size());
    for (int i = 0; i < material_p.size(); i++){
        material_p(1) = 1;
    }

    UCDUtilities UCD;
    UCD.ExportPolygons ("polygons_paraview.inp", points_p, polygon_vertices, points_properties, polygons_properties, material_p);


    //tracce
    // numero totale di estremi
    size_t n = t.size();
    size_t num = t.size() * 2;
    //salvo tutti gli estrmei in una matrice
    MatrixXd points_t(3, num) ;
    for (size_t i = 0; i < t.size(); i++){
        points_t.col(i*2) = t[i].vertices.col(0);
        points_t.col(i*2+1) = t[i].vertices.col(1);
    }

    //salvo con lo stesso Id gli estremi di una traccia
    MatrixXi index_edges(2, n);
    for (size_t i = 0; i < n; i++){
        index_edges (0, i) = i * 2;
        index_edges (1, i) = i * 2 + 1;
    }

    //inizializzo vettori propietà
    VectorXi material_t(n);
    for (unsigned int i = 0; i<n; i++) {
        material_t(i) = i;
    }
    UCD.ExportSegments("traces.inp",points_t, index_edges, points_properties, polygons_properties,material_t);

}


void Export_paraview2( vector<MatrixXd> &f){
    cout << "okay" << endl;
    // numero totale di punti
    int N = 0;
    for (const auto&  polygon: f){
        N += (polygon.cols()-2)*3;
    }

    // metto tutti i punti in una grande matrice
    MatrixXd points_p(3,N);
    unsigned int col = 0;
    for (const auto& polygon: f) {
        Vector3d fixed = polygon.col(0);
        for (unsigned int i = 1; i < polygon.cols()-1; i++){
            points_p.col(col) = fixed;
            col++;
            points_p.col(col) = polygon.col(i);;
            col++;
            points_p.col(col) = polygon.col(i+1);
            col++;
        }
    }


    // creo vettore con id raggruppati per poligono
    vector<vector<unsigned int>> polygon_vertices;
    unsigned int start_index = 0;
    for (const auto&  polygon: f){

        //triangolazione del poligono
        vector<Matrix3d> triangolation;
        Vector3d fixed = polygon.col(0);
        for (unsigned int i = 1; i < polygon.cols()-1; i++){

            Matrix3d temp(3,3);
            temp.col(0) = fixed;
            temp.col(1) = polygon.col(i);
            temp.col(2) = polygon.col(i+1);
            triangolation.push_back(temp);
        }

        //stesso id per stesso
        for(unsigned int i = 0; i < triangolation.size(); i++){
            vector<unsigned int> polygon_ids;
            Matrix3d temp = triangolation[i];
            for(unsigned j = 0; j< 3; j++){
                polygon_ids.push_back(start_index);
                start_index ++;
            }
            polygon_vertices.push_back(polygon_ids);
        }
    }

    //inizializzo vettori propietà
    vector<UCDProperty<double>> points_properties;
    vector<UCDProperty<double>> polygons_properties;

    VectorXi material_p(polygon_vertices.size());
    for (int i = 0; i < material_p.size(); ++i){
        material_p(1) = 1;
    }

    UCDUtilities UCD;
    UCD.ExportPolygons ("cutted_paraview.inp", points_p, polygon_vertices, points_properties, polygons_properties, material_p);
}

