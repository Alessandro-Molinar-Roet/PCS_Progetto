#pragma once

#include <gtest/gtest.h>
#include <math.h>
#include <Eigen/Eigen>
#include <vector>
#include "Geometry.hpp"
#include "Utils.hpp"
#include "FracturesNetworkLibrary.hpp"
#include "tol.hpp"

using namespace Eigen;
using namespace std;
namespace FractureNetwork
{
// ***************************************************************************

TEST(lineFractIntersectTest, lineFractIntersectTest_intersection_Test)
{
    Matrix3Xd fracture1(3,4);
    fracture1 << 0.0, 1.0, 1.0, 0.0,
        0.0, 0.0, 1.0, 1.0,
        0.0, 0.0, 0.0, 0.0;
    Vector3d pointOnLine(0.8, 0.0, 0.0);
    Vector3d direction(0.0, -1.0, 0.0);
    vector<Vector3d> intersection1;
    Vector3d firstPointIntersect (0.8, 0.0, 0.0);
    Vector3d secondPointIntersect (0.8, 1.0, 0.0);

    bool result = lineFractIntersect(fracture1, pointOnLine, direction, intersection1);

    ASSERT_TRUE(result);
    ASSERT_TRUE(intersection1[0] == firstPointIntersect && intersection1[1] == secondPointIntersect);
}

// ***************************************************************************
// ***************************************************************************

TEST(normalPTest, normalPTest_norma_Test)
{
    Matrix3Xd fracture(3,4);
    fracture << 0.0, 1.0, 1.0, 0.0,
        0.0, 0.0, 1.0, 1.0,
        0.0, 0.0, 0.0, 0.0;
    Vector3d normal (0.0, 0.0, 1.0);
    Vector3d err_normal (1.0, 0.0, 1.0);
    ASSERT_TRUE(normal == normalP(fracture));
    ASSERT_FALSE(err_normal == normalP(fracture));
}

// ***************************************************************************
// ***************************************************************************

TEST(near1Test, near1Test_boundingbox_Test)
{
    Matrix3Xd fracture1(3,4);
    fracture1 << 0.0, 1.0, 1.0, 0.0,
        0.0, 0.0, 1.0, 1.0,
        0.0, 0.0, 0.0, 0.0;
    Matrix3Xd fracture2(3,4);
    fracture2 << 0.8, 0.8, 0.8, 0.8,
        0.0, 0.0, 1.0, 1.0,
        -0.1, 0.3, 0.3, -0.1;
    Matrix3Xd fracture3(3,4);
    fracture3 << -0.237778, 0.3161837, 0.3161837, -0.237778,
        0.5, 0.5, 0.5, 0.5,
        -0.34444, -0.34444, 0.4528389, 0.4528389;

    ASSERT_TRUE(near1(fracture1, fracture2));
    ASSERT_TRUE(near1(fracture1, fracture3));
    ASSERT_FALSE(near1(fracture2, fracture3));

    Matrix3Xd fracture1toll(3,4);
    fracture1toll << 0.0 +tol, 1.0 +tol, 1.0 +tol, 0.0 +tol,
        0.0 +tol, 0.0 +tol, 1.0 +tol, 1.0 +tol,
        0.0 +tol, 0.0 +tol, 0.0 +tol, 0.0 +tol;

    ASSERT_TRUE(near1(fracture1toll, fracture2));
    ASSERT_TRUE(near1(fracture1toll, fracture3));
}

// ***************************************************************************
// ***************************************************************************

TEST(near2Test, near2Test_circumference_Test)
{
    Matrix3Xd fracture1(3,4);
    fracture1 << 0.0, 1.0, 1.0, 0.0,
        0.0, 0.0, 1.0, 1.0,
        0.0, 0.0, 0.0, 0.0;
    Matrix3Xd fracture2(3,4);
    fracture2 << 0.8, 0.8, 0.8, 0.8,
        0.0, 0.0, 1.0, 1.0,
        -0.1, 0.3, 0.3, -0.1;
    Matrix3Xd fracture3(3,4);
    fracture3 << -0.237778, 0.3161837, 0.3161837, -0.237778,
        0.5, 0.5, 0.5, 0.5,
        -0.34444, -0.34444, 0.4528389, 0.4528389;

    ASSERT_TRUE(near2(fracture1, fracture2));
    ASSERT_TRUE(near2(fracture1, fracture3));

    Matrix3Xd fracture1toll(3,4);
    fracture1toll << 0.0 +tol, 1.0 +tol, 1.0 +tol, 0.0 +tol,
        0.0 +tol, 0.0 +tol, 1.0 +tol, 1.0 +tol,
        0.0 +tol, 0.0 +tol, 0.0 +tol, 0.0 +tol;

    ASSERT_TRUE(near2(fracture1toll, fracture2));
    ASSERT_TRUE(near2(fracture1toll, fracture3));
}
// ***************************************************************************

// Test in cui si controlla che alcune intersezioni posso essere viste da near1 ma non da near 2
TEST(nearDifferenceTest, nearDifferenceTest_near_Test)
{
    Matrix3Xd first_fracture(3,4);
    first_fracture << 100.0, -100.0, 0.0, 0.0,
        0.0, 0.0, 100.0, -100.0,
        0.0, 0.0, 0.0, 0.0;
    Matrix3Xd second_fracture(3,4);
    second_fracture << 100.0, -100.0, 0.0, 0.0,
        0.0, 0.0, 100.0, -100.0,
        0.5, 0.5, 0.5, 0.5;

    ASSERT_FALSE(near1(first_fracture, second_fracture));
    ASSERT_TRUE(near2(first_fracture, second_fracture));
}
// ***************************************************************************
void split(const MatrixXd& polygon, const vector<unsigned int>& all, const vector<Trace>& tracce, unsigned int counter, list<MatrixXd>& cutted){
    //inizializzazione variabili
    unsigned int current;
    Vector3d point1 = {};
    Vector3d point2 = {};
    Vector3d direzione = {};
    bool full = false;

    vector<Vector3d> intersezioni = {};
    vector<unsigned int> lati = {};

    if(counter < all.size()){
        current = all[counter]; // ID traccia tagliante
        point1 = tracce[current].vertices.col(0); // estremo 1
        point2 = tracce[current].vertices.col(1); // estremo 2
        direzione = point2-point1; // direzione dell traccia

        full = IsInside(polygon ,point1 ,direzione ,intersezioni, lati); // true se la frattura contiene l'estremo uno della traccia

        counter++;
        if( !full || lati[0] == lati[1]){
            split(polygon,all,tracce,counter,cutted);
        }
        else{
            if(intersezioni.size() != 2){
                cout << "error: " << intersezioni.size() << endl;
            }
            else{
                unsigned int dim = polygon.cols();
                MatrixXd temp_sx(3,dim+2);
                MatrixXd temp_dx(3,dim+2);

                bool change = false;
                unsigned int counter1 = 0;
                unsigned int counter2 = 0;

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
        }
    }
    else{
        cutted.push_back(polygon);
    }
}
// ***************************************************************************
bool IsInside(const MatrixXd& frattura, const Vector3d& puntoRetta,const Vector3d& direzione,vector<Vector3d>& intersezioni, vector<unsigned int>& lato){
    unsigned int counter = 1;
    unsigned int meet1 = 0;
    unsigned int meet2 = 0;

    for(unsigned int i= 0; i<frattura.cols();i++){
        Vector3d puntoipoly = frattura.col(i);
        Vector3d latoipoly =frattura.col((i+1) % frattura.cols())-frattura.col(i);

        if(!(latoipoly.cross(direzione)).isZero(tol)){
            double t = ((puntoRetta.cross(direzione)).dot((latoipoly.cross(direzione)))-(puntoipoly.cross(direzione)).dot((latoipoly.cross(direzione))))/(latoipoly.cross(direzione).norm()*latoipoly.cross(direzione).norm());

            if ( (t > (-tol)) && (t < (1.0+tol)) ){
                intersezioni.push_back(puntoipoly + t * latoipoly);
                lato.push_back(counter);

                double k = (puntoipoly + t*latoipoly -puntoRetta).dot(direzione)/(direzione.norm()*direzione.norm());
                if(k > tol){
                    meet1++;
                    if(k<1.0+tol){
                        meet2++;
                    }
                }
            }
        }
        else{
            if((puntoipoly-puntoRetta).cross(direzione).isZero(tol)){
                //ATTENZIONE CONTROLLO COMPRESO
                cout << "fakes" << endl;
                lato.push_back(counter);
            }
        }
        counter++; //passo al lato successivo
    }
    if(meet1 == 1 || meet2 == 1){
        return true;
    }
    return false;
}
// ***************************************************************************
TEST(splitTest, splitTest_split_Test)
{
    // FractureNetwork::Fracture fracture1;
    // unsigned int i = 0;
    // fracture1.ID = i;
    Matrix3Xd fracture_split(3,4);
    fracture_split << 0.0, 1.0, 0.0, 1.0,
                      0.0, 1.0, 1.0, 0.0,
                      0.0, 0.0, 0.0, 0.0;
    vector<unsigned int> all;
    all.push_back(0);

    vector<Trace> tracce;
    Trace tr;
    tr.ID = 0;
    tr.length = 1.0;
    tr.vertices << 0.5, 0.5,
                   0.0, 0.5,
                   0.0, 0.0;
    tr.first_generator = 0;
    tr.second_generator = 1;

    tracce.push_back(tr);

    unsigned int counter = 0;
    list<MatrixXd> cutted;

    list<MatrixXd> result;
    MatrixXd matrice1;
    matrice1 << 0.0, 0.5, 0.5, 0.0,
                0.0, 0.0, 0.5, 1.0,
                0.0, 0.0, 0.0, 0.0;
    result.push_back(matrice1);
    MatrixXd matrice2;
    matrice2 << 0.5, 1.0, 1.0, 0.5,
                0.0, 0.0, 1.0, 0.5,
                0.0, 0.0, 0.0, 0.0;
    result.push_back(matrice1);

    split(fracture_split, all, tracce, counter, cutted);

    ASSERT_TRUE(result == cutted);

}
}
