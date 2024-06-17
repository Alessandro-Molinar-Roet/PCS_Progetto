#pragma once

#include <gtest/gtest.h>
#include <Eigen/Eigen>
#include <vector>
#include "../src/Geometry.hpp"
#include "../src/FracturesNetworkLibrary.hpp"
#include "../src/tol.hpp"

using namespace Eigen;
using namespace std;
namespace FractureNetwork{

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
    fracture1toll << 0.0 +tol1, 1.0 +tol1, 1.0 +tol1, 0.0 +tol1,
        0.0 +tol1, 0.0 +tol1, 1.0 +tol1, 1.0 +tol1,
        0.0 +tol1, 0.0 +tol1, 0.0 +tol1, 0.0 +tol1;

    ASSERT_TRUE(near1(fracture1toll, fracture2));
    ASSERT_TRUE(near1(fracture1toll, fracture3));
}


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
    fracture1toll << 0.0 +tol1, 1.0 +tol1, 1.0 +tol1, 0.0 +tol1,
        0.0 +tol1, 0.0 +tol1, 1.0 +tol1, 1.0 +tol1,
        0.0 +tol1, 0.0 +tol1, 0.0 +tol1, 0.0 +tol1;

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
TEST(splitTest, splitTest_split_Test)
{
    Matrix3Xd fracture(3,4);
    fracture << 0.0, 1.0, 1.0, 0.0,
                      0.0, 0.0, 1.0, 1.0,
                      0.0, 0.0, 0.0, 0.0;
    vector<unsigned int> all = {0};

    Trace tr;
    tr.ID = 0;
    tr.length = 1.0;
    tr.vertices << 0.5, 0.5,
                   0.0, 0.5,
                   0.0, 0.0;
    tr.first_generator = 0;
    tr.second_generator = 1;
    vector<Trace> tracce = {tr};

    unsigned int counter = 0;
    list<MatrixXd> cutted = {};

    MatrixXd matrice1(3,4);
    matrice1 << 0.0, 0.5, 0.5, 0.0,
                0.0, 0.0, 1.0, 1.0,
                0.0, 0.0, 0.0, 0.0;
    MatrixXd matrice2(3,4);
    matrice2 << 1.0, 1.0, 0.5, 0.5,
                0.0, 1.0, 1.0, 0.0,
                0.0, 0.0, 0.0, 0.0;
    list<MatrixXd> result = {matrice1, matrice2};

    split(fracture, all, tracce, counter, cutted);
    EXPECT_EQ(cutted, result);

    //Caso con 2 tracce non passanti
    vector<unsigned int> all_2 = {0, 1};

    Trace tr1;
    Trace tr2;
    tr1.ID = 0;
    tr2.ID = 1;
    tr1.length = 0.5;
    tr2.length = 0.25;

    tr1.vertices << 0.25, 0.75,
        0.5, 0.5,
        0.0, 0.0;
    tr2.vertices << 0.875, 0.875,
        0.375, 0.625,
        0.0, 0.0;

    tr1.first_generator = 0;
    tr1.second_generator = 1;
    tr2.first_generator = 2;
    tr2.second_generator = 3;

    vector<Trace> tracce_2 = {tr1, tr2};
    list<MatrixXd> cutted_2 = {};

    MatrixXd matrice3(3,4);
    matrice3 << 0.0, 0.875, 0.875, 0.0,
        0.0, 0.0, 0.5, 0.5,
        0.0, 0.0, 0.0, 0.0;
    MatrixXd matrice4(3,4);
    matrice4 << 1.0, 1.0, 0.875, 0.875,
        0.0, 0.5, 0.5, 0.0,
        0.0, 0.0, 0.0, 0.0;
    MatrixXd matrice5(3,4);
    matrice5 << 1.0, 0.875, 0.875, 1.0,
        1.0, 1.0, 0.5, 0.5,
        0.0, 0.0, 0.0, 0.0;
    MatrixXd matrice6(3,4);
    matrice6 << 0.0, 0.0, 0.875, 0.875,
        1.0, 0.5, 0.5, 1.0,
        0.0, 0.0, 0.0, 0.0;

    list<MatrixXd> result_2 = {matrice3, matrice4, matrice5, matrice6};

    split(fracture, all_2, tracce_2, counter, cutted_2);

    EXPECT_EQ(cutted_2, result_2);
}

}
