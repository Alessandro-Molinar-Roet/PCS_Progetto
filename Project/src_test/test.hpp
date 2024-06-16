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
    // FractureNetwork::Fracture fracture1;
    // unsigned int i = 0;
    // fracture1.ID = i;
    Matrix3Xd fracture_split(3,4);
    fracture_split << 0.0, 1.0, 1.0, 0.0,
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

    split(fracture_split, all, tracce, counter, cutted);
    EXPECT_EQ(cutted, result);
}

}
