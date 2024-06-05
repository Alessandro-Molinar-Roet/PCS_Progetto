//ATTENZIONE ESEMPI VICINI DA CANCELLARE:
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef __example_unit_test_HPP__
#define __example_unit_test_HPP__

#include <gtest/gtest.h>
#include "../src/Geometry.hpp"
#include "../src/FracturesNetworkLibrary.hpp"
#include "../src/tol.hpp"
#include "../src/Utils.hpp"
#include <math.h>
#include <Eigen/Dense>
#include <vector>
using namespace std;
using namespace Eigen;
using namespace FractureNetwork;

// ***************************************************************************
double squareRoot(const double a)
{
    const double b = sqrt(a);

    if (b != b) // NaN check
        return -1.0;
    else
        return b;
}
// ***************************************************************************
// bool lineFractIntersect(const MatrixXd& fracture, const Vector3d& pointOnLine,const Vector3d& direction){
//     vector<Vector3d> intersections;
//     for(unsigned int i=0;i<fracture.cols();i++){
//         //intersezione lato-retta
//         Vector3d ithFractureVertex = fracture.col(i);
//         Vector3d ithFractureEdge =fracture.col((i+1)%fracture.cols())-fracture.col(i);

//         double t = ((pointOnLine.cross(direction)).dot((ithFractureEdge.cross(direction)))-(ithFractureVertex.cross(direction)).dot((ithFractureEdge.cross(direction))))/(ithFractureEdge.cross(direction).norm()*ithFractureEdge.cross(direction).norm());
//         if (t >= 0 && t <= 1.0){
//             intersections.push_back(ithFractureVertex + t * ithFractureEdge);
//         }
//         if(intersections.size()==2){
//             return true;
//         }
//     }
//     return false;
// }
//****************************************************************************

// TEST(lineFractIntersectTest, intersectionTest) {
//     Matrix3Xd fracture1(3, 4);
//     fracture1 << 0.0, 1.0, 1.0, 0.0,
//         0.0, 0.0, 1.0, 1.0,
//         0.0, 0.0, 0.0, 0.0;
//     Vector3d pointOnLine(0.8, 0.0, 0.0);
//     Vector3d direction(0.0, -1.0, 0.0);
//     vector<Vector3d> intersections;

//     // Call the function and check the result
//     bool result = lineFractIntersect(fracture1, pointOnLine, direction,intersections);
//     ASSERT_TRUE(result);  // Use ASSERT_TRUE for boolean checks
// }


// TEST(lineFractIntersectTest, intersecati){
//     Matrix3Xd fracture1(3, 4)<< 0.0, 1.0, 1.0, 0.0,
//                  0.0, 0.0, 1.0, 1.0,
//                  0.0, 0.0, 0.0, 0.0;
//     ASSERT_EQ(true, lineFractIntersect(fracture1,[0.8,0.0,0.0],[0.0,-1,0.0]));
// }
//****************************************************************************
TEST(SquareRootTest, PositiveNos)
{
    ASSERT_EQ(6, squareRoot(36.0));
    ASSERT_EQ(18.0, squareRoot(324.0));
    ASSERT_EQ(25.4, squareRoot(645.16));
    ASSERT_EQ(0, squareRoot(0.0));
}
// ***************************************************************************
TEST(SquareRootTest, NegativeNos)
{
    ASSERT_EQ(-1.0, squareRoot(-15.0));
    ASSERT_EQ(-1.0, squareRoot(-0.2));
}
// ***************************************************************************

#endif // __example_unit_test_HPP__
