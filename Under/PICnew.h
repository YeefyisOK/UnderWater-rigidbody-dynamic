#pragma once
#include <Eigen/Dense>
#include <Eigen/Cholesky>  
#include <Eigen/LU>  
#include <Eigen/QR>  
#include <Eigen/SVD>  
#include <vector>
#include "PIC.h"
#include "iostream"
using namespace std;
using namespace Eigen;
struct VertexandNormalST//顶点坐标和顶点法向
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	Vector3d coordinate;
	Vector3d vertexNormal;
};
struct FaceandNormalST//面索引和面法向
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	int vertexIndex[3];
	Vector3d faceNormal;
};
class PICnew {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	vector <VertexandNormalST> vertexandnormal;
	vector <FaceandNormalST> faceandnormal;
	PICnew(PIC *m_PIC);
	//PICnew();	

};