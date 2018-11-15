#include "KirchhoffF.h"
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <math.h>
#include <vector>
using namespace Eigen;
using namespace std;
CKirchhoffF::CKirchhoffF(PIC m_pic)
{
	numPoints = m_pic.V.size();
	numFaces = m_pic.F.size();
	for (int i = 0; i < numPoints; i++){
		vertex.row[i] = m_pic.V[i];//赋值顶点矩阵
	}
	for (int i = 0; i < numPoints; i++){
		normal.row[i] = m_pic.VN[i];//赋值法向矩阵 vector赋值矩阵的一行
	}
	for (int k = 0; k < 3; k++){
		for (int i = 0; i < numPoints; i++){
			vector <double>zuobiao;
			zuobiao.push_back(m_pic.V[m_pic.F[i].V[k]].X);
			zuobiao.push_back(m_pic.V[m_pic.F[i].V[k]].Y);
			zuobiao.push_back(m_pic.V[m_pic.F[i].V[k]].Z);
			face[k].row[i] =zuobiao ;//赋值面矩阵 把索引改成顶点
		}
	}
}
MatrixXd CKirchhoffF::compute(double offset){
	MatrixXd K(6, 6);
	MatrixXd S = vertex - offset * normal;
	MatrixXd M = solid_angle(S);
	MatrixXd MF = motion_flux();
	MatrixXd sigma = MF.ldlt().solve(M);
	//MatrixXd sigma = division(MF, M);
	MatrixXd S =
}
double CKirchhoffF::surf_area(){

}
MatrixXd CKirchhoffF::motion_flux(){

}
MatrixXd CKirchhoffF::one_point_quadrature(){

}
MatrixXd CKirchhoffF::angular_vector(){

}
MatrixXd CKirchhoffF::area_vector(){
	MatrixXd AV(numFaces,3);
	AV = 0.5*(face[0].cross(face[1]) + face[1].cross(face[2]) + face[2].cross(face[0]));
	return AV;
}
MatrixXd CKirchhoffF::solid_angle(MatrixXd src){//numPoints*3
	MatrixXd res(numFaces, numPoints);
//	MatrixXd FacePoint(numFaces, 3);
	for (int i = 0; i < numFaces; i++){
		for (int j = 0; j < numPoints; j++){
			VectorXd R1 = face[0].row(i) - src.row(j);
			VectorXd R2 = face[1].row(i) - src.row(j);
			VectorXd R3 = face[2].row(i) - src.row(j);
			MatrixXd temp(3, 3);
			temp.row[0] = R1;
			temp.row[1] = R2;
			temp.row[2] = R3;
			double N = R1.dot(R2.cross(R3));
			double l1 = sqrt(R1.dot(R1));
			double l2 = sqrt(R2.dot(R2));
			double l3 = sqrt(R3.dot(R3));
			double Den = l1*l2*l3 + l1*R2.dot(R3) + l2*R3.dot(R1) + l3*R1.dot(R2);
			//den = l1.*l2.*l3 + l1.*dot(p2', p3') + l2.*dot(p3', p1') + l3.*dot(p1', p2');
			res(i, j) = 2 * atan(N / Den);

		}
	}
	return res;
}