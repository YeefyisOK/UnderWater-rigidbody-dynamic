#include "Kirchhoff.h"
using namespace Eigen;
using namespace std;
CKirchhoff::CKirchhoff(PIC m_pic)
{
	numPoints = m_pic.V.size();
	numFaces = m_pic.F.size();
	vertex.resize(numPoints, 3);
	normal.resize(numFaces, 3);
	face[0].resize(numFaces, 3);
	face[1].resize(numFaces, 3);
	face[2].resize(numFaces, 3);
	if (m_pic.V.size() > 0) {
		for (int i = 0; i < numPoints; i++) {
			VectorXd temp(3);
			temp(0) = m_pic.V[i].X;
			temp(1) = m_pic.V[i].Y;
			temp(2) = m_pic.V[i].Z;		
			vertex.row(i) = temp;//赋值顶点矩阵
		}
	}
	if (m_pic.F.size() > 0) {
		for (int i = 0; i < numPoints; i++){
			
			VectorXd temp(3);
			temp(0) = m_pic.VN[m_pic.F[i].N[0]].NX ;
			temp(1) = m_pic.VN[m_pic.F[i].N[0]].NY;
			temp(2) = m_pic.VN[m_pic.F[i].N[0]].NZ;
			normal.row(i) = temp;//赋值法向矩阵 vector赋值矩阵的一行
		}
	}
	if (m_pic.F.size() > 0) {
		for (int k = 0; k < 3; k++){//第几行
			for (int i = 0; i < numPoints; i++){
				/*
				vector <double>zuobiao;
				zuobiao.push_back(m_pic.V[m_pic.F[i].V[k]].X);
				zuobiao.push_back(m_pic.V[m_pic.F[i].V[k]].Y);
				zuobiao.push_back(m_pic.V[m_pic.F[i].V[k]].Z);*/
				VectorXd temp(3);
				temp(0) = m_pic.V[m_pic.F[i].V[k]].X;
				temp(1) = m_pic.V[m_pic.F[i].V[k]].Y;
				temp(2) = m_pic.V[m_pic.F[i].V[k]].Z;
				face[k].row(i) = temp;//赋值面矩阵 把索引改成顶点
			}
		}
	}

}
MatrixXd CKirchhoff::computeKF(double offset){
	MatrixXd K(6, 6);
	MatrixXd C = face_center();
	MatrixXd S = C - offset * normal;
	MatrixXd M = solid_angle(S);
	MatrixXd FL = motion_flux();
	MatrixXd sigma(numFaces, 6);
	for (int i = 0; i< 6;i++) {

		sigma.col(i) = M.colPivHouseholderQr().solve(FL.col(i));//sigma=strength NAN!
	}
	//MatrixXd sigma = division(MF, M);
	MatrixXd SL = single_layer(S ,C);//inf！
	MatrixXd phi = SL* sigma;// numFaces * SL * sigma;
	MatrixXd Q = one_point_quadrature();
	K = Q*phi;
	return K;
}

MatrixXd CKirchhoff::single_layer(MatrixXd S , MatrixXd C){
	MatrixXd res(numFaces, numFaces);
	for (int i = 0; i < numFaces; i++){//C Z
		for (int j = 0;j < numFaces;j++) {//S
			VectorXd temp = C.row(i) - S.row(j);
			double isZero = temp.dot(temp);
			if (isZero<0.05 || isZero>-0.05)
				isZero = 0.05;
			res(i,j) = 1.0/sqrt(isZero);
		}
	}
	return res;
}


MatrixXd CKirchhoff::face_normal(){
	MatrixXd AV = area_vector();
	MatrixXd areas = triangle_area();
	for (int i = 0; i < 3; i++){
		AV.col(i) = AV.col(i).cwiseQuotient(areas);
	}
	return AV;
}
MatrixXd CKirchhoff::face_center(){
	MatrixXd res(numFaces, 3);
	res = (face[0] + face[1] + face[2])/3;
	return res;
}
double CKirchhoff::area(){
	MatrixXd temp = triangle_area();
	double res= temp.colwise().sum()(0,0);
	return res;
}
MatrixXd CKirchhoff::motion_flux(){
	MatrixXd res(numFaces, 6);
	res.leftCols(3) = angular_vector();
	res.rightCols(3) = area_vector();
	return res;
}

MatrixXd CKirchhoff::triangle_area(){
	MatrixXd res(numFaces, 1);
	MatrixXd AV = area_vector();//(numFaces, 3);
	for (int i = 0; i < numPoints; i++){
		res(i, 0) = sqrt(AV.row(i).dot( AV.row(i) ));
	}
	return res;
}
MatrixXd CKirchhoff::one_point_quadrature(){
	MatrixXd areas = triangle_area();
	MatrixXd VF = face_center();
	MatrixXd NF = face_normal();
	MatrixXd CR(numFaces, 3);
	for (int i = 0;i < numFaces;i++) {
		Vector3d tempVF = VF.row(i);
		Vector3d tempNF = NF.row(i);
		CR.row(i) = tempVF.cross(tempNF);
	}
	MatrixXd Q(6, numFaces);
	Q.row(0) = (areas.array()*CR.col(0).array()).matrix().transpose();
	Q.row(1) = (areas.array()*CR.col(1).array()).matrix().transpose();
	Q.row(2) = (areas.array()*CR.col(2).array()).matrix().transpose();
	Q.row(3) = (areas.array()*NF.col(0).array()).matrix().transpose();
	Q.row(4) = (areas.array()*NF.col(1).array()).matrix().transpose();
	Q.row(5) = (areas.array()*NF.col(2).array()).matrix().transpose();
	return Q;
}
MatrixXd CKirchhoff::angular_vector(){
	MatrixXd p1 = face[0];
	MatrixXd p2 = face[1];
	MatrixXd p3 = face[2];
	MatrixXd pp1(numFaces, 1);
	MatrixXd pp2(numFaces, 1);
	MatrixXd pp3(numFaces, 1);
	for (int i = 0; i < numFaces; i++){
		pp1(i, 0) = face[0].row(i).dot(face[0].row(i));
		pp2(i, 0) = face[1].row(i).dot(face[1].row(i));
		pp3(i, 0) = face[2].row(i).dot(face[2].row(i));
	}
	MatrixXd pp12(numFaces, 1);
	MatrixXd pp23(numFaces, 1);
	MatrixXd pp31(numFaces, 1);
	for (int i = 0; i < numFaces; i++){
		pp12(i, 0) = face[0].row(i).dot(face[1].row(i));
		pp23(i, 0) = face[2].row(i).dot(face[1].row(i));
		pp31(i, 0) = face[2].row(i).dot(face[0].row(i));
	}
	//?
	/*
	MatrixXd p12(numFaces, 1);
	MatrixXd p23(numFaces, 1);
	MatrixXd p31(numFaces, 1);*/
	MatrixXd p12= pp1 + pp2 + pp12;
	MatrixXd p23= pp2 + pp3 + pp23;
	MatrixXd p31= pp3 + pp1 + pp31;

	MatrixXd p2p1 = p2 - p1;
	MatrixXd p3p2 = p3 - p2;
	MatrixXd p1p3 = p1 - p3;
	MatrixXd t1(numFaces, 3);
	MatrixXd t2(numFaces, 3);
	MatrixXd t3(numFaces, 3);
	for (int i = 0;i < numFaces;i++) {

		t1(i, 2) = p12(i, 0)*p2p1.col(2)(i, 0);//abort!!!
		t2(i, 2) = p23(i, 0)*p3p2.col(2)(i, 0);
		t3(i, 2) = p31(i, 0)*p1p3.col(2)(i, 0);

		t1(i, 0) = p12(i, 0)*p2p1.col(0)(i, 0);
		t2(i, 0) = p23(i, 0)*p3p2.col(0)(i, 0);
		t3(i, 0) = p31(i, 0)*p1p3.col(0)(i, 0);

		t1(i, 1) = p12(i, 0)*p2p1.col(1)(i, 0);
		t2(i, 1) = p23(i, 0)*p3p2.col(1)(i, 0);
		t3(i, 1) = p31(i, 0)*p1p3.col(1)(i, 0);
	}
	MatrixXd res= -(t1 + t2 + t3) / 6;
	return res;
}
MatrixXd CKirchhoff::area_vector(){
	MatrixXd AV(numFaces,3);
	for (int i = 0;i < numFaces;i++) {
		Vector3d tempface0 = face[0].row(i);
		Vector3d tempface1 = face[1].row(i);
		Vector3d tempface2 = face[2].row(i);

		AV.row(i) = 0.5*(tempface0.cross(tempface1) + tempface0.cross(tempface2)
			+ tempface2.cross(tempface0));
	}
	return AV;
}
MatrixXd CKirchhoff::solid_angle(MatrixXd src){//numPoints*3
	MatrixXd res(numFaces, numFaces);
//	MatrixXd FacePoint(numFaces, 3);
	for (int i = 0; i < numFaces; i++){
		for (int j = 0; j < numFaces; j++){
			Vector3d R1 = face[0].row(i) - src.row(j);
			Vector3d R2 = face[1].row(i) - src.row(j);
			Vector3d R3 = face[2].row(i) - src.row(j);//向量还是矩阵可以dot
			Matrix3d temp;
			temp.row(0) = R1;
			temp.row(1) = R2;
			temp.row(2) = R3;
			double N = R1.dot(R2.cross(R3));
			double l1 = sqrt(R1.dot(R1));
			double l2 = sqrt(R2.dot(R2));
			double l3 = sqrt(R3.dot(R3));
			double Den = l1*l2*l3 + l1*R2.dot(R3) + l2*R3.dot(R1) + l3*R1.dot(R2);
			res(i, j) = 2 * atan(N / Den);
		}
	}
	return res;
}


void CKirchhoff::Subexpressions(double &w0, double &w1, double &w2, double &f1, double &f2, double &f3, double &g0, double &g1, double &g2) {
	double temp0 = w0 + w1;
	f1 = temp0 + w2;
	double temp1 = w0*w0;
	double temp2 = temp1 + w1*temp0;
	f2 = temp2 + w2*f1;
	f3 = w0*temp1 + w1*temp2 + w2*f2;
	g0 = f2 + w0*(f1 + w0);
	g1 = f2 + w1*(f1 + w1);
	g2 = f2 + w2*(f1 + w2);
}
MatrixXd CKirchhoff::computeKB() {
	double mass = 0;
	CPoint3D cm(0, 0, 0);
	MatrixXd inertia(6, 6);
	const double mult[10] = { 1 / 6, 1 / 24, 1 / 24, 1 / 24, 1 / 60, 1 / 60, 1 / 60, 1 / 120, 1 / 120, 1 / 120 };
	double intg[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; // order: 1, x, y, z, x^2, y^2, z^2, xy, yz, zx 
	for (int t = 0; t < numFaces; t++) { // get vertices of triangle t 
		/*
		i0 = index[3*t]; 	i1 = index[3*t+1];		i2 = index[3*t+2];
		x0 = p[i0].x; 		y0 = p[i0].y; 		z0 = p[i0].z;
		x1 = p[i1].x;       y1 = p[i1].y;       z1 = p[i1].z;
		x2 = p[i2].x;       y2 = p[i2].y;       z2 = p[i2].z;
		*/
		double x0 = face[0](t,0);
		double y0 = face[0](t,1);
		double z0 = face[0](t,2);

		double x1 = face[1](t,0);
		double y1 = face[1](t,1);
		double z1 = face[1](t,2);

		double x2 = face[2](t,0);
		double y2 = face[2](t,1);
		double z2 = face[2](t,2);
		// get edges and cross product of edges 
		double a1 = x1 - x0; double  b1 = y1 - y0; double  c1 = z1 - z0;
		double a2 = x2 - x0; double b2 = y2 - y0; double c2 = z2 - z0;
		double d0 = b1*c2 - b2*c1; double d1 = a2*c1 - a1*c2; double d2 = a1*b2 - a2*b1;
		// compute integral terms
		double f1x, f2x, f3x;
		double f1y, f2y, f3y;
		double f1z, f2z, f3z;
		double g0x, g1x, g2x;
		double g0y, g1y, g2y;
		double g0z, g1z, g2z;
		Subexpressions(x0, x1, x2, f1x, f2x, f3x, g0x, g1x, g2x);
		Subexpressions(y0, y1, y2, f1y, f2y, f3y, g0y, g1y, g2y);
		Subexpressions(z0, z1, z2, f1z, f2z, f3z, g0z, g1z, g2z);
		// update integrals 
		//const double mult[10] = { 1 / 6, 1 / 24, 1 / 24, 1 / 24, 1 / 60, 1 / 60, 1 / 60, 1 / 120, 1 / 120, 1 / 120 };

		intg[0] += d0*f1x/6; intg[1] += d0*f2x/24; intg[2] += d1*f2y / 24; intg[3] += d2*f2z / 24;
		intg[4] += d0*f3x / 60; intg[5] += d1*f3y / 60; intg[6] += d2*f3z / 60;
		intg[7] += d0*(y0*g0x + y1*g1x + y2*g2x) / 120; 
		intg[8] += d1*(z0*g0y + z1*g1y + z2*g2y) / 120; 
		intg[9] += d2*(x0*g0z + x1*g1z + x2*g2z) / 120;
	}
	mass = intg[0];
	// 质心
	cm.setData(intg[1] / mass, intg[2] / mass, intg[3] / mass);
	// 相对于质心的惯性张量
	inertia(0, 0) = intg[5] + intg[6] - mass*(cm.getm_data(1)*cm.getm_data(1) + cm.getm_data(2)*cm.getm_data(2));
	inertia(1, 1) = intg[4] + intg[6] - mass*(cm.getm_data(2)*cm.getm_data(2) + cm.getm_data(0)*cm.getm_data(0));
	inertia(2, 2) = intg[4] + intg[5] - mass*(cm.getm_data(0)*cm.getm_data(0) + cm.getm_data(1)*cm.getm_data(1));
	inertia(0, 1) = inertia(1, 0) = -(intg[7] - mass*cm.getm_data(0)*cm.getm_data(1));
	inertia(1, 2) = inertia(2, 1) = -(intg[8] - mass*cm.getm_data(1)*cm.getm_data(2));
	inertia(2, 0) = inertia(0, 2) - (intg[9] - mass*cm.getm_data(2)*cm.getm_data(0));
	return inertia;
}
CVector6D CKirchhoff::computeK(){
	MatrixXd KB = computeKB();//mass
	//MatrixXd KF = computeKF(0.4);//offest
	MatrixXd temp1 = KB.rowwise().sum();
	//MatrixXd temp2 = KF.rowwise().sum();
	//temp1 = temp1 + temp2;
	double a = temp1(0, 0);
	double b = temp1(1, 0);
	double c = temp1(2, 0);
	double d = temp1(3, 0);
	double e = temp1(4, 0);
	double f = temp1(5, 0);
	CVector6D Kirchhoff(a,b,c,d,e,f);
	return Kirchhoff;
}