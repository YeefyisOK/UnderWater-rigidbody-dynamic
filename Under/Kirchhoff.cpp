#include "Kirchhoff.h"
using namespace Eigen;
using namespace std;
CKirchhoff::CKirchhoff(PIC m_pic)
{
	numPoints = m_pic.V.size();
	numFaces = m_pic.F.size();
	for (int i = 0; i < numPoints; i++){
		vertex.row[i] = m_pic.V[i];//��ֵ�������
	}
	for (int i = 0; i < numPoints; i++){
		normal.row[i] = m_pic.VN[i];//��ֵ������� vector��ֵ�����һ��
	}
	for (int k = 0; k < 3; k++){
		for (int i = 0; i < numPoints; i++){
			vector <double>zuobiao;
			zuobiao.push_back(m_pic.V[m_pic.F[i].V[k]].X);
			zuobiao.push_back(m_pic.V[m_pic.F[i].V[k]].Y);
			zuobiao.push_back(m_pic.V[m_pic.F[i].V[k]].Z);
			face[k].row[i] =zuobiao ;//��ֵ����� �������ĳɶ���
		}
	}
}
MatrixXd CKirchhoff::computeKF(double offset){
	MatrixXd K(6, 6);
	MatrixXd S = vertex - offset * normal;
	MatrixXd M = solid_angle(S);
	MatrixXd MF = motion_flux();
	MatrixXd sigma = MF.ldlt().solve(M);//sigma=strength
	//MatrixXd sigma = division(MF, M);
	MatrixXd C = face_center();
	MatrixXd SL = single_layer(S ,C);
	MatrixXd phi = numPoints*SL*sigma;
	MatrixXd Q = one_point_quadrature();
	K = Q*phi;
}

MatrixXd CKirchhoff::single_layer(MatrixXd S , MatrixXd C){
	MatrixXd res(numPoints, 1);
	MatrixXd rr = C - S;
	for (int i = 0; i < numPoints; i++){
		res(i,0) = 1/sqrt( rr.row(i).dot(rr.row(i) ) );
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
	MatrixXd res(numPoints, 3);
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
}

MatrixXd CKirchhoff::triangle_area(){
	MatrixXd res(numFaces, 1);
	MatrixXd AV = area_vector();
	for (int i = 0; i < numPoints; i++){
		res(i, 0) = sqrt(AV.dot(AV));
	}
	return res;
}
MatrixXd CKirchhoff::one_point_quadrature(){
	MatrixXd areas = area_vector();
	MatrixXd VF = face_center();
	MatrixXd NF = face_normal();
	MatrixXd CR(numFaces, 3);
	CR = VF.cross(NF);
	MatrixXd Q(6, numFaces);
	Q.setZero(6, numFaces);
	Q.row(0) = areas*CR.col(0);
	Q.row(1) = areas*CR.col(1);
	Q.row(2) = areas*CR.col(2);
	Q.row(3) = areas*NF.col(3);
	Q.row(4) = areas*NF.col(4);
	Q.row(5) = areas*NF.col(5);
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
		pp23(i, 0) = face[1].row(i).dot(face[2].row(i));
		pp31(i, 0) = face[2].row(i).dot(face[0].row(i));
	}
	//?
	MatrixXd p12(numFaces, 1);
	MatrixXd p23(numFaces, 1);
	MatrixXd p31(numFaces, 1);
	p12 = pp1 + pp2 + pp12;
	p23 = pp2 + pp3 + pp23;
	p31 = pp3 + pp1 + pp31;
	MatrixXd p2p1 = p2-p1;
	MatrixXd p3p2 = p3 - p1;
	MatrixXd p1p3 = p1 - p3;
	MatrixXd t1(numFaces, 1);
	MatrixXd t2(numFaces, 1);
	MatrixXd t3(numFaces, 1);
	t1.col(2) = p12*p2p1.col(2);
	t2.col(2) = p23*p3p2.col(2);
	t3.col(2) = p31*p1p3.col(2);

	t1.col(0) = p12*p2p1.col(0);
	t2.col(0) = p23*p3p2.col(0);
	t3.col(0) = p31*p1p3.col(0);

	t1.col(1) = p12*p2p1.col(1);
	t2.col(1) = p23*p3p2.col(1);
	t3.col(1) = p31*p1p3.col(1);
	MatrixXd res;
	res = -(t1 + t2 + t3) / 6;
	return res;
}
MatrixXd CKirchhoff::area_vector(){
	MatrixXd AV(numFaces,3);
	AV = 0.5*(face[0].cross(face[1]) + face[1].cross(face[2]) + face[2].cross(face[0]));
	return AV;
}
MatrixXd CKirchhoff::solid_angle(MatrixXd src){//numPoints*3
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
MatrixXd CKirchhoff::computeKB(MatrixXd face[], int numFaces, int index[], double mass) {
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
		double x0 = face[0][0];
		double y0 = face[0][1];
		double z0 = face[0][2];

		double x1 = face[1][0];
		double y1 = face[1][1];
		double z1 = face[1][2];

		double x2 = face[2][0];
		double y2 = face[2][1];
		double z2 = face[2][2];
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
		intg[0] += d0*f1x; intg[1] += d0*f2x; intg[2] += d1*f2y; intg[3] += d2*f2z;
		intg[4] += d0*f3x; intg[5] += d1*f3y; intg[6] += d2*f3z;
		intg[7] += d0*(y0*g0x + y1*g1x + y2*g2x); intg[8] += d1*(z0*g0y + z1*g1y + z2*g2y); intg[9] += d2*(x0*g0z + x1*g1z + x2*g2z);
	}
	for (int i = 0; i < 10; i++)
		intg[i] *= mult[i];
	mass = intg[0];
	// center of mass
	cm.setData(intg[1] / mass, intg[2] / mass, intg[3] / mass);
	// inertia tensor relative to center of mass 
	inertia(0, 0) = intg[5] + intg[6] - mass*(cm.getm_data(1)*cm.getm_data(1) + cm.getm_data(2)*cm.getm_data(2));
	inertia(1, 1) = intg[4] + intg[6] - mass*(cm.getm_data(2)*cm.getm_data(2) + cm.getm_data(0)*cm.getm_data(0));
	inertia(2, 2) = intg[4] + intg[5] - mass*(cm.getm_data(0)*cm.getm_data(0) + cm.getm_data(1)*cm.getm_data(1));
	inertia(0, 1) = inertia(1, 0) = -(intg[7] - mass*cm.getm_data(0)*cm.getm_data(1));
	inertia(1, 2) = inertia(2, 1) = -(intg[8] - mass*cm.getm_data(1)*cm.getm_data(2));
	inertia(2, 0) = inertia(0, 2) - (intg[9] - mass*cm.getm_data(2)*cm.getm_data(0));
	return inertia;
}
CVector6D CKirchhoff::computeK(){
	MatrixXd KF, KB;
	MatrixXd temp1 = KF.rowwise().sum();
	MatrixXd temp2 = KB.rowwise().sum();
	temp1 = temp1 + temp2;
	double a = temp1(0, 0);
	double b = temp1(1, 0);
	double c = temp1(2, 0);
	double d = temp1(3, 0);
	double e = temp1(4, 0);
	double f = temp1(5, 0);
	CVector6D Kirchhoff(a,b,c,d,e,f);
	return Kirchhoff;
}