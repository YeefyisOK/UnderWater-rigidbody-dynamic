#include "Kirchhoff.h"
#include<iostream>
using namespace Eigen;
using namespace std;
CKirchhoff::CKirchhoff(PIC m_pic)
{
	numPoints = m_pic.V.size();
	numFaces = m_pic.F.size();
	vertex.resize(numPoints, 3);
	normal.resize(numPoints, 3);
	face[0].resize(numFaces, 3);
	face[1].resize(numFaces, 3);
	face[2].resize(numFaces, 3);
	if (m_pic.V.size() > 0) {
#pragma omp parallel
		for (int i = 0; i < numPoints; i++) {
			Vector3f temp(3);
			temp(0) = m_pic.V[i].X;
			temp(1) = m_pic.V[i].Y;
			temp(2) = m_pic.V[i].Z;		
			vertex.row(i) = temp.transpose();//赋值顶点矩阵
		}
	}
	if (m_pic.F.size() > 0) {
#pragma omp parallel
		for (int j = 0;j < numPoints; j++) {
			Vector3f sum(0,0,0);
			int k = 0;
#pragma omp parallel
			for (int i = 0; i < numFaces; i++){
				if (m_pic.F[i].V[0] == j) {
					int index_N0 = m_pic.F[i].N[0];
					Vector3f temp(m_pic.VN[index_N0].NX, m_pic.VN[index_N0].NY, m_pic.VN[index_N0].NZ);
					sum = sum + temp;
					k++;
				}
				else if (m_pic.F[i].V[1] == j) {
					int index_N1 = m_pic.F[i].N[1];
					Vector3f temp(m_pic.VN[index_N1].NX, m_pic.VN[index_N1].NY, m_pic.VN[index_N1].NZ);
					sum = sum + temp;
					k++;
				}
				else if(m_pic.F[i].V[2] == j){
					int index_N2 = m_pic.F[i].N[2];
					Vector3f temp(m_pic.VN[index_N2].NX, m_pic.VN[index_N2].NY, m_pic.VN[index_N2].NZ);
					sum = sum + temp;
					k++;
				}
				/*
				if (m_pic.F[i].V[0] == j|| m_pic.F[i].V[1] == j || m_pic.F[i].V[2] == j) {
					int index_N0 = m_pic.F[i].N[0];
					//int index_N1 = m_pic.F[i].N[1];
					//int index_N2 = m_pic.F[i].N[2];
					Vector3f temp(m_pic.VN[index_N0].NX, m_pic.VN[index_N0].NY, m_pic.VN[index_N0].NZ);
					sum = sum + temp;
					k++;
				}*/				
			}
			sum(0) = sum(0) / k;
			sum(1) = sum(1) / k;
			sum(2) = sum(2) / k;
			normal.row(j) = sum.transpose();
		}		
	}
	if (m_pic.F.size() > 0) {
		for (int k = 0; k < 3; k++){//第几行
			for (int i = 0; i < numFaces; i++){
				/*
				vector <float>zuobiao;
				zuobiao.push_back(m_pic.V[m_pic.F[i].V[k]].X);
				zuobiao.push_back(m_pic.V[m_pic.F[i].V[k]].Y);
				zuobiao.push_back(m_pic.V[m_pic.F[i].V[k]].Z);*/
				Vector3f temp(3);
				temp(0) = m_pic.V[m_pic.F[i].V[k]].X;
				temp(1) = m_pic.V[m_pic.F[i].V[k]].Y;
				temp(2) = m_pic.V[m_pic.F[i].V[k]].Z;
				face[k].row(i) = temp.transpose();//赋值面矩阵 把索引改成顶点
			}
		}
	}
	//清除内存
	m_pic.V.swap(vector<Vertex>());
	m_pic.VN.swap(vector< FaXiangLiang>());
	//m_pic.VT.swap(vector< WenLi>());
	m_pic.F.swap(vector<Mian>());
	/*
	cout << "face[0]=" << face[0] << endl;
	cout << "face[1]=" << face[1] << endl;
	cout << "face[2]=" << face[2] << endl;
	cout << "normal=" << normal << endl;
	cout << "vertex=" << vertex << endl;*/
}
MatrixXf CKirchhoff::computeKF(float offset){
	MatrixXf KF(6, 6);
	MatrixXf C = face_center();
	//cout<< "C=" << C << endl;//对
	MatrixXf S = vertex - offset * normal;
	//cout << "S=" << S << endl;//对
	MatrixXf M = solid_angle(S);
	//cout << "M=" << M << endl;//对
	MatrixXf FL = motion_flux();
	//cout << " FL=" << FL << endl;//对
	MatrixXf sigma(numPoints, 6);
	for (int i = 0; i< 6;i++) {
		sigma.col(i) = M.fullPivHouseholderQr().solve(FL.col(i));//sigma=strength NAN!
	}
//	cout << "sigma=" << sigma << endl;
//	cout << "M*sigma=" << M*sigma << endl;
	MatrixXf SL = single_layer(S ,C);
//	cout << "SL=" << SL << endl;//对
	MatrixXf phi = numPoints *SL* sigma;//  SL * sigma;
//	cout << "phi=" << phi << endl;
	MatrixXf Q = one_point_quadrature();
//	cout << "Q=" << Q << endl;	
	KF = Q* phi;
	cout << "KF=" << KF << endl;
	return KF;
}
	
MatrixXf CKirchhoff::single_layer(MatrixXf S , MatrixXf C){
	MatrixXf res(numFaces, numPoints);
	for (int i = 0; i < numFaces; i++){//C Z
#pragma omp parallel
		for (int j = 0;j < numPoints;j++) {//S
			Vector3f temp = (C.row(i) - S.row(j)).transpose();
			float isZero = sqrt(temp.dot(temp));
			/*
			if (isZero<0.005 && isZero>-0.005)//防止除以0
				isZero = 0.005;
			*/
			res(i,j) = 1.0/isZero;
		}
	}
	return res;
}

MatrixXf CKirchhoff::face_normal(){
	MatrixXf AV = area_vector();
	MatrixXf areas = triangle_area();
	for (int i = 0; i < 3; i++){
		AV.col(i) = AV.col(i).cwiseQuotient(areas);// AV/w对应相除
	}
	return AV;
}

MatrixXf CKirchhoff::face_center(){
	MatrixXf res(numFaces, 3);
	res = (face[0] + face[1] + face[2])/3;
	return res;
}
float CKirchhoff::area(){
	MatrixXf temp = triangle_area();
	float res= temp.colwise().sum()(0,0);
	return res;
}
MatrixXf CKirchhoff::motion_flux(){
	MatrixXf res(numFaces, 6);
	res.leftCols(3) = angular_vector();
	res.rightCols(3) = area_vector();
	return res;
}

MatrixXf CKirchhoff::triangle_area(){
	MatrixXf res(numFaces, 1);
	MatrixXf p01 = face[1] - face[0];
	MatrixXf p02 = face[2] - face[0];
#pragma omp parallel
	for (int i = 0; i < numFaces; i++){
		Vector3f p01i = p01.row(i).transpose();
		Vector3f p02i = p02.row(i).transpose();
		Vector3f temp = p01i.cross(p02i);
		float length = sqrt(temp(0)*temp(0) + temp(1)*temp(1) + temp(2)*temp(2));
		//计算叉乘的模
		res(i, 0) = 0.5*length;
	}
	return res;
}
MatrixXf CKirchhoff::one_point_quadrature(){
	MatrixXf areas = triangle_area();
	MatrixXf VF = face_center();
	MatrixXf NF = face_normal();
	MatrixXf CR(numFaces, 3);
#pragma omp parallel
	for (int i = 0;i < numFaces;i++) {
		Vector3f tempVF = VF.row(i).transpose();
		Vector3f tempNF = NF.row(i).transpose();
		//cout << "tempVF=" << tempVF << endl;
		CR.row(i) = tempVF.cross(tempNF);
	}
	MatrixXf Q(6, numFaces);
	Q.row(0) = (areas.array()*CR.col(0).array()).matrix().transpose();
	Q.row(1) = (areas.array()*CR.col(1).array()).matrix().transpose();
	Q.row(2) = (areas.array()*CR.col(2).array()).matrix().transpose();
	Q.row(3) = (areas.array()*NF.col(0).array()).matrix().transpose();
	Q.row(4) = (areas.array()*NF.col(1).array()).matrix().transpose();
	Q.row(5) = (areas.array()*NF.col(2).array()).matrix().transpose();
	return Q;
}
MatrixXf CKirchhoff::angular_vector(){
	
	MatrixXf p1 = face[0];
	MatrixXf p2 = face[1];
	MatrixXf p3 = face[2];
	MatrixXf pp1(numFaces, 1);
	MatrixXf pp2(numFaces, 1);
	MatrixXf pp3(numFaces, 1);
#pragma omp parallel
	for (int i = 0; i < numFaces; i++){
		pp1(i, 0) = face[0].row(i).dot(face[0].row(i));
		pp2(i, 0) = face[1].row(i).dot(face[1].row(i));
		pp3(i, 0) = face[2].row(i).dot(face[2].row(i));
	}

	MatrixXf pp12(numFaces, 1);
	MatrixXf pp23(numFaces, 1);
	MatrixXf pp31(numFaces, 1);
#pragma omp parallel
	for (int i = 0; i < numFaces; i++){
		pp12(i, 0) = face[0].row(i).dot(face[1].row(i));
		pp23(i, 0) = face[2].row(i).dot(face[1].row(i));
		pp31(i, 0) = face[2].row(i).dot(face[0].row(i));
	}
	MatrixXf p12= pp1 + pp2 + pp12;
	MatrixXf p23= pp2 + pp3 + pp23;
	MatrixXf p31= pp3 + pp1 + pp31;

	MatrixXf p2p1 = p2 - p1;
	MatrixXf p3p2 = p3 - p2;
	MatrixXf p1p3 = p1 - p3;
	MatrixXf t1(numFaces, 3);
	MatrixXf t2(numFaces, 3);
	MatrixXf t3(numFaces, 3);
#pragma omp parallel
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
	MatrixXf res = -(t1 + t2 + t3) / 6;
	/*
	MatrixXf res(numFaces,3);
	for (int i = 0;i < numFaces;i++) {
		Vector3f omega1(1, 0, 0);
		Vector3f omega2(0, 1, 0);
		Vector3f omega3(0, 0, 1);
		Matrix3f temp;
		temp.row(0) = omega1.transpose();
		temp.row(1) = normal.row(i);
		temp.row(2) = face_center().row(i);
		res(i, 0) = -0.5*temp.determinant();
		temp.row(0) = omega2.transpose();
		res(i, 1) = -0.5*temp.determinant();
		temp.row(0) = omega3.transpose();
		res(i, 2) = -0.5*temp.determinant();
	}*/
	return res;
}
MatrixXf CKirchhoff::area_vector(){
	MatrixXf AV(numFaces,3);
	for (int i = 0;i < numFaces;i++) {
		Vector3f p1 = face[0].row(i);
		Vector3f p2 = face[1].row(i);
		Vector3f p3 = face[2].row(i);


		Vector3f temp = 0.5*(p1.cross(p2) + p2.cross(p3) + p3.cross(p1));
		AV.row(i) = temp.transpose();
	}
	/*
	MatrixXf area = triangle_area();//向量
	
	for (int i = 0;i < numFaces;i++) {
		AV(i, 0) = area(i, 0) * normal(i, 0);
		AV(i, 1) = area(i, 0) * normal(i, 1);
		AV(i, 2) = area(i, 0) * normal(i, 2);
	}*/

	return AV;
}
MatrixXf CKirchhoff::solid_angle(MatrixXf src){//numPoints*3
	MatrixXf res(numFaces, numPoints);
	cout << "S=" << src << endl;
//	MatrixXf FacePoint(numFaces, 3);
	for (int i = 0; i < numFaces; i++){
#pragma omp parallel
		for (int j = 0; j < numPoints; j++){
			Vector3f R1 = (face[0].row(i) - src.row(j)).transpose();
			Vector3f R2 = (face[1].row(i) - src.row(j)).transpose();
			Vector3f R3 = (face[2].row(i) - src.row(j)).transpose();
			MatrixXf temp(3,3);
			temp.row(0) = R1.transpose();
			temp.row(1) = R2.transpose();
			temp.row(2) = R3.transpose();
			float N = temp.determinant();
			//float test= R1.dot(R2.cross(R3));
			float l1 = sqrt(R1.dot(R1));
			float l2 = sqrt(R2.dot(R2));
			float l3 = sqrt(R3.dot(R3));
			float Den = l1*l2*l3 + l1*R2.dot(R3) + l2*R1.dot(R3) + l3*R1.dot(R2);
			res(i, j) = 2 * numPoints*atan2(N , Den);
		}
	}
	for (int i = 0;i < numPoints;i++) {
		float p=res.col(i).sum();
//		cout << "是不是4?" << p / 3.1415926 << endl;
	}
	return res;
}
void CKirchhoff::Subexpressions(float &w0, float &w1, float &w2, float &f1, float &f2, float &f3, float &g0, float &g1, float &g2) {
	float temp0 = w0 + w1;
	f1 = temp0 + w2;
	float temp1 = w0*w0;
	float temp2 = temp1 + w1*temp0;
	f2 = temp2 + w2*f1;
	f3 = w0*temp1 + w1*temp2 + w2*f2;
	g0 = f2 + w0*(f1 + w0);
	g1 = f2 + w1*(f1 + w1);
	g2 = f2 + w2*(f1 + w2);
}
Matrix3f CKirchhoff::comuputeJ() {
	float mass = 0;
	Matrix3f inertia;
	const float mult[10] = { 1.0 / 6, 1.0 / 24, 1.0 / 24, 1.0 / 24, 
		1.0 / 60, 1.0 / 60, 1.0 / 60, 1.0 / 120, 1.0 / 120, 1.0 / 120 };
	float intg[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; // order: 1, x, y, z, x^2, y^2, z^2, xy, yz, zx 

	float f1x = 0, f2x = 0, f3x = 0;
	float f1y = 0, f2y = 0, f3y = 0;
	float f1z = 0, f2z = 0, f3z = 0;
	float g0x = 0, g1x = 0, g2x = 0;
	float g0y = 0, g1y = 0, g2y = 0;
	float g0z = 0, g1z = 0, g2z = 0;
#pragma omp parallel
	for (int t = 0; t < numFaces; t++) { // get vertices of triangle t 
		/*
		i0 = index[3*t]; 	i1 = index[3*t+1];		i2 = index[3*t+2];
		x0 = p[i0].x; 		y0 = p[i0].y; 		z0 = p[i0].z;
		x1 = p[i1].x;       y1 = p[i1].y;       z1 = p[i1].z;
		x2 = p[i2].x;       y2 = p[i2].y;       z2 = p[i2].z;
		*/
		float x0 = face[0](t,0);
		float y0 = face[0](t,1);
		float z0 = face[0](t,2);

		float x1 = face[1](t,0);
		float y1 = face[1](t,1);
		float z1 = face[1](t,2);

		float x2 = face[2](t,0);
		float y2 = face[2](t,1);
		float z2 = face[2](t,2);
		// get edges and cross product of edges 
		float a1 = x1 - x0; float b1 = y1 - y0; float c1 = z1 - z0;
		float a2 = x2 - x0; float b2 = y2 - y0; float c2 = z2 - z0;
		float d0 = b1*c2 - b2*c1; float d1 = a2*c1 - a1*c2; float d2 = a1*b2 - a2*b1;
		// compute integral terms
		Subexpressions(x0, x1, x2, f1x, f2x, f3x, g0x, g1x, g2x);
		Subexpressions(y0, y1, y2, f1y, f2y, f3y, g0y, g1y, g2y);
		Subexpressions(z0, z1, z2, f1z, f2z, f3z, g0z, g1z, g2z);
		// update integrals 
		intg[0] += d0*f1x ;
		intg[1] += d0*f2x ; intg[2] += d1*f2y ; intg[3] += d2*f2z ;
		intg[4] += d0*f3x ; intg[5] += d1*f3y ; intg[6] += d2*f3z ;
		intg[7] += d0*(y0*g0x + y1*g1x + y2*g2x) ; 
		intg[8] += d1*(z0*g0y + z1*g1y + z2*g2y) ; 
		intg[9] += d2*(x0*g0z + x1*g1z + x2*g2z) ;
	}
	for (int i = 0;i < 10;i++) {
		intg[i] *= mult[i];
	}
	mass = intg[0];
	// 质心
	Vector3f cm(intg[1] / mass, intg[2] / mass, intg[3] / mass);
	// 相对于质心的惯性张量   0x 1y 2z
	inertia(0, 0) = intg[5] + intg[6] - mass*(cm(1)*cm(1) + cm(2)*cm(2));//yz 12
	inertia(1, 1) = intg[4] + intg[6] - mass*(cm(2)*cm(2) + cm(0)*cm(0));//zx 20
	inertia(2, 2) = intg[4] + intg[5] - mass*(cm(0)*cm(0) + cm(1)*cm(1));//xy 01
	inertia(0, 1) = inertia(1, 0) = -(intg[7] - mass*cm(0)*cm(1));
	inertia(1, 2) = inertia(2, 1) = -(intg[8] - mass*cm(1)*cm(2));
	inertia(0, 2) = inertia(2, 0) = -(intg[9] - mass*cm(2)*cm(0));
	return inertia;
}
MatrixXf CKirchhoff::computeKB(float m) {
	MatrixXf res(6, 6);
	res.setZero(6, 6);
	Matrix3f J = comuputeJ();
	res.block(0, 0, 3, 3) = J;
	Matrix3f identity;
	identity.setIdentity(3, 3);
	res.block(3, 3, 3, 3) = m*identity;//标量乘以矩阵
	return res;
}
MatrixXf CKirchhoff::computeK(){
	MatrixXf KB = computeKB(5.0f);//mass
	MatrixXf KF = computeKF(0.01);//offest 不同的模型需要修改 避免源点跑出去
	//MatrixXf temp = KF.rowwise().sum();
	/*
	float a = KB[0]; //temp(0, 0);
	float b = KB[1];//temp(1, 0);
	float c = KB[2];//temp(2, 0);
	float d = KB[3];//temp(3, 0);
	float e = KB[4];//temp(4, 0);
	float f = KB[5];//temp(5, 0);*/
	//KB.setIdentity();//先用单位阵试试
	/*
	KF.setZero();
	KF(0, 3)= KF(3, 0) = 0.9565;
	KF(1, 4) = KF(4, 1) = -1.9204;
	KF(2, 5) = KF(5, 2) = 0.9565;

	KF(0, 0) = 2.6242;
	KF(1, 1) = 3.5184;
	KF(2, 2) = 2.6239;
	KF(3, 3) = 1.8576;
	KF(4, 4) = 4.7054;
	KF(5, 5) = 1.8585;
	cout << "KB:" << KB << endl;*/
	MatrixXf Kirchhoff = KB + KF;
	cout << "KB:" << KB << endl;
	cout << "Kirchhoff:" << Kirchhoff << endl;
	return Kirchhoff;
}