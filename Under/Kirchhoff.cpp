#include "Kirchhoff.h"
#include<iostream>
using namespace Eigen;
using namespace std;
CKirchhoff::CKirchhoff(PICnew *m_picnew, double m_bodyDensity, double m_fluidDensity)
{
	numPoints = m_picnew->vertexandnormal.size();
		// m_pic.V.size();
	numFaces = m_picnew->faceandnormal.size();
		//m_pic.F.size();
	vertex.resize(numPoints, 3);
	normal.resize(numPoints, 3);
	face[0].resize(numFaces, 3);
	face[1].resize(numFaces, 3);
	face[2].resize(numFaces, 3);
	if (numPoints> 0) {
		for (int i = 0; i < numPoints; i++) {
			Vector3d temp(3);
			temp(0) = m_picnew->vertexandnormal[i].coordinate(0);
				//m_pic.V[i].X;
			temp(1) = m_picnew->vertexandnormal[i].coordinate(1);// m_pic.V[i].Y;
			temp(2) = m_picnew->vertexandnormal[i].coordinate(2);// m_pic.V[i].Z;		
			vertex.row(i) = temp.transpose();//赋值顶点矩阵
		}
	}
	if (numFaces> 0) {
		/*
		for (int j = 0;j < numPoints; j++) {
			Vector3d sum(0,0,0);
			int k = 0;
#pragma omp parallel for
			for (int i = 0; i < numFaces; i++){
				if (m_pic.F[i].V[0] == j) {
					int index_N0 = m_pic.F[i].N[0];
					Vector3d temp(m_pic.VN[index_N0].NX, m_pic.VN[index_N0].NY, m_pic.VN[index_N0].NZ);
					sum = sum + temp;
					k++;
				}
				else if (m_pic.F[i].V[1] == j) {
					int index_N1 = m_pic.F[i].N[1];
					Vector3d temp(m_pic.VN[index_N1].NX, m_pic.VN[index_N1].NY, m_pic.VN[index_N1].NZ);
					sum = sum + temp;
					k++;
				}
				else if(m_pic.F[i].V[2] == j){
					int index_N2 = m_pic.F[i].N[2];
					Vector3d temp(m_pic.VN[index_N2].NX, m_pic.VN[index_N2].NY, m_pic.VN[index_N2].NZ);
					sum = sum + temp;
					k++;
				}
				
				//if (m_pic.F[i].V[0] == j|| m_pic.F[i].V[1] == j || m_pic.F[i].V[2] == j) {
				//	int index_N0 = m_pic.F[i].N[0];
				//	//int index_N1 = m_pic.F[i].N[1];
				//	//int index_N2 = m_pic.F[i].N[2];
				//	Vector3d temp(m_pic.VN[index_N0].NX, m_pic.VN[index_N0].NY, m_pic.VN[index_N0].NZ);
				//	sum = sum + temp;
				//	k++;
				//}			
			}
			sum(0) = sum(0) / k;
			sum(1) = sum(1) / k;
			sum(2) = sum(2) / k;
			normal.row(j) = sum.transpose();
		}		*/
		for (int i = 0;i < numPoints;i++) {
			Vector3d temp(3);
			temp(0) = m_picnew->vertexandnormal[i].vertexNormal(0);
			//m_pic.VN[i].NX;
			temp(1) = m_picnew->vertexandnormal[i].vertexNormal(1);
			temp(2) = m_picnew->vertexandnormal[i].vertexNormal(2);
			//Vector3d temp(m_pic.VN[i].NX, m_pic.VN[i].NY, m_pic.VN[i].NZ);
			normal.row(i) = temp.transpose();
		}
	}
	if (numFaces> 0) {
		for (int k = 0; k < 3; k++){//第几行
			for (int i = 0; i < numFaces; i++){
				/*
				vector <double>zuobiao;
				zuobiao.push_back(m_pic.V[m_pic.F[i].V[k]].X);
				zuobiao.push_back(m_pic.V[m_pic.F[i].V[k]].Y);
				zuobiao.push_back(m_pic.V[m_pic.F[i].V[k]].Z);*/
				Vector3d temp(3);
				int a=m_picnew->faceandnormal[i].vertexIndex[k];//顶点索引
				temp(0) = m_picnew->vertexandnormal[a].coordinate(0);
					//m_pic.V[m_pic.F[i].V[k]].X;
				temp(1) = m_picnew->vertexandnormal[a].coordinate(1);
				temp(2) = m_picnew->vertexandnormal[a].coordinate(2);
				face[k].row(i) = temp.transpose();//赋值面矩阵 把索引改成顶点
			}
		}
	}
	//清除内存
	//m_pic.V.swap(vector<Vertex>());
	//m_pic.VN.swap(vector< FaXiangLiang>());
	////m_pic.VT.swap(vector< WenLi>());
	//m_pic.F.swap(vector<Mian>());
	fluidDensity = m_fluidDensity;
	bodyDensity = m_bodyDensity;
	
	/*cout << "face[0]=" << face[0] << endl;
	cout << "face[1]=" << face[1] << endl;
	cout << "face[2]=" << face[2] << endl;
	cout << "normal=" << normal << endl;
	cout << "vertex=" << vertex << endl;*/
}
MatrixXd CKirchhoff::computeKF(double offset){
	MatrixXd KF(6, 6);
	MatrixXd C = face_center();
	//cout<< "C=" << C << endl;//对
	MatrixXd S = vertex - offset * normal;
	//cout << "S=" << S << endl;//对
	MatrixXd M = solid_angle(S);
	//cout << "M=" << M << endl;//对
	MatrixXd FL = motion_flux();
	//cout << " FL=" << FL << endl;//对
	MatrixXd sigma(numPoints, 6);
	for (int i = 0; i< 6;i++) {
		sigma.col(i) = M.fullPivHouseholderQr().solve(FL.col(i));//sigma=strength NAN!
	}
	//cout << "sigma=" << sigma << endl;
	//cout << "M*sigma=" << M*sigma << endl;
	MatrixXd SL = single_layer(S ,C);
	//cout << "SL=" << SL << endl;//对
	MatrixXd phi = numPoints *SL* sigma;//  SL * sigma;
	//cout << "phi=" << phi << endl;
	MatrixXd Q = one_point_quadrature();
	//cout << "Q=" << Q << endl;	
	KF = Q* phi;
	return KF;
}
	
MatrixXd CKirchhoff::single_layer(MatrixXd S , MatrixXd C){
	MatrixXd res(numFaces, numPoints);
	for (int i = 0; i < numFaces; i++){//C Z
		for (int j = 0;j < numPoints;j++) {//S
			Vector3d temp = (C.row(i) - S.row(j)).transpose();
			double isZero = sqrt(temp.dot(temp));
			/*
			if (isZero<0.005 && isZero>-0.005)//防止除以0
				isZero = 0.005;
			*/
			res(i,j) = 1.0/isZero;
		}
	}
	return res;
}

MatrixXd CKirchhoff::face_normal(){
	MatrixXd AV = area_vector();
	MatrixXd areas = triangle_area();
	for (int i = 0; i < 3; i++){
		AV.col(i) = AV.col(i).cwiseQuotient(areas);// AV/w对应相除
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
	MatrixXd p01 = face[1] - face[0];
	MatrixXd p02 = face[2] - face[0];
	for (int i = 0; i < numFaces; i++){
		Vector3d p01i = p01.row(i).transpose();
		Vector3d p02i = p02.row(i).transpose();
		Vector3d temp = p01i.cross(p02i);
		double length = sqrt(temp(0)*temp(0) + temp(1)*temp(1) + temp(2)*temp(2));
		//计算叉乘的模
		res(i, 0) = 0.5*length;
	}
	return res;
}
MatrixXd CKirchhoff::one_point_quadrature(){
	MatrixXd areas = triangle_area();
	MatrixXd VF = face_center();
	MatrixXd NF = face_normal();
	MatrixXd CR(numFaces, 3);
	for (int i = 0;i < numFaces;i++) {
		Vector3d tempVF = VF.row(i).transpose();
		Vector3d tempNF = NF.row(i).transpose();
		//cout << "tempVF=" << tempVF << endl;
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
	MatrixXd res = -(t1 + t2 + t3) / 6;
	/*
	MatrixXd res(numFaces,3);
	for (int i = 0;i < numFaces;i++) {
		Vector3d omega1(1, 0, 0);
		Vector3d omega2(0, 1, 0);
		Vector3d omega3(0, 0, 1);
		Matrix3d temp;
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
MatrixXd CKirchhoff::area_vector(){
	MatrixXd AV(numFaces,3);
	for (int i = 0;i < numFaces;i++) {
		Vector3d p1 = face[0].row(i);
		Vector3d p2 = face[1].row(i);
		Vector3d p3 = face[2].row(i);

		Vector3d temp = 0.5*(p1.cross(p2) + p2.cross(p3) + p3.cross(p1));
		AV.row(i) = temp.transpose();
	}
	/*
	MatrixXd area = triangle_area();//向量
	
	for (int i = 0;i < numFaces;i++) {
		AV(i, 0) = area(i, 0) * normal(i, 0);
		AV(i, 1) = area(i, 0) * normal(i, 1);
		AV(i, 2) = area(i, 0) * normal(i, 2);
	}*/
	return AV;
}
MatrixXd CKirchhoff::solid_angle(MatrixXd src){//numPoints*3
	MatrixXd res(numFaces, numPoints);
//	cout << "S=" << src << endl;
//	MatrixXd FacePoint(numFaces, 3);
	for (int i = 0; i < numFaces; i++){
		for (int j = 0; j < numPoints; j++){
			Vector3d R1 = (face[0].row(i) - src.row(j)).transpose();
			Vector3d R2 = (face[1].row(i) - src.row(j)).transpose();
			Vector3d R3 = (face[2].row(i) - src.row(j)).transpose();
			MatrixXd temp(3,3);
			temp.row(0) = R1.transpose();
			temp.row(1) = R2.transpose();
			temp.row(2) = R3.transpose();
			double N = temp.determinant();
			//double test= R1.dot(R2.cross(R3));
			double l1 = sqrt(R1.dot(R1));
			double l2 = sqrt(R2.dot(R2));
			double l3 = sqrt(R3.dot(R3));
			double Den = l1*l2*l3 + l1*R2.dot(R3) + l2*R1.dot(R3) + l3*R1.dot(R2);
			res(i, j) = 2 *numPoints*atan2(N , Den);
		}
	}
	for (int i = 0;i < numPoints;i++) {
		double p=res.col(i).sum();
//		cout << "是不是4?" << p / 3.1415926 << endl;
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
Matrix3d CKirchhoff::computeJ() {//之前单词打错了
	Matrix3d inertia;
	const double mult[10] = { 1.0 / 6, 1.0 / 24, 1.0 / 24, 1.0 / 24, 
		1.0 / 60, 1.0 / 60, 1.0 / 60, 1.0 / 120, 1.0 / 120, 1.0 / 120 };
	double intg[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; // order: 1, x, y, z, x^2, y^2, z^2, xy, yz, zx 

	double f1x = 0, f2x = 0, f3x = 0;
	double f1y = 0, f2y = 0, f3y = 0;
	double f1z = 0, f2z = 0, f3z = 0;
	double g0x = 0, g1x = 0, g2x = 0;
	double g0y = 0, g1y = 0, g2y = 0;
	double g0z = 0, g1z = 0, g2z = 0;
	for (int t = 0; t < numFaces; t++) { // get vertices of triangle t 
		/*
		i0 = index[3*t]; 	i1 = index[3*t+1];		i2 = index[3*t+2];
		x0 = p[i0].x; 		y0 = p[i0].y; 		z0 = p[i0].z;
		x1 = p[i1].x;       y1 = p[i1].y;       z1 = p[i1].z;
		x2 = p[i2].x;       y2 = p[i2].y;       z2 = p[i2].z;
		*/
		double x0 = face[0](t,0);//该面第一个索引点的坐标
		double y0 = face[0](t,1);
		double z0 = face[0](t,2);

		double x1 = face[1](t,0);
		double y1 = face[1](t,1);
		double z1 = face[1](t,2);

		double x2 = face[2](t,0);
		double y2 = face[2](t,1);
		double z2 = face[2](t,2);
		// get edges and cross product of edges 
		double a1 = x1 - x0; double b1 = y1 - y0; double c1 = z1 - z0;
		double a2 = x2 - x0; double b2 = y2 - y0; double c2 = z2 - z0;
		double d0 = b1*c2 - b2*c1; double d1 = a2*c1 - a1*c2; double d2 = a1*b2 - a2*b1;
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
	volume = intg[0];
	cout << "volume" << volume << endl;
	fluidMass = intg[0]*fluidDensity;
	for (int i = 0;i < 10;i++) {
		intg[i] *=  bodyDensity;
	}
	bodyMass = intg[0];
	// 质心
	Vector3d cm(intg[1] / bodyMass, intg[2] / bodyMass, intg[3] / bodyMass);
	Cm = cm;
	// 相对于质心的惯性张量   0x 1y 2z
	inertia(0, 0) = intg[5] + intg[6] - bodyMass *(cm(1)*cm(1) + cm(2)*cm(2));//yz 12
	inertia(1, 1) = intg[4] + intg[6] - bodyMass *(cm(2)*cm(2) + cm(0)*cm(0));//zx 20
	inertia(2, 2) = intg[4] + intg[5] - bodyMass *(cm(0)*cm(0) + cm(1)*cm(1));//xy 01
	inertia(0, 1) = inertia(1, 0) = -(intg[7] - bodyMass *cm(0)*cm(1));
	inertia(1, 2) = inertia(2, 1) = -(intg[8] - bodyMass *cm(1)*cm(2));
	inertia(0, 2) = inertia(2, 0) = -(intg[9] - bodyMass *cm(2)*cm(0));
	return inertia;
}
MatrixXd CKirchhoff::computeKB() {
	MatrixXd res(6, 6);
	res.setZero(6, 6);
	Matrix3d J = computeJ();//计算结束得到mass，注意并行问题
	res.block(0, 0, 3, 3) = J;
	Matrix3d identity;
	identity.setIdentity(3, 3);
	res.block(3, 3, 3, 3) = bodyMass *identity;//标量乘以矩阵
	return res;
}
MatrixXd CKirchhoff::computeK(){
	MatrixXd KB = computeKB();
	MatrixXd KF = computeKF(0.02);//offest 不同的模型需要修改 避免源点跑出去

	cout << "KF=" << KF << endl;
	//MatrixXd temp = KF.rowwise().sum();
	/*
	double a = KB[0]; //temp(0, 0);
	double b = KB[1];//temp(1, 0);
	double c = KB[2];//temp(2, 0);
	double d = KB[3];//temp(3, 0);
	double e = KB[4];//temp(4, 0);
	double f = KB[5];//temp(5, 0);*/
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
	//KB.setIdentity();
	MatrixXd Kirchhoff = KB + KF;//
	cout << "KB:" << KB << endl;
	cout << "Kirchhoff:" << Kirchhoff << endl;
	return Kirchhoff;
}
//不用在这里计算空间坐标的合外力和合外力矩
VectorXd CKirchhoff::computetsfs() {
	VectorXd tsfs(6);
	Vector3d g(0, -9.8, 0);
	Vector3d tg(0, 0, 0);//质量均匀分布
	Vector3d fg(0, 0, 0);//=(bodyMass - volume*fluidDensity)*g;// 
	tsfs.block(0, 0, 3, 1) = tg;
	tsfs.block(3, 0, 3, 1) = fg;
	return tsfs;
}
