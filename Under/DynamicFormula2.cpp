#include "DynamicFormula2.h"
#include <cmath>
double DynamicFormula2::mu=0.001;//

//计算每个面上的应力
VectorXd DynamicFormula2::computetraction(vector<Body*> m_body){
	int n = m_body[0]->idnum;//获取总面片数,idnum是静态变量
	MatrixXd coefficient(3*n, 3*n);
	coefficient.setZero();
	VectorXd traction(3 * n);
	MatrixXd H(3*n, 3*n);
	H.setZero();
	VectorXd u(3*n);		
	Matrix3d identity;
	identity.setIdentity();
	for (int p = 0;p < m_body.size();p++) {
		for (int q = 0;q < m_body[p]->v_onepoint.size();q++) {//field point
			Matrix3d R = m_body[p]->g.block(0, 0, 3, 3);
			Vector3d y_ = m_body[p]->g.block(0, 3, 3, 1);//物体在世界坐标的位置
			Vector3d y = m_body[p]->v_onepoint[q].midpoint;//物体坐标系下field点的位置
			//cout << "看看y" << y << endl; 不需要y，面的中点了
			// 得到y， q面中点在世界坐标系的位置
			int b= m_body[p]->v_onepoint[q].id;//下标符号
			Vector3d ny = R*m_body[p]->v_onepoint[q].normal;//得到y的法向
		/*	Vector3d np0 = R * m_body[p]->v_onepoint[q].normal3[0];
			Vector3d np1 = R * m_body[p]->v_onepoint[q].normal3[1];
			Vector3d np2 = R * m_body[p]->v_onepoint[q].normal3[2];*/
			Vector3d w = m_body[p]->epsilon.block(0, 0, 3, 1);
			Vector3d v = m_body[p]->epsilon.block(3, 0, 3, 1);
			u.block(3 * b,0,3,1) = R*(w.cross(y) + v);
			//三个顶点
			Vector3d p0 = R * m_body[p]->v_onepoint[q].vertex[0] + y_;
			Vector3d p1 = R * m_body[p]->v_onepoint[q].vertex[1] + y_;
			Vector3d p2 = R * m_body[p]->v_onepoint[q].vertex[2] + y_;
			//三个顶点的中点
			Vector3d p01 = (p0 + p1) / 2;
			Vector3d p02 = (p0 + p2) / 2;
			Vector3d p12 = (p1 + p2) / 2;
			//S0
			Vector3d mid0 = (p0 + p01 + p02) / 3;
			double S0 = trianglearea(p0, p01, p02);
			//S1
			Vector3d mid1 = (p1 + p01 + p12) / 3;
			double S1 = trianglearea(p1, p01, p12);
			//S2
			Vector3d mid2 = (p2 + p02 + p12) / 3;
			double S2 = trianglearea(p2, p02, p12);
			for (int i = 0;i < m_body.size();i++) {//source point
				for (int j = 0;j < m_body[i]->v_onepoint.size();j++) {
					Matrix3d R2 = m_body[i]->g.block(0, 0, 3, 3);
					Vector3d y_2 = m_body[i]->g.block(0, 3, 3, 1);//物体在 世界坐标的位置
					Vector3d x = R2*m_body[i]->v_onepoint[j].midpoint+y_2;
					//cout << "看看x" << x << endl;
					//得到sourcepoint x
					int a = m_body[i]->v_onepoint[j].id;//下标符号
					Vector3d nx = R*m_body[i]->v_onepoint[j].normal;//得到x的法向
					if (a==b) {//对角线
						int xifen_num = 7;
						double area = trianglearea(p0, p1, p2)/pow(4,xifen_num);
						Matrix3d KS1 = digui(p0, p1, p2, xifen_num, area, 0, x, nx, ny);
						//cout << "对角线的KS是" << KS1 << endl;
						coefficient.block(3 * a, 3 * b, 3, 3) = KS1 - 0.5*identity;	
						/*Matrix3d HS1= digui(p0, p1, p2, xifen_num, area, 1, x, nx, ny);
						H.block(3 * a, 3 * b, 3, 3) = HS1;*/
						/*
						//中间三角形继续划分
						//顶点中点的中点
						Vector3d p0102 = (p01 + p02) / 2;
						Vector3d p0112 = (p01 + p12) / 2;
						Vector3d p0212 = (p02 + p12) / 2;
						//S01
						Vector3d mid01 = (p01 + p0102 + p0112) / 3;
						double S01 = trianglearea(p01, p0102, p0112);
						//S02
						Vector3d mid02 = (p02 + p0102 + p0212) / 3;
						double S02 = trianglearea(p02, p0102, p0212);
						//S12
						Vector3d mid12 = (p12 + p0112 + p0212) / 3;
						double S12 = trianglearea(p12, p0112, p0212);
						//中间三角形继续划分
						//顶点中点的中点
						Vector3d p01020112 = (p0102 + p0112) / 2;
						Vector3d p01020212 = (p0102 + p0212) / 2;
						Vector3d p01120212 = (p0112 + p0212) / 2;
						//S0102
						Vector3d mid0102 = (p0102 + p01020112 + p01020212) / 3;
						double S0102 = trianglearea(p0102, p01020112, p01020212);
						//S0112
						Vector3d mid0112 = (p0112 + p01020112 + p01120212) / 3;
						double S0112 = trianglearea(p0112, p01020112, p01120212);
						//S0212
						Vector3d mid0212 = (p0212 + p01020212 + p01120212) / 3;
						double S0212 = trianglearea(p0212, p01020212, p01120212);
						Matrix3d KS = computeKij(x, nx, mid0, ny)*S0
							+ computeKij(x, nx, mid1, ny)*S1
							+ computeKij(x, nx, mid2, ny)*S2
							+ computeKij(x, nx, mid01, ny)*S01
							+ computeKij(x, nx, mid02, ny)*S02
							+ computeKij(x, nx, mid12, ny)*S12
							+ computeKij(x, nx, mid0102, ny)*S0102
							+ computeKij(x, nx, mid0112, ny)*S0112
							+ computeKij(x, nx, mid0212, ny)*S0212;
						cout << "KS" << endl <<KS << endl;						
						coefficient.block(3 * a, 3 * b, 3, 3) = KS - 0.5*identity;
						//Matrix3d HS = computeHij(x, nx, mid0, ny)*S0
						//	+ computeHij(x, nx, mid1, ny)*S1
						//	+ computeHij(x, nx, mid2, ny)*S2
						//	+ computeHij(x, nx, mid01, ny)*S01
						//	+ computeHij(x, nx, mid02, ny)*S02
						//	+ computeHij(x, nx, mid12, ny)*S12
						//	+ computeHij(x, nx, mid0102, ny)*S0102
						//	+ computeHij(x, nx, mid0112, ny)*S0112
						//	+ computeHij(x, nx, mid0212, ny)*S0212;
						////cout << "H" << endl << computeHij(x, nx, y, ny) << endl;
						//H.block(3 * a, 3 * b, 3, 3) = HS;
						*/
					}
					else{//非对角线
						int xifen_num2 = 3;
						double area = trianglearea(p0, p1, p2) / pow(4, xifen_num2);
						Matrix3d KS1 = digui(p0, p1, p2, xifen_num2, area, 0, x, nx, ny);
						//cout << "KS1!" << KS1 << endl;
						coefficient.block(3 * a, 3 * b, 3, 3) = KS1;
						Matrix3d HS1 = digui(p0, p1, p2, xifen_num2, area, 1, x, nx, ny);
						H.block(3 * a, 3 * b, 3, 3) = HS1;
						/*
						double weight[] = { 0.308641975308642,0.493827160493828,0.308641975308642
						,0.493827160493828 ,0.308641975308642 ,0.493827160493828
						,0.308641975308642 ,0.493827160493828 ,0.790123456790124 };
						//和我手写的图的顺序不太一样
						double point[] = { -0.774596669241483 ,-0.774596669241483
							,-0.774596669241483 	,0
							,-0.774596669241483 ,0.774596669241483
							,0 ,0.774596669241483
							,0.774596669241483 ,0.774596669241483
							,0.774596669241483 ,0
							,0.774596669241483 ,-0.774596669241483
							,0 ,-0.774596669241483
							,0 ,0
						};
						Matrix3d KS1;
						KS1.setZero();
						Matrix3d HS1;
						HS1.setZero();
						for (int i = 0;i < 9;i++) {
							KS1 += weight[i] * computef(p0, p1, p2
								, point[2 * i], point[2 * i + 1], x, nx, ny, 0);
							HS1 += weight[i] * computef(p0, p1, p2
								, point[2 * i], point[2 * i + 1], x, nx, ny, 1);
						}
						cout << "KS1" << KS1 << endl;
						coefficient.block(3 * a, 3 * b, 3, 3) = KS1;
						H.block(3 * a, 3 * b, 3, 3) = HS1;*/
						
				/*		//S3
						Vector3d mid3 = (p01 + p02 + p12) / 3;
						double S3 = trianglearea(p01, p02, p12);
						Matrix3d KS = computeKij(x, nx, mid0, ny)*S0
							+ computeKij(x, nx, mid1, ny)*S1
							+ computeKij(x, nx, mid2, ny)*S2
							+ computeKij(x, nx, mid3, ny)*S3;
						coefficient.block(3 * a, 3 * b, 3, 3) = KS;
						Matrix3d HS = computeHij(x, nx, mid0, ny)*S0
							+ computeHij(x, nx, mid1, ny)*S1
							+ computeHij(x, nx, mid2, ny)*S2
							+ computeHij(x, nx, mid3, ny)*S3;				
						cout << "KS" << KS << endl;
							//cout << "H"<<endl << computeHij(x, nx, y, ny) << endl;
						H.block(3 * a, 3 * b, 3, 3) = HS;*/
					}				
				}
			}
		}
	}
	Matrix3d zero;
	zero.setZero();
	//Matrix3d HSum;
	//HSum.setZero();
	//for (int p = 0;p < m_body.size();p++) {//遍历所有模型
	//	for (int i = 0;i < m_body[p]->v_onepoint.size();i++) {//模型的局部行
	//		int row = m_body[p]->v_onepoint[i].id;//行数
	//		HSum.setZero();
	//		for (int j = 0;j < m_body[p]->v_onepoint.size();j++) {//模型的局部列
	//			int col= m_body[p]->v_onepoint[j].id;//列数
	//			HSum += H.block(row * 3, col * 3, 3, 3);
	//		}
	//		H.block(row * 3, row * 3, 3, 3) = zero - HSum;
	//	}
	//}

	for (int i = 0;i < n;i++) {
		Matrix3d HSum;
		HSum.setZero();
		for (int j = 0;j < n&&i!=j;j++) {
			HSum += H.block(3 * i, 3 * j, 3, 3);
		}
		H.block(3 * i, 3 * i, 3, 3) = zero-HSum;//
	}
	//看一下对角线矩阵,前1/3看看
	cout << "系数矩阵" << endl;
	for (int i = 0;i < n/3;i++) {
		cout << coefficient.block(3 * i, 3 * i, 3, 3) << endl;
	}	
	cout << "H矩阵" << endl;
	for (int i = 0;i < n / 3;i++) {
		cout <<  H.block(3 * i, 3 * i, 3, 3) << endl;
	}
	//cout <<"coefficient"<<endl<< coefficient << endl;
	//解线性方程组
	VectorXd b = H * u;
	//cout << "H:" << endl << H << endl;
	//cout << "u:" << endl << u << endl;
	//cout << "b:" << endl << b << endl;
	traction = coefficient.fullPivHouseholderQr().solve(b);//sigma=strength NAN!
	int k = 0;//taction索引
	//计算整个面上的力
	for (int p = 0;p < m_body.size();p++) {
		for (int q = 0;q < m_body[p]->v_onepoint.size();q++) {
			//遍历所有的物体所有面
			double S = m_body[p]->v_onepoint[q].area;
			//cout << q << "个三角形的面积是" << S << endl;
			traction(k) *= S;
			traction(k+1) *= S;
			traction(k+2) *= S;
			k += 3;
		}
	}
	//cout << "乘回去的结果b:" << endl << coefficient*traction << endl;
	//cout <<"traction:"<< traction << endl;
	return traction;
}

Matrix3d DynamicFormula2::computeKij(Vector3d x, Vector3d nx, Vector3d y, Vector3d ny) {
	Matrix3d res;
	res.setZero();
	double r2 = pow(x(0) - y(0), 2) +
		pow(x(1) - y(1), 2) + pow(x(2) - y(2), 2);
	double r =sqrt(r2);//r的平方
	for (int i = 0;i < 3;i++) {
		for (int j = 0;j < 3;j++) {
			double rknk = 0;
			for (int k = 0;k < 3;k++) {
				//cout << "yk-xk" << y(k) - x(k) << endl;
				//cout << "nx(k)" << nx(k) << endl;
				rknk += (y(k)-x(k) ) / r * nx(k);
			}
			res(i, j) =( 3 / (4 * 3.1415926*r2))*( (y(i) - x(i) )*( y(j) - x(j)) )/ r2 *rknk;
			//r*r合并增加精度
			//cout << "i:" << i << "j:" << j <<":"<< res(i, j)<<endl;
		}
	}
	return res;
}
double DynamicFormula2::dirac(int a, int b) {
	if (a == b) {
		return 1;
	}
	else {
		return 0;
	}
}

Matrix3d DynamicFormula2::computeHij(Vector3d x, Vector3d nx, Vector3d y, Vector3d ny) {
	Matrix3d res;
	res.setZero();
	double r2 = pow(x(0) - y(0), 2) +
		pow(x(1) - y(1), 2) + pow(x(2) - y(2), 2);
	double r = sqrt(r2);//r的平方
	//Vector3d n = (y-x) / r;//y到x的单位向量
	double rlnl=0;
	for (int l = 0;l < 3;l++) {
		rlnl += (y(l) - x(l)) * ny(l);
	}
	rlnl /= r;//提高精度
	for (int i = 0;i < 3;i++) {
		for (int j = 0;j < 3 ;j++) {
			for (int k = 0;k < 3;k++) {
				res(i, j) += mu / (4 * 3.1415926)*
				(
					(3 * (dirac(i, j)* (y(k) - x(k)) / (r2 * r2) + dirac(j, k) * (y(i) - x(i)) / (r2 * r2) )
						- 30 * (y(i) - x(i))  * (y(j) - x(j)) * (y(k) - x(k)) / (r2*r2*r2)
					)  *rlnl
						+ 3 * ( ny(i)*(y(j) - x(j))  * (y(k) - x(k)) / (r2*r2*r) 
										+ ny(k) *(y(i) - x(i)) * (y(j) - x(j) ) / ( r2*r2*r) 
								 )
						+ 2 * dirac(i, k)*ny(j)/(r2*r)
				)*nx(k);
				//cout << "H i =" << i << "j=" << j << H(i, j) << endl;
			}
		}
	}
	return res;
}
Vector3d DynamicFormula2::chazhi(Vector3d p1, Vector3d p2, Vector3d p3, double yipu, double yita) {
	double N1 = (1 - yipu)*(1 - yita) / 4;
	double N2 = (1 + yipu)*(1 - yita) / 4;
	double N3 = (1 + yita) / 2;
	Vector3d res;
	res = N1 * p1 + N2 * p2 + N3 * p3;
	return res;
}
Matrix3d DynamicFormula2::computef(Vector3d p1, Vector3d p2, Vector3d p3
	, double yipu, double yita, Vector3d x, Vector3d nx
	, Vector3d ny, int flag) {
	//cout <<"p1!!" <<p1 << endl << p2 << endl << p3 << endl << yipu << endl << yita;
	//cout << "x" << x << endl << nx << endl << ny << endl;
	double piany1yipu =  (-1 + yita)*p1(0)/ 4 +  (1 - yita) *p2(0)/ 4; 
	double piany2yipu =  (-1 + yita)*p1(1)/ 4 +  (1 - yita) *p2(1)/ 4;
	double piany3yipu =  (-1 + yita)*p1(2)/ 4 +  (1 - yita) *p2(2)/ 4; 

	double piany1yita =  (-1 + yipu)*p1(0)/ 4 +  (-1 - yipu) *p2(0) / 4 + p3(0)/ 2;
	double piany2yita =  (-1 + yipu)*p1(1)/ 4  + (-1 - yipu) *p2(1) / 4 + p3(1)/ 2 ;
	double piany3yita =  (-1 + yipu)*p1(2)/ 4 +  (-1 - yipu) *p2(2) / 4 + p3(2)/ 2;
	MatrixXd jacobiMatrix(2,3);
	jacobiMatrix.setZero();
	Vector3d a;
	Vector3d b;
	a(0) = piany1yipu;	a(1) = piany2yipu;	a(2) = piany3yipu;
	b(0) = piany1yita;	b(1) = piany2yita;	b(2) = piany3yita;
	double jacobi = a.cross(b).norm();
	//jacobiMatrix(0, 0) = piany1yipu;	jacobiMatrix(0, 1) = piany2yipu;	jacobiMatrix(0, 2) = piany3yipu;
	//jacobiMatrix(1, 0) = piany1yita;	jacobiMatrix(1, 1) = piany2yita;	jacobiMatrix(1, 2) = piany3yita;
	////cout << "H矩阵jacobiMatrix" << endl<<jacobiMatrix << endl;
	//Matrix2d jMatrix = jacobiMatrix*jacobiMatrix.transpose();
	//double jacobi2 = jMatrix.determinant();
	////Vector3d a(piany1yipu, piany2yipu, piany3yipu);
	////Vector3d b(piany2yita, piany2yita, piany3yita);
	//double jacobi = sqrt(jacobi2);
	//cout << "jacobi" << jacobi << endl;
	Vector3d y = chazhi(p1, p2, p3, yipu, yita);
//	Vector3d ny = chazhi(np1, np2, np3, yipu, yita);
	Matrix3d res;
	res.setZero();
	if (flag == 0) {
		res = computeKij(x, y, nx, ny)*jacobi;
	}
	else if (flag == 1) {
		res = computeHij(x, y, nx, ny)*jacobi;
	}
	//cout << "f函数值是：" << endl << res << endl;
	return res;
}
Matrix3d DynamicFormula2::digui(Vector3d p0, Vector3d p1, Vector3d p2,
	int n, double area,int flag,Vector3d x,Vector3d nx,Vector3d ny) {//area是单个小三角形的面积
	if (n == 0) {
		Vector3d mid = (p0 + p1 + p2) / 3;
		if (flag == 0) {
			if (fabs(mid(0) - x(0)) < yipusilon && fabs(mid(1) - x(1)) < yipusilon &&fabs(mid(2) - x(2)) < yipusilon) {
				Matrix3d zero;
				zero.setZero();
				/*cout << "zero" << endl;
				cout << endl;
				cout << endl;
				cout << endl;*/
				return zero;
			}
			else {
				//Matrix3d temp = area * computeKij(x, nx, mid, ny);
				/*
				if (temp(0, 0) > 10000 || temp(0, 1) > 10000 || temp(0, 2) > 10000
					|| temp(1, 0) > 10000 || temp(1, 1) > 10000 || temp(1, 2) > 10000
					|| temp(2, 0) > 10000 || temp(2, 1) > 10000 || temp(2, 2) > 10000) {
					cout << "mid" << mid<< endl<<"x"<<x<<endl;
					if (mid == x) {
						cout << "相等" << endl;
					}
					else {
						cout << "不相等" << endl;
					}						
				}*/
				return area * computeKij(x, nx, mid, ny);
			}
		}
		else
		{
			if (fabs(mid(0) - x(0)) < yipusilon && fabs(mid(1) - x(1)) < yipusilon &&fabs(mid(2) - x(2)) < yipusilon) {
				Matrix3d zero;
				zero.setZero();
				return zero;
			}
			else {
				return area * computeHij(x, nx, mid, ny);
			}
		}
	}
	else {
		Vector3d p01 = (p0 + p1) / 2;
		Vector3d p02 = (p0 + p2) / 2;
		Vector3d p12 = (p1 + p2) / 2;		
		n--;
		Matrix3d res1 = digui(p0, p01, p02, n, area, flag, x, nx, ny);
		Matrix3d res2 = digui(p1, p01, p12, n, area, flag, x, nx, ny);
		Matrix3d res3 = digui(p2, p02, p12, n, area, flag, x, nx, ny);
		Matrix3d res4 = digui(p01, p12, p02, n, area, flag, x, nx, ny);
		return res1 + res2 + res3 + res4;
	}
	
}

