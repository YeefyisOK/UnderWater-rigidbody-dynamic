#include <Eigen/Dense>
#include"Point3D.h"
#include"KirchhoffF.h"
using namespace Eigen;
using namespace std;
void Subexpressions(double &w0, double &w1, double &w2, double &f1, double &f2, double &f3, double &g0, double &g1, double &g2) {
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
void Compute(MatrixXd face[], int numFaces, int index[], double mass) {
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
		Subexpressions(x0,x1,x2,f1x,f2x,f3x,g0x,g1x,g2x); 
		Subexpressions(y0,y1,y2,f1y,f2y,f3y,g0y,g1y,g2y);
		Subexpressions(z0,z1,z2,f1z,f2z,f3z,g0z,g1z,g2z);
	// update integrals 
		intg[0] += d0*f1x; intg[1] += d0*f2x; intg[2] += d1*f2y; intg[3] += d2*f2z;
		intg[4] += d0*f3x; intg[5] += d1*f3y; intg[6] += d2*f3z; 
		intg[7] += d0*(y0*g0x+y1*g1x+y2*g2x); intg[8] += d1*(z0*g0y+z1*g1y+z2*g2y); intg[9] += d2*(x0*g0z+x1*g1z+x2*g2z);
	} 
	for (int i = 0; i < 10; i++) 
		intg[i] *= mult[i];
	mass = intg[0];
	// center of mass
	cm.setData(intg[1]/mass ,intg[2]/mass ,intg[3]/mass);
	// inertia tensor relative to center of mass 
	inertia(0,0) = intg[5] + intg[6] - mass*(cm.getm_data(1)*cm.getm_data(1) + cm.getm_data(2)*cm.getm_data(2));
	inertia(1,1) = intg[4] + intg[6] - mass*(cm.getm_data(2)*cm.getm_data(2) + cm.getm_data(0)*cm.getm_data(0));
	inertia(2,2) = intg[4] + intg[5] - mass*(cm.getm_data(0)*cm.getm_data(0) + cm.getm_data(1)*cm.getm_data(1));
	inertia(0,1) = inertia(1,0) = -(intg[7] - mass*cm.getm_data(0)*cm.getm_data(1));
	inertia(1,2) = inertia(2,1) = -(intg[8] - mass*cm.getm_data(1)*cm.getm_data(2));
	inertia(2,0) = inertia(0,2) -(intg[9] - mass*cm.getm_data(2)*cm.getm_data(0));
	}
