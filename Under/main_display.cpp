#include<windows.h>
#include<iostream>
#include<sstream> 
#include<GL/glut.h>
#include<stdio.h>
#include"Kirchhoff.h"
#include "PIC.h"
#include"DynamicFormula.h"
#include <fstream>
#include<iostream>
using namespace std;
/*
定时器：glutTimerFunc(unsigned int millis, void (*func)(int value), int value);
重绘标志：glutPosRedisplay();
定时器第一个参数是每隔millis毫秒便调用func函数，并且创一个value参数进去
因为一个定时器只被调用一次，所以需要多次调用定时器
*/

//obj读取
string name = "H:\\MeshData\\cube.obj";
PIC m_pic;

void drawScene();
//窗口的大小
GLfloat windowWidth;
GLfloat windowHeight;
/*
//偏移量 旋转量
CVector3D temp_deltay(0,0,0);
CQuaternion q(0, 0, 0, 0);*/
Vector3d omega(0,1,0);
Vector3d velocity(0, -1, 0);
Matrix3d R = Matrix3d::Identity();

Vector3d y(0,0,0);
Vector3d ts(0,0,0);
Vector3d fs(0,-10,0);
MatrixXd K;
double delta_t=0.5;

DynamicFormula m_DF(omega,velocity,R,y,ts,fs,K,delta_t);

//double delta_t = 0.5;
void ReadPIC()
{
	ifstream ifs(name);//cube bunny Eight
	string s;
	Mian *f;
	Vertex *v;
	FaXiangLiang *vn;
	WenLi  *vt;
	while (getline(ifs, s))
	{
		if (s.length()<2)continue;
		if (s[0] == 'v'){
			if (s[1] == 't'){//vt 0.581151 0.979929 纹理
				istringstream in(s);
				vt = new WenLi();
				string head;
				in >> head >> vt->TU >> vt->TV;
				m_pic.VT.push_back(*vt);
			}
			else if (s[1] == 'n'){//vn 0.637005 -0.0421857 0.769705 法向量
				istringstream in(s);
				vn = new FaXiangLiang();
				string head;
				in >> head >> vn->NX >> vn->NY >> vn->NZ;
				m_pic.VN.push_back(*vn);
			}
			else{//v -53.0413 158.84 -135.806 点
				istringstream in(s);
				v = new Vertex();
				string head;
				in >> head >> v->X >> v->Y >> v->Z;
				m_pic.V.push_back(*v);
			}
		}
		else if (s[0] == 'f'){//f 2443//2656 2442//2656 2444//2656 面
			for (int k = s.size() - 1; k >= 0; k--){
				if (s[k] == '/')s[k] = ' ';
			}
			istringstream in(s);
			f = new Mian();
			string head;
			in >> head;
			int i = 0;
			while (i<3)
			{//
				if (m_pic.V.size() != 0)
				{
					in >> f->V[i];
					f->V[i] -= 1;
				}
				if (m_pic.VT.size() != 0)
				{
					in >> f->T[i];
					f->T[i] -= 1;
				}
				if (m_pic.VN.size() != 0)
				{
					in >> f->N[i];
					f->N[i] -= 1;
				}
				i++;
			}
			m_pic.F.push_back(*f);
		}
	}
}
void GLDraw()
{
	GLfloat ambient[] = { 0.3f, 0.3f, 0.3f, 1.0f };  // 环境强度
	GLfloat diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };  // 散射强度
	GLfloat specular[] = { 1.0f, 1.0f, 1.0f, 1.0f }; // 镜面强度
													 // 平行光源, GL_POSITION属性的最后一个参数为0
	GLfloat direction[] = { -3.0f, -3.4f, -8.8f, 0.0f };
	glLightfv(GL_LIGHT0, GL_POSITION, direction);
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular);

	glShadeModel(GL_SMOOTH);//中多边形内部各点的光照采用双线性插值的方法得到
	//设置多边形材质
	GLfloat mat_ambient[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat mat_diffuse[] = { 0.1, 0.5, 0.8, 1.0 };
	GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat no_mat[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat shininess = 50.0f;

	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_EMISSION, no_mat);
	glMaterialf(GL_FRONT, GL_SHININESS, shininess);

	//使能GL_COLOR_MATERIAL，而后利用glColorMaterial()函数指明glColor*()函数将影响到的材质属性
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_DIFFUSE);
	
	for (int i = 0; i<m_pic.F.size(); i++)
	{
		glBegin(GL_TRIANGLE_FAN);                            // 绘制三角形GL_TRIANGLES;GL_LINE_LOOP;GL_LINES;GL_POINTS
		if (m_pic.VT.size() != 0)glTexCoord2f(m_pic.VT[m_pic.F[i].T[0]].TU, m_pic.VT[m_pic.F[i].T[0]].TV);  //纹理    
		if (m_pic.VN.size() != 0)glNormal3f(m_pic.VN[m_pic.F[i].N[0]].NX, m_pic.VN[m_pic.F[i].N[0]].NY, m_pic.VN[m_pic.F[i].N[0]].NZ);//法向量
		glVertex3f(m_pic.V[m_pic.F[i].V[0]].X , m_pic.V[m_pic.F[i].V[0]].Y , m_pic.V[m_pic.F[i].V[0]].Z );        // 上顶点

		if (m_pic.VT.size() != 0)glTexCoord2f(m_pic.VT[m_pic.F[i].T[1]].TU, m_pic.VT[m_pic.F[i].T[1]].TV);  //纹理
		if (m_pic.VN.size() != 0)glNormal3f(m_pic.VN[m_pic.F[i].N[1]].NX, m_pic.VN[m_pic.F[i].N[1]].NY, m_pic.VN[m_pic.F[i].N[1]].NZ);//法向量
		glVertex3f(m_pic.V[m_pic.F[i].V[1]].X , m_pic.V[m_pic.F[i].V[1]].Y , m_pic.V[m_pic.F[i].V[1]].Z );        // 左下

		if (m_pic.VT.size() != 0)glTexCoord2f(m_pic.VT[m_pic.F[i].T[2]].TU, m_pic.VT[m_pic.F[i].T[2]].TV);  //纹理
		if (m_pic.VN.size() != 0)glNormal3f(m_pic.VN[m_pic.F[i].N[2]].NX, m_pic.VN[m_pic.F[i].N[2]].NY, m_pic.VN[m_pic.F[i].N[2]].NZ);//法向量
		glVertex3f(m_pic.V[m_pic.F[i].V[2]].X , m_pic.V[m_pic.F[i].V[2]].Y , m_pic.V[m_pic.F[i].V[2]].Z );        // 右下
		glEnd();

		glColor3f(0.0, 1.0, 0.0);     
		// 三角形绘制结束    
		/*if(m_pic.VN.size()!=0){
		glBegin(GL_LINES);                            // 绘制三角形
		glVertex3f(m_pic.V[m_pic.F[i].V[0]].X/YU,m_pic.V[m_pic.F[i].V[0]].Y/YU, m_pic.V[m_pic.F[i].V[0]].Z/YU);        // 上顶点
		glVertex3f(m_pic.V[m_pic.F[i].V[0]].X/YU+m_pic.VN[m_pic.F[i].N[0]].NX
		,m_pic.V[m_pic.F[i].V[0]].Y/YU+m_pic.VN[m_pic.F[i].N[0]].NY
		, m_pic.V[m_pic.F[i].V[0]].Z/YU+m_pic.VN[m_pic.F[i].N[0]].NZ);                    // 左下
		glEnd();                                // 三角形绘制结束
		}*/
	}
}
void init() {
	ReadPIC();
	//基尔霍夫张量
	CKirchhoff m_K(m_pic);
	K = m_K.computeK();//初始时得到K矩阵
	m_DF.setK(K);
	//glClearColor(0.0, 0.0, 0.0, 1.0);           //设置背景颜色  
	glClearColor(0.75f, 0.75f, 0.75f, 0.0f);
	//深度测试的相关设置 
	glClearDepth(1.0);                    //设置深度缓存的初始值 
	glDepthFunc(GL_LEQUAL);           //深度测试的方法 
	glEnable(GL_DEPTH_TEST);          //启用深度测试
	//光照
//材质反光性设置
	GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };  //镜面反射参数
	GLfloat mat_shininess[] = { 50.0 };               //高光指数
	GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };
	GLfloat white_light[] = { 1.0, 1.0, 1.0, 1.0 };   //灯位置(1,1,1), 最后1-开关
	GLfloat Light_Model_Ambient[] = { 0.2, 0.2, 0.2, 1.0 }; //环境光参数

	glClearColor(0.0, 0.0, 0.0, 0.0);  //背景色
	glShadeModel(GL_SMOOTH);           //多变性填充模式

	//材质属性
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

	//灯光设置
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, white_light);   //散射光属性
	glLightfv(GL_LIGHT0, GL_SPECULAR, white_light);  //镜面反射光
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Light_Model_Ambient);  //环境光参数

	glEnable(GL_LIGHTING);   //开关:使用光
	glEnable(GL_LIGHT0);     //打开0#灯
	glEnable(GL_DEPTH_TEST); //打开深度测试

}
//窗口重绘时的响应函数 
void display() {
	//用前面设置的背景色和深度值清除颜色缓存和深度缓存 
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//初始化模型视图矩阵 
	drawScene();            //几何图元的绘制 
	glFlush();              //绘制结束，Flush 当前渲染流水线 
	glutSwapBuffers();       //交换前后缓存（只用于双缓冲的模式） 
}
void DrawCylinder(double r, double h, int nSlice)
{
	//glTranslated(0, 0, -h / 2);

#define  PI   3.1415926
	double delta = PI * 2 / nSlice;
	double angle = delta;

	glColor3f(1.0f, 0.0f, 0.0f);

	//top
	glBegin(GL_TRIANGLE_FAN);
	glNormal3f(0, 0, 1);
	glVertex3d(0, 0, h);
	glVertex3d(r, 0, h);
	for (int i = 1; i < nSlice; ++i)
	{
		glVertex3d(r * cos(angle), r * sin(angle), h);
		angle += delta;
	}
	glVertex3d(r, 0, h);
	glEnd();

	//bottom
	angle = delta;
	glColor3f(0.0f, 1.0f, 0.0f);

	glBegin(GL_TRIANGLE_FAN);
	glNormal3f(0, 0, -1);
	glVertex3d(0, 0, 0);
	glVertex3d(r, 0, 0);
	for (int i = 1; i < nSlice; ++i)
	{
		glVertex3d(r * cos(angle), r * sin(angle), 0);
		angle += delta;
	}
	glVertex3d(r, 0, 0);
	glEnd();

	//cylinder
	angle = delta;
	glColor3f(0.0f, 0.0f, 1.0f);

	glBegin(GL_QUAD_STRIP);
	glNormal3f(1, 0, 0);
	glVertex3d(r, 0, h);
	glVertex3d(r, 0, 0);
	for (int i = 1; i < nSlice; ++i)
	{
		double c = cos(angle);
		double s = sin(angle);
		glNormal3d(c, s, 0);

		double x = r * c;
		double y = r * s;
		glVertex3d(x, y, h);
		glVertex3d(x, y, 0);

		angle += delta;
	}
	glNormal3f(1, 0, 0);
	glVertex3d(r, 0, h);
	glVertex3d(r, 0, 0);
	glEnd();
}
void drawScene()           //绘制
{	
	//glRotated()
	//glTranslated(m_translate[0], m_translate[1], m_translate[2]);
	glTranslated(m_DF.temp_deltay(0), m_DF.temp_deltay(1), m_DF.temp_deltay(2));
	glRotated(m_DF.q.w(), m_DF.q.x(), m_DF.q.y(), m_DF.q.z());
	glColor3f(0.0, 1.0, 0.0);     //绿
	GLDraw();
	//DrawCylinder(1,2, 32);
}
//窗口大小发生变化时的响应函数 
void reshape(int width, int height) {
	glViewport(0, 0, width, height);      //设置视窗大小 
	//设置视景体大小 
	glMatrixMode(GL_PROJECTION);
	double ratio = (double)width / height;
	glLoadIdentity();
	//gluPerspective(60, ratio, 1, 1000);
	glOrtho(-5, 5, -5, 5, -10, 10);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(4, 0, -2, 0, 0, -2, 0, 0, 1);
}

void TimerFunction(int value)
{
	m_DF.nextTime();
	/*
	glutPostRedisplay 标记当前窗口需要重新绘制。通过glutMainLoop下一次循环时，
	窗口显示将被回调以重新显示窗口的正常面板。多次调用glutPostRedisplay，
	在下一个显示回调只产生单一的重新显示回调
	*/
	glutPostRedisplay(); //标志重新绘制
	glutTimerFunc(1000, TimerFunction, 1);
}

int main(int argc, char* argv[]) {
	glutInit(&argc, argv);    //GLUT 库的初始化 
	//显示模式初始化：颜色格式——GLUT_RGBA 
	//           单缓冲或双缓冲——GLUT_SINGLE 或者 GLUT_DOUBLE //          是否使用深度缓存——GLUT_DEPTH  
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowPosition(100, 100);          //窗口起始位置 
	glutInitWindowSize(500, 500);             //窗口大小 
	glutCreateWindow("UnderWaterRidgebody");       //创建窗口并指定窗口的名称 

	init();                                 //OpenGL 的初始化设置 

	//设置消息响应的回调函数 
	glutDisplayFunc(display);               //设置窗口重绘时的响应函数 
	glutReshapeFunc(reshape);              //设置窗口大小发生变化时的响应函数 

	//定时器  每500毫秒触发一次
	glutTimerFunc(500, TimerFunction, 1);

	glutMainLoop();                   //消息循环：获取消息，分发消息，响应消息 
	return 0;
} //
