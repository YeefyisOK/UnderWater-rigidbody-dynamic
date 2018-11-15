#include "stdafx.h"
#include "OpenGLRender.h"
#include "Model.h"
#include "Mesh.h"


COpenGLRender COpenGLRender::m_render;

COpenGLRender::COpenGLRender()
{
	m_dScale = 1.0;
	m_displayMode = SmoothMode;
}


COpenGLRender::~COpenGLRender()
{
}

void COpenGLRender::RenderScene()
{

	//coordinate system
	glViewport(0, 0, 100, 100);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 10.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0, 0, 3, 0, 0, 0, 0, 1, 0);

	//rotate
	double *pData = m_arcball.GetRotationData();
	glMultMatrixd(pData);

	DrawCoordSystem();

	//cylinder
	glViewport(0, 0, m_nWidth, m_nHeight);
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	//scale and translation

	glScaled(m_dScale, m_dScale, m_dScale);
	glTranslated(m_translate[0], m_translate[1], m_translate[2]);

	//DrawCylinder(0.5, 1.0, 32);	
	DrawModel();

	glPushMatrix();
	glScaled(0.1 / m_dScale, 0.1 / m_dScale, 0.1 / m_dScale);
	glDisable(GL_DEPTH_TEST);
	DrawCoordSystem();
	glEnable(GL_DEPTH_TEST);
	glPopMatrix();
}

void COpenGLRender::ZoomToFit()
{
	const CBBox &bbox = CModel::GetInstance().GetBBox();
	if (bbox.SquareDiagonalLength() < 1.0e-16)
	{
		return;
	}
	m_translate = CPoint3D(0, 0, 0) - bbox.Center();
	m_dScale = 1.0 / bbox.DiagonalLength();
}

void COpenGLRender::DrawModel()
{
	CModel &model = CModel::GetInstance();

	for (int i = 0; i < model.GetNumOfEntities(); ++i)
	{
		CEntity *pEntity = model.GetEntityAt(i);

		if (pEntity->GetObjectType() == MESH_TYPE)
		{
			DrawMesh((CMesh *)pEntity);
		}
	}
}



void COpenGLRender::DrawMesh(CMesh *pMesh)
{
	if (!pMesh->IsShown())
	{
		return;
	}
	glEnable(GL_LIGHTING);	

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, pMesh->GetColor().Data());

	if (m_displayMode == SmoothMode)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
	else if (m_displayMode == WireframeMode)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}
	else
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
		glPointSize(2.0f);
	}
	

	if (pMesh->HasTexture())
	{
		//
		glEnable(GL_TEXTURE_2D);
		if (!pMesh->IsTextureIntited())
		{			
			CreateTexture(pMesh);			
		}

		glBindTexture(GL_TEXTURE_2D, pMesh->GetTextureId());		
		glBegin(GL_TRIANGLES);
		for (size_t i = 0; i < pMesh->GetNumOfFaces(); ++i)
		{
			int *vId = pMesh->GetFaceAt(i).GetVerticeId();

			CVertex &v0 = pMesh->GetVertexAt(vId[0]);
			glNormal3dv(v0.GetNormal().Data());
			glTexCoord2dv(v0.GetTexCoord());
			glVertex3dv(v0.Data());

			CVertex &v1 = pMesh->GetVertexAt(vId[1]);
			glNormal3dv(v1.GetNormal().Data());
			glTexCoord2dv(v1.GetTexCoord());
			glVertex3dv(v1.Data());

			CVertex &v2 = pMesh->GetVertexAt(vId[2]);
			glNormal3dv(v2.GetNormal().Data());
			glTexCoord2dv(v2.GetTexCoord());
			glVertex3dv(v2.Data());
		}
		glEnd();

		glDisable(GL_TEXTURE_2D);
	} 
	else
	{
		glBegin(GL_TRIANGLES);				
		for (size_t i = 0; i < pMesh->GetNumOfFaces(); ++i)
		{
			int *vId = pMesh->GetFaceAt(i).GetVerticeId();

			CVertex &v0 = pMesh->GetVertexAt(vId[0]);
			glNormal3dv(v0.GetNormal().Data());
			glVertex3dv(v0.Data());

			CVertex &v1 = pMesh->GetVertexAt(vId[1]);
			glNormal3dv(v1.GetNormal().Data());
			glVertex3dv(v1.Data());

			CVertex &v2 = pMesh->GetVertexAt(vId[2]);
			glNormal3dv(v2.GetNormal().Data());
			glVertex3dv(v2.Data());
		}
		glEnd();

		//edges
// 		glDisable(GL_LIGHTING);
// 		glLineWidth(1.0f);
// 		glColor3f(0.0f, 1.0f, 0.0f);
// 		glBegin(GL_LINES);
// 		for (size_t i = 0; i < pMesh->GetNumOfEdges(); ++i)
// 		{
// 			int *vId = pMesh->GetEdgeAt(i).GetVerticeId();
// 			glVertex3dv(pMesh->GetVertexAt(vId[0]).Data());
// 			glVertex3dv(pMesh->GetVertexAt(vId[1]).Data());
// 		}
// 		glEnd();			
	}
	
}

void COpenGLRender::CreateTexture(CMesh *pMesh)
{
	//
	unsigned int textureId = 0;
	glGenTextures(1, &textureId);

	glBindTexture(GL_TEXTURE_2D, textureId);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	cv::Mat &image = pMesh->GetTexture();
	int nChannels = image.channels();
	GLenum format = GL_BGR_EXT;
	if (nChannels == 4)
	{
		format = GL_BGRA_EXT;
	}
	gluBuild2DMipmaps(GL_TEXTURE_2D, nChannels , image.cols, image.rows, format, GL_UNSIGNED_BYTE, image.data);

	pMesh->SetTextureId(textureId);
}


void COpenGLRender::DrawCylinder(double r, double h, int nSlice)
{
//#define  PI   3.1415926
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

void COpenGLRender::DrawCoordSystem()
{
	glDisable(GL_LIGHTING);

	double len = 0.6, radius = 0.1, height = 0.4;

	//z axis
	glColor3f(0, 0, 1.0f);

	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 0, len);
	glEnd();

	GLUquadric *pQuad = gluNewQuadric();

	glPushMatrix();
	glTranslated(0, 0, len);
	gluCylinder(pQuad, radius, 0, height, 32, 32);
	glPopMatrix();

	//x axis
	glColor3f(1.0f, 0, 0.0f);

	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(len, 0, 0);
	glEnd();

	glPushMatrix();
	glTranslated(len, 0, 0);
	glRotated(90, 0, 1, 0);
	gluCylinder(pQuad, radius, 0, height, 32, 32);
	glPopMatrix();

	//y axis
	glColor3f(0.0f, 1.0f, 0.0f);

	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(0, len, 0);
	glEnd();

	glPushMatrix();
	glTranslated(0, len, 0);
	glRotated(-90, 1, 0, 0);
	gluCylinder(pQuad, radius, 0, height, 32, 32);
	glPopMatrix();

	gluDeleteQuadric(pQuad);

	glEnable(GL_LIGHTING);
}