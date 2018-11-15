#pragma once
#include "ArcBall.h"

class CPointCloud;
class CTMesh;
class CBBox;
class CSkeleton;
class CJoint;
class CMesh;
class CLineSeg;
class CObb;
class CBlendShape;

enum DisplayMode_e
{
	WireframeMode = 0,
	PointSetMode,
	SmoothMode
};

class COpenGLRender
{
public:
	COpenGLRender();
	~COpenGLRender();

	static COpenGLRender& GetInstance()
	{
		return m_render;
	}

	void RenderScene();

	void DrawModel();

	void DrawCoordSystem();

	void DrawCylinder(double r, double h, int nSlice);

	void DrawMesh(CMesh *pMesh);

	void SetWindowSize(int width, int height)
	{
		m_nWidth = width;
		m_nHeight = height;

		m_arcball.SetBounds(width, height);
	}

	CArcBall &GetArcball()
	{
		return m_arcball;
	}
	
	void Translate(const CVector3D &offset)
	{
		m_translate += offset;
	}

	void Scale(double scale)
	{
		m_dScale *= scale;
	}

	/////////

	void ZoomToFit();

	void SetDisplayMode(DisplayMode_e mode)
	{
		m_displayMode = mode;
	}

	DisplayMode_e GetDisplayMode() const
	{
		return m_displayMode;
	}
private:
	void CreateTexture(CMesh *pMesh);

private:
	CArcBall m_arcball;

	double m_dScale;
	CVector3D m_translate;

	int m_nWidth, m_nHeight;

	DisplayMode_e m_displayMode;

	static COpenGLRender m_render;
};

