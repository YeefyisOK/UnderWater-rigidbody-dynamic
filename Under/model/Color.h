#pragma once
class CColor
{
public:
	CColor(byte red = 255, byte green = 255, byte blue = 255);
	~CColor();

	void Set(byte red, byte green, byte blue)
	{
		m_red = red;
		m_green = green;
		m_blue = blue;
	}

	byte Red() const
	{
		return m_red;
	}

	byte Green() const
	{
		return m_red;
	}

	byte Blue() const
	{
		return m_red;
	}

	float * Data() const
	{
		static float color[4];
		color[0] = m_red / 255.0f;
		color[1] = m_green / 255.0f;
		color[2] = m_blue / 255.0f;
		color[3] = 1.0f;
		return color;
	}
private:
	byte m_red, m_green, m_blue;
};

