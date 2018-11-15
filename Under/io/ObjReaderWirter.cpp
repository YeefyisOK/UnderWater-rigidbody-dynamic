#include "stdafx.h"
#include "ObjReaderWirter.h"

#include "Mesh.h"

class Vertex
{
public:
	Vertex()
	{
		m_coord[0] = m_coord[1] = m_coord[2] = 0.0;
	}
	Vertex(double coord[3])
	{
		m_coord[0] = coord[0];
		m_coord[1] = coord[1];
		m_coord[2] = coord[2];
	}

	double m_coord[3];
};

class TextCoord
{
public:
	TextCoord()
	{
		m_texCoord[0] = m_texCoord[1] = 0;
	}
	TextCoord(double coord[2])
	{
		m_texCoord[0] = coord[0];
		m_texCoord[1] = coord[1];
	}

	double m_texCoord[2];
};

class Face
{
public:
	Face(int vid[3])
	{
		m_vid[0] = vid[0];
		m_vid[1] = vid[1];
		m_vid[2] = vid[2];
	}

	int m_vid[3];
};

CObjReaderWirter::CObjReaderWirter()
{
}


CObjReaderWirter::~CObjReaderWirter()
{
}

CMesh *CObjReaderWirter::Read(const std::string &path)
{
	std::string strLine;

	std::list<Vertex> vList;
	std::list<Face> fList;
	std::list<Face> tList;   //texture index list
	std::vector<TextCoord> vecTextCoord;

	DWORD start = GetTickCount();

	Vertex vertex;
	TextCoord textCoord;
	std::string texturePath;

	std::ifstream file(path);
	if (!file.is_open())
	{
		return false;
	}

	std::string rootPath = path;
	auto nPos = path.rfind("\\");
	if (nPos == std::string::npos)
	{
		nPos = path.rfind("/");
	}
	if (nPos != std::string::npos)
	{
		rootPath = path.substr(0, nPos + 1);
	}
	while (ReadLine(file, strLine))
	{
		if (strLine.substr(0, 6) == "mtllib")
		{
			ReadMtl(rootPath + strLine.substr(7), texturePath);
			continue;
		}
		switch (strLine[0])
		{
		case 'v':
			switch (strLine[1])
			{
			case 'n': //normal
				break;
			case 't':  //texture
				sscanf_s(strLine.c_str(), "%*s %lf %lf", &textCoord.m_texCoord[0],
					&textCoord.m_texCoord[1]);
				vecTextCoord.push_back(textCoord);
				break;
			case ' ': //vertex
			case '\t':
				sscanf_s(strLine.c_str(), "%*s %lf %lf %lf", &vertex.m_coord[0],
					&vertex.m_coord[1], &vertex.m_coord[2]);
				vList.push_back(vertex);
				break;
			default:
				std::cerr << "Invalid parameter" << std::endl;
				break;
			}
			break;
		case 'f':
			if (!ReadFace(strLine, fList, tList))
				return NULL;
			break;
		default:
			break;
		}
	}

	DWORD dur = GetTickCount() - start;
	std::cout << dur << std::endl;

	if (vList.empty() || fList.empty())
		return NULL;


	//convert to mesh
	CMesh *pMesh = new CMesh();

	//vertex
	pMesh->SetNumOfVertices(vList.size());
	std::list<Vertex>::iterator iter = vList.begin();
	int i = 0;
	for (; iter != vList.end(); ++iter)
	{
		pMesh->SetVertexAt(i, iter->m_coord);
		++i;
	}

	//read the image
	cv::Mat image;
	if (!texturePath.empty())
	{
		image = cv::imread(rootPath + texturePath);
		cv::flip(image, image, 0);

	}
	pMesh->SetTexture(image);

	//set the texture coordinate

	//face
	pMesh->SetNumOfFaces(fList.size());
	std::list<Face>::iterator fIter = fList.begin();
	std::list<Face>::iterator tIter = tList.begin();

	bool bHasTexture = !image.empty() && (fList.size() == tList.size());

	i = 0;
	for (; fIter != fList.end(); ++fIter)
	{
		pMesh->SetFaceAt(i, fIter->m_vid);

		if (bHasTexture)
		{
			CVertex &v0 = pMesh->GetVertexAt(fIter->m_vid[0]);
			v0.SetTexCoord(vecTextCoord[tIter->m_vid[0]].m_texCoord);

			CVertex &v1 = pMesh->GetVertexAt(fIter->m_vid[1]);
			v1.SetTexCoord(vecTextCoord[tIter->m_vid[1]].m_texCoord);

			CVertex &v2 = pMesh->GetVertexAt(fIter->m_vid[2]);
			v2.SetTexCoord(vecTextCoord[tIter->m_vid[2]].m_texCoord);

			++tIter;
		}
		++i;
	}

	pMesh->CalculateEdgesFromFacesAndVertices();

	pMesh->ComputeFaceVertexNormal();
	return pMesh;
}

bool CObjReaderWirter::ReadFace(const std::string &strLine, std::list<Face> &fList, std::list<Face> &tList)
{
	int vid[3], tid[3];
	std::string::size_type start = 0, offset = 0, pos = 0;

	//ignore the keyword of 'f'
	start = strLine.find_first_of(' ') + 1;
	int i = 0;
	while (i < 3 && (offset = strLine.find_first_of(' ', start)) !=
		std::string::npos)
	{
		std::string str = strLine.substr(start, offset - start);
		start = offset + 1;

		if ((pos = str.find('/')) == std::string::npos)
		{
			vid[i] = atoi(str.c_str()) - 1;
		}
		else
		{
			std::string strV = str.substr(0, pos);
			vid[i] = atoi(strV.c_str()) - 1;

			//texture
			offset = pos + 1;
			pos = str.find_first_of('/', offset);
			if (pos > offset)
			{
				strV = str.substr(offset, pos - offset);
				//texture 
				tid[i] = atoi(strV.c_str()) - 1;
			}

		}
		++i;

	}

	if (i < 3)
	{
		std::string str = strLine.substr(start);
		if ((pos = str.find('/')) == std::string::npos)
		{
			vid[i] = atoi(str.c_str()) - 1;
		}
		else
		{
			std::string strV = str.substr(0, pos);
			vid[i] = atoi(strV.c_str()) - 1;

			//texture
			offset = pos + 1;
			pos = str.find_first_of('/', offset);
			if (pos > offset)
			{
				strV = str.substr(offset, pos - offset);
				//texture 
				tid[i] = atoi(strV.c_str()) - 1;
			}

		}
	}

	Face face(vid);
	fList.push_back(vid);

	Face textid(tid);
	tList.push_back(textid);

	return true;
}

void CObjReaderWirter::ReadMtl(const std::string& mtlPath, std::string &texturePath)
{
	//open the file
	std::ifstream file(mtlPath);
	if (!file.is_open())
	{
		std::cerr << "Failed to open material file!" << std::endl;
		return;
	}

	//find the key word "map_kd"
	std::string strLine;
	while (true)
	{
		do
		{
			if (!getline(file, strLine))   //end of file
				return;

			//trim left spaces
			Trim(strLine);
		} while (strLine.empty() || strLine[0] == '#');

		if (strLine.substr(0, 6) == "map_Kd")
		{
			
			texturePath = strLine.substr(7);
			return;
		}
	}

}


bool CObjReaderWirter::ReadLine(std::ifstream &file, std::string &strLine, char comment /*= '#'*/)
{
	do
	{
		if (!getline(file, strLine))   //end of file
			return false;

		//trim left spaces
		Trim(strLine);
	} while (strLine.empty() || strLine[0] == comment);

	return true;
}

inline bool CObjReaderWirter::TrimLeft(std::string &str)
{
	std::string::size_type pos = str.find_first_not_of(" \t");
	if (pos != std::string::npos)
		str.erase(0, pos);

	return true;
}

inline bool CObjReaderWirter::TrimRight(std::string &str)
{
	std::string::size_type pos = str.find_last_not_of(" \t");
	if (pos != std::string::npos)
		str.erase(pos + 1);
	else
		str.clear();

	return true;
}

inline bool CObjReaderWirter::Trim(std::string &str)
{
	std::string::size_type pos = str.find_last_not_of(" \t");
	if (pos != std::string::npos)
	{
		str.erase(pos + 1);
		pos = str.find_first_not_of(" \t");
		if (pos != std::string::npos)
			str.erase(0, pos);
	}
	else
		str.clear();

	return true;
}

bool CObjReaderWirter::SplitTerms(const std::string &strLine,
	std::vector< std::string > &terms,
	char delimeter /*= ' '*/)
{
	std::string::size_type start = 0, offset = 0;
	while ((offset = strLine.find_first_of(delimeter, start))
		!= std::string::npos)
	{
		terms.push_back(strLine.substr(start, offset - start));

		start = offset + 1;
	}

	terms.push_back(strLine.substr(start, strLine.length() - start));

	return true;
}

bool CObjReaderWirter::Write(CMesh *pMesh, const std::string &path)
{
	std::ofstream fs(path);
	if (!fs.is_open())
	{
		return false;
	}

	for (size_t i = 0; i < pMesh->GetNumOfVertices(); ++i)
	{
		CVertex &v = pMesh->GetVertexAt(i);
		fs << "v " << v[0] * 0.01 << " " << v[1] * 0.01 << " " << v[2] * 0.01 << std::endl;
	}

	for (size_t i = 0; i < pMesh->GetNumOfFaces(); ++i)
	{
		int *vId = pMesh->GetFaceAt(i).GetVerticeId();
		fs << "f " << vId[0] + 1 << " " << vId[1] + 1 << " " << vId[2] + 1 << std::endl;
	}

	fs.close();

	return true;
}
