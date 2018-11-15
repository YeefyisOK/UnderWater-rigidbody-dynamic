#pragma once
class CMesh;
class Face;

class CObjReaderWirter
{
public:
	CObjReaderWirter();
	~CObjReaderWirter();

	static CMesh * Read(const std::string &path);
	static bool Write(CMesh *pMesh, const std::string &path);

private:
	static bool ReadFace(const std::string &strLine, std::list<Face> &fList, std::list<Face> &tList);
	static void ReadMtl(const std::string& mtlPath, std::string &texturePath);

	static bool ReadLine(std::ifstream &file, std::string &strLine, char comment = '#');
	static bool TrimLeft(std::string &str);
	static bool TrimRight(std::string &str);
	static bool Trim(std::string &str);

	static bool SplitTerms(const std::string &strLine,
		std::vector< std::string > &terms,
		char delimeter = ' ');
};

