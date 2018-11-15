#include "stdafx.h"
#include "FileManager.h"


#include "ObjReaderWirter.h"
#include "Mesh.h"


#include "Model.h"

CFileManager::CFileManager()
{
}


CFileManager::~CFileManager()
{
}

std::string CFileManager::GetFileExt(const std::string &path)
{
	//extension
	std::string ext;
	auto pos = path.rfind(".");
	if (pos != std::string::npos)
	{
		ext = path.substr(pos + 1);
		std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
	}

	return ext;
}

bool CFileManager::Write(const std::string &path)
{
	//extension
	std::string ext = GetFileExt(path);

	CModel &model = CModel::GetInstance();
	if (ext == "obj")
	{
		CMesh *pMesh = model.GetMesh();
		if (pMesh == nullptr)
			return false;

		return CObjReaderWirter::Write(pMesh, path);
	}

	return false;
}

CEntity *CFileManager::Read(const std::string &path)
{
	//extension
	std::string ext = GetFileExt(path);

	if (ext == "obj")
	{
		CMesh *pMesh = CObjReaderWirter::Read(path);
		return pMesh;
	}
	

	return nullptr;
}
