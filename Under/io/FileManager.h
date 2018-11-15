#pragma once
class CEntity;

class CFileManager
{
public:
	CFileManager();
	~CFileManager();

	static bool Write(const std::string &path);

	static CEntity *Read(const std::string &path);

private:
	static std::string GetFileExt(const std::string &path);
};

