#pragma once
#include <vector>
#include <glm/glm.hpp>

using namespace std;

class Loader
{
public:
	Loader();
	~Loader();

	vector< unsigned int > vertexIndices, uvIndices, normalIndices;
	vector< glm::vec3 > tmp_vertices;
	vector< glm::vec2 > tmp_uvs;
	vector< glm::vec3 > tmp_normals;

	bool LoadOBJ(const char* paths, vector<glm::vec3> & out_Vert, vector<glm::vec2> & out_UV, vector<glm::vec3> & out_Normals);

};

Loader::Loader()
{
}

Loader::~Loader()
{
}
