#pragma once
#include <vector>
#include <glm/glm.hpp>

using namespace std;

bool LoadOBJ(const char* paths, vector<glm::vec3> & out_Vert, vector<glm::vec2> & out_UV, vector<glm::vec3> & out_Normals)
