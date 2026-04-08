#include "bromesh/io/obj.h"

#include <cstdio>
#include <cstring>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace bromesh {

MeshData loadOBJ(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "r");
    if (!f) return {};

    std::vector<float> tempPos;     // x,y,z triples
    std::vector<float> tempNorm;    // x,y,z triples
    std::vector<float> tempUV;      // u,v pairs

    MeshData mesh;

    // Map from (posIdx, uvIdx, normIdx) to output vertex index.
    // Use -1 for "not present".
    std::map<std::tuple<int,int,int>, uint32_t> vertexMap;

    char line[1024];
    while (std::fgets(line, sizeof(line), f)) {
        // Skip comments and empty lines
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;

        if (line[0] == 'v' && line[1] == ' ') {
            float x, y, z;
            if (std::sscanf(line + 2, "%f %f %f", &x, &y, &z) == 3) {
                tempPos.push_back(x);
                tempPos.push_back(y);
                tempPos.push_back(z);
            }
        } else if (line[0] == 'v' && line[1] == 'n' && line[2] == ' ') {
            float x, y, z;
            if (std::sscanf(line + 3, "%f %f %f", &x, &y, &z) == 3) {
                tempNorm.push_back(x);
                tempNorm.push_back(y);
                tempNorm.push_back(z);
            }
        } else if (line[0] == 'v' && line[1] == 't' && line[2] == ' ') {
            float u, v;
            if (std::sscanf(line + 3, "%f %f", &u, &v) == 2) {
                tempUV.push_back(u);
                tempUV.push_back(v);
            }
        } else if (line[0] == 'f' && line[1] == ' ') {
            // Parse face vertices
            struct FaceVert {
                int pos = -1, uv = -1, norm = -1;
            };
            std::vector<FaceVert> faceVerts;

            std::istringstream iss(line + 2);
            std::string token;
            while (iss >> token) {
                FaceVert fv;
                // Possible formats: v, v/vt, v/vt/vn, v//vn
                int vals[3] = {0, 0, 0};
                bool hasDoubleSlash = token.find("//") != std::string::npos;

                if (hasDoubleSlash) {
                    // v//vn
                    if (std::sscanf(token.c_str(), "%d//%d", &vals[0], &vals[2]) >= 2) {
                        fv.pos = vals[0] - 1;
                        fv.norm = vals[2] - 1;
                    }
                } else {
                    // Count slashes
                    int slashes = 0;
                    for (char c : token) if (c == '/') slashes++;

                    if (slashes == 0) {
                        // v
                        if (std::sscanf(token.c_str(), "%d", &vals[0]) == 1) {
                            fv.pos = vals[0] - 1;
                        }
                    } else if (slashes == 1) {
                        // v/vt
                        if (std::sscanf(token.c_str(), "%d/%d", &vals[0], &vals[1]) == 2) {
                            fv.pos = vals[0] - 1;
                            fv.uv = vals[1] - 1;
                        }
                    } else {
                        // v/vt/vn
                        if (std::sscanf(token.c_str(), "%d/%d/%d", &vals[0], &vals[1], &vals[2]) == 3) {
                            fv.pos = vals[0] - 1;
                            fv.uv = vals[1] - 1;
                            fv.norm = vals[2] - 1;
                        }
                    }
                }

                if (fv.pos >= 0) {
                    faceVerts.push_back(fv);
                }
            }

            // Resolve face vertices to output indices and triangulate as fan
            auto resolveVertex = [&](const FaceVert& fv) -> uint32_t {
                auto key = std::make_tuple(fv.pos, fv.uv, fv.norm);
                auto it = vertexMap.find(key);
                if (it != vertexMap.end()) return it->second;

                uint32_t idx = static_cast<uint32_t>(mesh.positions.size() / 3);
                vertexMap[key] = idx;

                // Position
                if (fv.pos >= 0 && (fv.pos * 3 + 2) < (int)tempPos.size()) {
                    mesh.positions.push_back(tempPos[fv.pos * 3 + 0]);
                    mesh.positions.push_back(tempPos[fv.pos * 3 + 1]);
                    mesh.positions.push_back(tempPos[fv.pos * 3 + 2]);
                } else {
                    mesh.positions.push_back(0); mesh.positions.push_back(0); mesh.positions.push_back(0);
                }

                // Normal
                if (fv.norm >= 0 && (fv.norm * 3 + 2) < (int)tempNorm.size()) {
                    mesh.normals.push_back(tempNorm[fv.norm * 3 + 0]);
                    mesh.normals.push_back(tempNorm[fv.norm * 3 + 1]);
                    mesh.normals.push_back(tempNorm[fv.norm * 3 + 2]);
                }

                // UV
                if (fv.uv >= 0 && (fv.uv * 2 + 1) < (int)tempUV.size()) {
                    mesh.uvs.push_back(tempUV[fv.uv * 2 + 0]);
                    mesh.uvs.push_back(tempUV[fv.uv * 2 + 1]);
                }

                return idx;
            };

            // Fan triangulation: 0-1-2, 0-2-3, 0-3-4, ...
            if (faceVerts.size() >= 3) {
                uint32_t v0 = resolveVertex(faceVerts[0]);
                for (size_t i = 1; i + 1 < faceVerts.size(); ++i) {
                    uint32_t v1 = resolveVertex(faceVerts[i]);
                    uint32_t v2 = resolveVertex(faceVerts[i + 1]);
                    mesh.indices.push_back(v0);
                    mesh.indices.push_back(v1);
                    mesh.indices.push_back(v2);
                }
            }
        }
    }

    std::fclose(f);
    return mesh;
}

bool saveOBJ(const MeshData& mesh, const std::string& path) {
    if (mesh.empty()) return false;

    FILE* f = std::fopen(path.c_str(), "w");
    if (!f) return false;

    size_t vc = mesh.vertexCount();
    bool hasN = mesh.hasNormals();
    bool hasUV = mesh.hasUVs();

    // Write positions
    for (size_t i = 0; i < vc; ++i) {
        std::fprintf(f, "v %f %f %f\n",
            mesh.positions[i * 3 + 0],
            mesh.positions[i * 3 + 1],
            mesh.positions[i * 3 + 2]);
    }

    // Write normals
    if (hasN) {
        for (size_t i = 0; i < vc; ++i) {
            std::fprintf(f, "vn %f %f %f\n",
                mesh.normals[i * 3 + 0],
                mesh.normals[i * 3 + 1],
                mesh.normals[i * 3 + 2]);
        }
    }

    // Write UVs
    if (hasUV) {
        for (size_t i = 0; i < vc; ++i) {
            std::fprintf(f, "vt %f %f\n",
                mesh.uvs[i * 2 + 0],
                mesh.uvs[i * 2 + 1]);
        }
    }

    // Write faces (1-based indices)
    size_t tc = mesh.triangleCount();
    for (size_t i = 0; i < tc; ++i) {
        uint32_t a = mesh.indices[i * 3 + 0] + 1;
        uint32_t b = mesh.indices[i * 3 + 1] + 1;
        uint32_t c = mesh.indices[i * 3 + 2] + 1;

        if (hasN && hasUV) {
            std::fprintf(f, "f %u/%u/%u %u/%u/%u %u/%u/%u\n",
                a, a, a, b, b, b, c, c, c);
        } else if (hasN) {
            std::fprintf(f, "f %u//%u %u//%u %u//%u\n",
                a, a, b, b, c, c);
        } else if (hasUV) {
            std::fprintf(f, "f %u/%u %u/%u %u/%u\n",
                a, a, b, b, c, c);
        } else {
            std::fprintf(f, "f %u %u %u\n", a, b, c);
        }
    }

    std::fclose(f);
    return true;
}

} // namespace bromesh
