// 3D World - AssImp Reader Wrapper
// by Frank Gennari
// 10/23/2014
// Reference: https://github.com/assimp/assimp

#include "3DWorld.h"
#include "model3d.h"

//#include <assimp/cimport.h>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

class file_reader_assimp {
	model3d &model;

	void process_mesh(aiMesh *mesh, const aiScene *scene){
		for(unsigned int i = 0; i < mesh->mNumVertices; i++) {
			// process vertex positions, normals and texture coordinates
			// TODO
		}
		// process indices
		// TOO
		// process material
		if(mesh->mMaterialIndex >= 0) {
			// TODO
		}
	}  
	void process_node(aiNode *node, const aiScene *scene) {
		// process all the node's meshes (if any)
		for (unsigned int i = 0; i < node->mNumMeshes; i++) {process_mesh(scene->mMeshes[node->mMeshes[i]], scene);}
		// then do the same for each of its children
		for(unsigned int i = 0; i < node->mNumChildren; i++) {process_node(node->mChildren[i], scene);}
	} 
public:
	file_reader_assimp(model3d &model_) : model(model_) {}

	bool read(string const &fn, geom_xform_t const &xf, bool verbose) {
		Assimp::Importer importer;
		// aiProcess_OptimizeMeshes
		aiScene const* const scene(importer.ReadFile(fn, (aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_GenNormals)));
		
		if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
			cerr << "AssImp Import Error: " << importer.GetErrorString() << endl;
			return 0;
		}
		//directory = fn.substr(0, fn.find_last_of('/'));
		process_node(scene->mRootNode, scene);
		return 1;
	}
};

bool read_assimp_model(string const &filename, model3d &model, geom_xform_t const &xf, int use_vertex_normals, bool verbose) {
	timer_t timer("Read AssImp Model");
	file_reader_assimp reader(model);
	return reader.read(filename, xf, verbose);
}
