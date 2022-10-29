// 3D World - AssImp Reader Wrapper
// by Frank Gennari
// 10/23/2014
// Reference: https://github.com/assimp/assimp

#include "3DWorld.h"
#include "model3d.h"

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

vector3d aiVector3D_to_vector3d(aiVector3D const &v) {return vector3d(v.x, v.y, v.z);}

// For reference, see: https://learnopengl.com/Model-Loading/Model
class file_reader_assimp {
	model3d &model;

	int load_texture(aiMaterial const* const mat, aiTextureType const type, bool is_normal_map=0) {
		unsigned const count(mat->GetTextureCount(type));
		if (count == 0) return -1; // no texture
		// load only the first texture, as that's all we support
		aiString fn; // TODO: is this absolute, or relative to the model file?
		mat->GetTexture(type, 0, &fn);
		// is_alpha_mask=0, verbose=1, invert_alpha=0, wrap=1, mirror=0, force_grayscale=0
		return model.tmgr.create_texture(fn.C_Str(), 0, 1, 0, 1, 0, 0, is_normal_map);
	}
	void process_mesh(aiMesh *mesh, const aiScene *scene) {
		vector<vert_norm_tc> verts(mesh->mNumVertices);
		vector<unsigned> indices;
		indices.reserve(3*mesh->mNumFaces);

		for (unsigned i = 0; i < mesh->mNumVertices; i++) { // process vertices
			vert_norm_tc &v(verts[i]);
			v.v = aiVector3D_to_vector3d(mesh->mVertices[i]); // position
			v.n = aiVector3D_to_vector3d(mesh->mNormals [i]); // normals

			if (mesh->mTextureCoords[0]) { // texture coordinates are optional and default to (0,0); we only use the first of 8
				v.t[0] = mesh->mTextureCoords[0][i].x; 
				v.t[1] = mesh->mTextureCoords[0][i].y;
			}
		} // for i
		for (unsigned i = 0; i < mesh->mNumFaces; i++) { // process faces/indices
			aiFace const& face(mesh->mFaces[i]);
			assert(face.mNumIndices == 3); // must be triangles
			for (unsigned j = 0; j < face.mNumIndices; j++) {indices.push_back(face.mIndices[j]);}
		}
		//if (mesh->mMaterialIndex >= 0) {} // according to the tutorial, this check should be done; but mMaterialIndex is unsigned, so it can't fail?
		material_t &mat(model.get_material(mesh->mMaterialIndex, 1)); // alloc_if_needed=1
		bool const is_new_mat(mat.empty());
		mat.add_triangles(verts, indices, 1); // add_new_block=1
		
		if (is_new_mat) { // process material if this is the first mesh using it
			aiMaterial const* const material(scene->mMaterials[mesh->mMaterialIndex]);
			mat.a_tid    = load_texture(material, aiTextureType_AMBIENT);
			mat.d_tid    = load_texture(material, aiTextureType_DIFFUSE);
			mat.s_tid    = load_texture(material, aiTextureType_SPECULAR);
			mat.bump_tid = load_texture(material, aiTextureType_NORMALS, 1); // is_normal_map=1; or aiTextureType_HEIGHT?
			//mat.refl_tid = load_texture(material, aiTextureType_REFLECTION); // unused
			// I guess the colors remain at the defaults?
		}
	}  
	void process_node_recur(aiNode *node, const aiScene *scene) {
		// process all the node's meshes (if any), in tree order rather than simply iterating over mMeshes
		for (unsigned i = 0; i < node->mNumMeshes; i++) {process_mesh(scene->mMeshes[node->mMeshes[i]], scene);}
		// then do the same for each of its children
		for (unsigned i = 0; i < node->mNumChildren; i++) {process_node_recur(node->mChildren[i], scene);}
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
		process_node_recur(scene->mRootNode, scene);
		return 1;
	}
};

bool read_assimp_model(string const &filename, model3d &model, geom_xform_t const &xf, int use_vertex_normals, bool verbose) {
	timer_t timer("Read AssImp Model");
	file_reader_assimp reader(model);
	return reader.read(filename, xf, verbose);
}
