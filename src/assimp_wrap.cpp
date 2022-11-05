// 3D World - AssImp Reader Wrapper
// by Frank Gennari
// 10/23/2014
// Reference: https://github.com/assimp/assimp

#include "3DWorld.h"
#include "model3d.h"

#define ENABLE_ASSIMP

#ifdef ENABLE_ASSIMP

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

vector3d  aiVector3D_to_vector3d(aiVector3D const &v) {return vector3d (v.x, v.y, v.z);}
colorRGBA aiColor4D_to_colorRGBA(aiColor4D  const &c) {return colorRGBA(c.r, c.g, c.b, c.a);}

xform_matrix aiMatrix4x4_to_xform_matrix(aiMatrix4x4 const &m) {
	xform_matrix M;
	M[0][0] = m.a1; M[0][1] = m.b1; M[0][2] = m.c1; M[0][3] = m.d1;
	M[1][0] = m.a2; M[1][1] = m.b2; M[1][2] = m.c2; M[1][3] = m.d2;
	M[2][0] = m.a3; M[2][1] = m.b3; M[2][2] = m.c3; M[2][3] = m.d3;
	M[3][0] = m.a4; M[3][1] = m.b4; M[3][2] = m.c4; M[3][3] = m.d4;
	return M;
}


// For reference, see: https://learnopengl.com/Model-Loading/Model
class file_reader_assimp {
	// input/output variables
	model3d &model;
	geom_xform_t cur_xf;
	string model_dir;
	bool load_animations=0;
	// internal loader state
	bool had_vertex_error=0;
	int space_count=0;
	map<string, unsigned> bone_name_to_index_map;

	int load_texture(aiMaterial const* const mat, aiTextureType const type, bool is_normal_map=0) {
		unsigned const count(mat->GetTextureCount(type));
		if (count == 0) return -1; // no texture
		// load only the first texture, as that's all we support
		aiString fn; // absolute path, not relative to the model file
		if (mat->GetTexture(type, 0, &fn) != AI_SUCCESS) return -1;
		bool const invert_y = 1;
		// is_alpha_mask=0, verbose=0, invert_alpha=0, wrap=1, mirror=0, force_grayscale=0
		return model.tmgr.create_texture((model_dir + fn.C_Str()), 0, 0, 0, 1, 0, 0, is_normal_map, invert_y);
	}

	void print_space() {
		for (int i = 0 ; i < space_count ; i++) {printf(" ");}
	}
	void print_assimp_matrix(aiMatrix4x4 const &m) {
		//xform_matrix const M(aiMatrix4x4_to_xform_matrix(m));
		//M.print();
		print_space(); printf("%f %f %f %f\n", m.a1, m.a2, m.a3, m.a4);
		print_space(); printf("%f %f %f %f\n", m.b1, m.b2, m.b3, m.b4);
		print_space(); printf("%f %f %f %f\n", m.c1, m.c2, m.c3, m.c4);
		print_space(); printf("%f %f %f %f\n", m.d1, m.d2, m.d3, m.d4);
	}

	unsigned get_bone_id(const aiBone* bone) {
		string const bone_name(bone->mName.C_Str());
		auto it(bone_name_to_index_map.find(bone_name));
		if (it != bone_name_to_index_map.end()) {return it->second;}
		unsigned const bone_id(bone_name_to_index_map.size()); // allocate an index for a new bone
		bone_name_to_index_map[bone_name] = bone_id;
		return bone_id;
	}
	void parse_single_bone(int bone_index, const aiBone* pBone, mesh_bone_data_t &bone_data, unsigned first_vertex_offset) {
		printf("      Bone %d: '%s' num vertices affected by this bone: %d\n", bone_index, pBone->mName.C_Str(), pBone->mNumWeights);
		unsigned const bone_id(get_bone_id(pBone));
		//printf("bone id %d\n", bone_id);
		print_assimp_matrix(pBone->mOffsetMatrix);

		for (unsigned i = 0; i < pBone->mNumWeights; i++) {
			//if (i == 0) {printf("\n");}
			const aiVertexWeight& vw = pBone->mWeights[i];
			//printf("       %d: vertex id %d weight %.2f\n", i, vw.mVertexId, vw.mWeight);
			unsigned const vertex_id(first_vertex_offset + vw.mVertexId);
			//printf("Vertex id %d ", vertex_id);
			assert(vertex_id < bone_data.vertex_to_bones.size());
			bone_data.vertex_to_bones[vertex_id].add(bone_id, vw.mWeight, had_vertex_error);
		}
		//printf("\n");
	}
	void parse_mesh_bones(const aiMesh* mesh, mesh_bone_data_t &bone_data, unsigned first_vertex_offset) {
		for (unsigned int i = 0; i < mesh->mNumBones; i++) {parse_single_bone(i, mesh->mBones[i], bone_data, first_vertex_offset);}
	}
	void process_mesh(aiMesh *mesh, const aiScene *scene) {
		assert(mesh != nullptr);
		if (!(mesh->mPrimitiveTypes & aiPrimitiveType_TRIANGLE)) return; // not a triangle mesh - skip for now (can be removed using options)
		vector<vert_norm_tc> verts(mesh->mNumVertices);
		vector<unsigned> indices;
		indices.reserve(3*mesh->mNumFaces);
		cube_t mesh_bcube;

		for (unsigned i = 0; i < mesh->mNumVertices; i++) { // process vertices
			vert_norm_tc &v(verts[i]);
			assert(mesh->mVertices != nullptr); // vertices are required
			assert(mesh->mNormals  != nullptr); // we specified normal creation, so these shouldbe non-null
			v.v = aiVector3D_to_vector3d(mesh->mVertices[i]); // position
			v.n = aiVector3D_to_vector3d(mesh->mNormals [i]); // normals
			cur_xf.xform_pos   (v.v);
			cur_xf.xform_pos_rm(v.n);

			if (mesh->mTextureCoords != nullptr && mesh->mTextureCoords[0] != nullptr) { // TCs are optional and default to (0,0); we only use the first of 8
				v.t[0] = mesh->mTextureCoords[0][i].x; 
				v.t[1] = mesh->mTextureCoords[0][i].y;
			}
			if (i == 0) {mesh_bcube.set_from_point(v.v);} else {mesh_bcube.union_with_pt(v.v);}
		} // for i
		assert(mesh->mFaces != nullptr);
		assert(mesh->mNumFaces > 0); // if there were verts, there must be faces

		for (unsigned i = 0; i < mesh->mNumFaces; i++) { // process faces/indices
			aiFace const& face(mesh->mFaces[i]);
			assert(face.mNumIndices == 3); // must be triangles
			for (unsigned j = 0; j < face.mNumIndices; j++) {indices.push_back(face.mIndices[j]);}
		}
		if (!mesh_bcube.is_all_zeros()) {model.union_bcube_with(mesh_bcube);}
		//if (mesh->mMaterialIndex >= 0) {} // according to the tutorial, this check should be done; but mMaterialIndex is unsigned, so it can't fail?
		material_t &mat(model.get_material(mesh->mMaterialIndex, 1)); // alloc_if_needed=1
		bool const is_new_mat(mat.empty());
		unsigned const first_vertex_offset(mat.add_triangles(verts, indices, 1)); // add_new_block=1; should return 0
		//cout << TXT(mesh->mName.C_Str()) << TXT(mesh->mNumVertices) << TXT(mesh->mNumFaces) << TXT(mesh->mNumBones) << endl;
		
		if (load_animations && mesh->HasBones()) { // handle bones
			mesh_bone_data_t &bone_data(mat.get_bone_data_for_last_added_tri_mesh());
			bone_data.vertex_to_bones.resize(first_vertex_offset + mesh->mNumVertices);
			parse_mesh_bones(mesh, bone_data, first_vertex_offset);
		}
		if (is_new_mat) { // process material if this is the first mesh using it
			assert(scene->mMaterials != nullptr);
			aiMaterial const* const material(scene->mMaterials[mesh->mMaterialIndex]);
			assert(material != nullptr);
			// setup and load textures
			mat.a_tid    = load_texture(material, aiTextureType_AMBIENT);
			mat.d_tid    = load_texture(material, aiTextureType_DIFFUSE);
			mat.s_tid    = load_texture(material, aiTextureType_SPECULAR);
			mat.bump_tid = load_texture(material, aiTextureType_NORMALS, 1); // is_normal_map=1; or aiTextureType_HEIGHT?
			//mat.refl_tid = load_texture(material, aiTextureType_REFLECTION); // unused
			// setup colors
			aiColor4D color;
			if (aiGetMaterialColor(material, AI_MATKEY_COLOR_AMBIENT,  &color) == AI_SUCCESS) {mat.ka = aiColor4D_to_colorRGBA(color);}
			if (aiGetMaterialColor(material, AI_MATKEY_COLOR_DIFFUSE,  &color) == AI_SUCCESS) {mat.kd = aiColor4D_to_colorRGBA(color);}
			if (aiGetMaterialColor(material, AI_MATKEY_COLOR_SPECULAR, &color) == AI_SUCCESS) {mat.ks = aiColor4D_to_colorRGBA(color);}
			if (aiGetMaterialColor(material, AI_MATKEY_COLOR_EMISSIVE, &color) == AI_SUCCESS) {mat.ke = aiColor4D_to_colorRGBA(color);}
			unsigned max1(1), max2(1), max3(1), max4(1);
			float shininess(0.0), strength(0.0), alpha(1.0);
			
			if (aiGetMaterialFloatArray(material, AI_MATKEY_SHININESS,          &shininess, &max1) == AI_SUCCESS &&
				aiGetMaterialFloatArray(material, AI_MATKEY_SHININESS_STRENGTH, &strength,  &max2) == AI_SUCCESS)
			{
				mat.ns = shininess * strength;
			}
			// check for dissolve, but skip if it's 0; might also want to look at AI_MATKEY_COLOR_TRANSPARENT
			if (aiGetMaterialFloatArray(material, AI_MATKEY_OPACITY, &alpha, &max3) == AI_SUCCESS && alpha > 0.0) {mat.alpha = alpha;}
			aiGetMaterialFloatArray(material, AI_MATKEY_TRANSMISSION_FACTOR, &mat.tr, &max4);
			//if (aiGetMaterialIntegerArray(mtl, AI_MATKEY_ENABLE_WIREFRAME, &wireframe, &max) == AI_SUCCESS) {}
			//if (aiGetMaterialIntegerArray(mtl, AI_MATKEY_TWOSIDED,         &two_sided, &max) == AI_SUCCESS) {}
			// illum? tf?
		}
	}  
	void process_node_recur(aiNode *node, const aiScene *scene) {
		assert(node != nullptr);
		// process all the node's meshes (if any), in tree order rather than simply iterating over mMeshes
		for (unsigned i = 0; i < node->mNumMeshes; i++) {process_mesh(scene->mMeshes[node->mMeshes[i]], scene);}
		// then do the same for each of its children
		for (unsigned i = 0; i < node->mNumChildren; i++) {process_node_recur(node->mChildren[i], scene);}
	} 
public:
	file_reader_assimp(model3d &model_, bool load_animations_=0) : model(model_), load_animations(load_animations_) {}

	bool read(string const &fn, geom_xform_t const &xf, bool recalc_normals, bool verbose) {
		cur_xf = xf;
		Assimp::Importer importer;
		// aiProcess_OptimizeMeshes
		// aiProcess_ValidateDataStructure - for debugging
		// aiProcess_ImproveCacheLocality - optional, but already supported by the model3d class
		// aiProcess_FindDegenerates, aiProcess_FindInvalidData - optional
		unsigned flags(aiProcess_Triangulate | aiProcess_SortByPType | aiProcess_FlipUVs | aiProcess_JoinIdenticalVertices |
			           aiProcess_FixInfacingNormals | aiProcess_GenUVCoords | aiProcess_OptimizeMeshes);
		// Note: here we treat the recalc_normals flag as using smooth normals; if the model already contains normals, they're always used
		flags |= (recalc_normals ? aiProcess_GenSmoothNormals : aiProcess_GenNormals);
		if (!load_animations) {flags |= aiProcess_PreTransformVertices | aiProcess_RemoveRedundantMaterials;}
		aiScene const* const scene(importer.ReadFile(fn, flags));
		
		if (scene == nullptr || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
			cerr << "AssImp Import Error: " << importer.GetErrorString() << endl;
			return 0;
		}
		model_dir = fn;
		while (!model_dir.empty() && model_dir.back() != '/' && model_dir.back() != '\\') {model_dir.pop_back();} // remove filename from end, but leave the slash
		process_node_recur(scene->mRootNode, scene);
		model.finalize(); // optimize vertices, remove excess capacity, compute bounding sphere, subdivide, compute LOD blocks
		model.load_all_used_tids();
		if (verbose) {cout << "bcube: " << model.get_bcube().str() << endl << "model stats: "; model.show_stats();}
		return 1;
	}
};

bool read_assimp_model(string const &filename, model3d &model, geom_xform_t const &xf, int recalc_normals, bool verbose) {
	//timer_t timer("Read AssImp Model");
	bool const load_animations = 1; // Note: incomplete
	file_reader_assimp reader(model, load_animations);
	return reader.read(filename, xf, recalc_normals, verbose);
}

#else // ENABLE_ASSIMP

bool read_assimp_model(string const &filename, model3d &model, geom_xform_t const &xf, int recalc_normals, bool verbose) {
	cerr << "Error: Assimp model import has not been enabled at compile time" << endl;
	return 0;
}

#endif
