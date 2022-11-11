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

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <fstream>

vector3d  aiVector3D_to_vector3d(aiVector3D const &v) {return vector3d (v.x, v.y, v.z);}
colorRGBA aiColor4D_to_colorRGBA(aiColor4D  const &c) {return colorRGBA(c.r, c.g, c.b, c.a);}
glm::vec3 aiVector3D_to_glm_vec3(aiVector3D const &v) {return glm::vec3(v.x, v.y, v.z);}
xform_matrix aiMatrix4x4_to_xform_matrix(aiMatrix4x4  const &m) {return xform_matrix(glm::transpose(glm::make_mat4(&m.a1)));}
glm::mat3    aiMatrix3x3_to_glm_mat3    (aiMatrix3x3  const &m) {return glm::transpose(glm::make_mat3(&m.a1));}
glm::quat    aiQuaternion_to_glm_quat   (aiQuaternion const &q) {return glm::quat(q.w, q.x, q.y, q.z);}
void print_assimp_matrix(aiMatrix4x4 const &m) {aiMatrix4x4_to_xform_matrix(m).print();}


// For reference, see: https://learnopengl.com/Model-Loading/Model
// Also: https://github.com/emeiri/ogldev
// Also: http://www.xphere.me/2019/05/bones-animation-with-openglassimpglm/

void model_anim_t::anim_data_t::init(unsigned np, unsigned nr, unsigned ns) {
	assert(pos.empty() && scale.empty() && rot.empty()); // can only call init() once
	assert(np > 0);
	assert(nr > 0);
	assert(ns > 0);
	pos.reserve(np);
	rot.reserve(nr);
	scale.reserve(ns);
}

unsigned model_anim_t::get_bone_id(string const &bone_name) {
	auto it(bone_name_to_index_map.find(bone_name));
	if (it != bone_name_to_index_map.end()) {return it->second;}
	unsigned const bone_id(bone_name_to_index_map.size()); // allocate an index for a new bone
	bone_name_to_index_map[bone_name] = bone_id;
	return bone_id;
}
bool model_anim_t::update_bone_transform(string const &node_name, xform_matrix const &global_transform) {
	auto it(bone_name_to_index_map.find(node_name));
	if (it == bone_name_to_index_map.end()) return 0; // not found
	unsigned const bone_index(it->second);
	assert(bone_index < bone_transforms.size() && bone_index < bone_offset_matrices.size());
	bone_transforms[bone_index] = global_inverse_transform * global_transform * bone_offset_matrices[bone_index];
	return 1;
}
vector3d model_anim_t::calc_interpolated_position(float anim_time, anim_data_t const &A) const {
	assert(!A.pos.empty());
	if (A.pos.size() == 1) {return A.pos[0].v;} // single value, no interpolation

	for (unsigned i = 0; i+1 < A.pos.size(); ++i) {
		anim_vec3_val_t const &cur(A.pos[i]), &next(A.pos[i+1]);
		if (anim_time >= next.time) continue; // not yet
		float const t((anim_time - cur.time) / (next.time - cur.time));
		assert(t >= 0.0f && t <= 1.0f);
		return cur.v + t*(next.v - cur.v);
	} // for i
	assert(0);
	return zero_vector; // never gets here
}
glm::quat model_anim_t::calc_interpolated_rotation(float anim_time, anim_data_t const &A) const {
	assert(!A.rot.empty());
	if (A.rot.size() == 1) {return A.rot[0].q;} // single value, no interpolation

	for (unsigned i = 0; i+1 < A.rot.size(); ++i) {
		anim_quat_val_t const &cur(A.rot[i]), &next(A.rot[i+1]);
		if (anim_time >= next.time) continue; // not yet
		float const t((anim_time - cur.time) / (next.time - cur.time));
		assert(t >= 0.0f && t <= 1.0f);
		//aiQuaternion::Interpolate(out, start_rot, end_rot, factor);
		//out = out.Normalize();
		//return glm::normalize(cur.q + t*(next.q - cur.q));
		return glm::normalize(glm::slerp(cur.q, next.q, t));
	} // for i
	assert(0);
	return glm::quat(); // never gets here
}
vector3d model_anim_t::calc_interpolated_scale(float anim_time, anim_data_t const &A) const {
	assert(!A.scale.empty());
	if (A.scale.size() == 1) {return A.scale[0].v;} // single value, no interpolation

	for (unsigned i = 0; i+1 < A.scale.size(); ++i) {
		anim_vec3_val_t const &cur(A.scale[i]), &next(A.scale[i+1]);
		if (anim_time >= next.time) continue; // not yet
		float const t((anim_time - cur.time) / (next.time - cur.time));
		assert(t >= 0.0f && t <= 1.0f);
		return cur.v + t*(next.v - cur.v);
	} // for i
	assert(0);
	return zero_vector; // never gets here
}
void model_anim_t::transform_node_hierarchy_recur(float anim_time, animation_t const &animation, unsigned node_ix, xform_matrix const &parent_transform) {
	assert(node_ix < anim_nodes.size());
	anim_node_t const &node(anim_nodes[node_ix]);
	xform_matrix node_transform(node.transform); // defaults to node transform; may be overwritten below
	auto it(animation.anim_data.find(node.name));

	if (it != animation.anim_data.end()) { // found
		anim_data_t const &A(it->second);
		glm::mat4 const translation(glm::translate(glm::mat4(1.0), vec3_from_vector3d(calc_interpolated_position(anim_time, A))));
		glm::mat4 const rotation(glm::toMat4(calc_interpolated_rotation(anim_time, A)));
		glm::mat4 const scaling(glm::scale(glm::mat4(1.0), vec3_from_vector3d(calc_interpolated_scale(anim_time, A))));
		node_transform = translation * rotation * scaling;
	}
	xform_matrix const global_transform(parent_transform * node_transform);
	update_bone_transform(node.name, global_transform);
	for (unsigned i : node.children) {transform_node_hierarchy_recur(anim_time, animation, i, global_transform);}
}
void model_anim_t::get_bone_transforms(unsigned anim_id, float cur_time) {
	assert(anim_id < animations.size());
	animation_t const &animation(animations[anim_id]);
	float const time_in_ticks(cur_time * animation.ticks_per_sec);
	float const anim_time(fmod(time_in_ticks, animation.duration));
	transform_node_hierarchy_recur(anim_time, animation, 0, root_transform); // root node is 0
}

class file_reader_assimp {
	model3d &model;
	geom_xform_t cur_xf;
	string model_dir;
	bool load_animations=0, had_vertex_error=0;

	int load_texture(aiMaterial const* const mat, aiTextureType const type, bool is_normal_map=0) {
		unsigned const count(mat->GetTextureCount(type));
		if (count == 0) return -1; // no texture
		// load only the first texture, as that's all we support
		aiString fn; // absolute path, not relative to the model file
		if (mat->GetTexture(type, 0, &fn) != AI_SUCCESS) return -1;
		// is_alpha_mask=0, verbose=0, invert_alpha=0, wrap=1, mirror=0, force_grayscale=0
		return model.tmgr.create_texture((model_dir + fn.C_Str()), 0, 0, 0, 1, 0, 0, is_normal_map);
	}

	aiNodeAnim const *find_node_anim(aiAnimation const *const pAnimation, string const &node_name) {
		for (unsigned i = 0; i < pAnimation->mNumChannels; i++) {
			aiNodeAnim const *const anim(pAnimation->mChannels[i]);
			if (string(anim->mNodeName.data) == node_name) return anim;
		}
		return NULL;
	}
	unsigned extract_animation_data_recur(aiScene const *const scene, aiNode const *const node, model_anim_t &model_anim) {
		string const node_name(node->mName.data);
		unsigned const node_ix(model_anim.anim_nodes.size());
		model_anim.anim_nodes.emplace_back(node_name, aiMatrix4x4_to_xform_matrix(node->mTransformation));

		for (unsigned a = 0; a < scene->mNumAnimations; ++a) {
			aiAnimation const *const animation(scene->mAnimations[a]);
			aiNodeAnim  const *const node_anim(find_node_anim(animation, node_name));
			if (!node_anim) continue; // no animation for this node
			model_anim_t::anim_data_t& A(model_anim.animations[a].anim_data[node_name]);
			A.init(node_anim->mNumPositionKeys, node_anim->mNumRotationKeys, node_anim->mNumScalingKeys);
			// position
			for (unsigned i = 0; i < node_anim->mNumPositionKeys; i++) {
				A.pos.emplace_back(node_anim->mPositionKeys[i].mTime, aiVector3D_to_vector3d(node_anim->mPositionKeys[i].mValue));
			}
			// rotation
			for (unsigned i = 0; i < node_anim->mNumRotationKeys; i++) {
				A.rot.emplace_back(node_anim->mRotationKeys[i].mTime, aiQuaternion_to_glm_quat(node_anim->mRotationKeys[i].mValue));
			}
			// scaling
			for (unsigned i = 0; i < node_anim->mNumScalingKeys; i++) {
				A.scale.emplace_back(node_anim->mScalingKeys[i].mTime, aiVector3D_to_vector3d(node_anim->mScalingKeys[i].mValue));
			}
		} // for a
		for (unsigned i = 0; i < node->mNumChildren; i++) {
			unsigned const child_ix(extract_animation_data_recur(scene, node->mChildren[i], model_anim));
			model_anim.anim_nodes[node_ix].children.push_back(child_ix);
		}
		return node_ix;
	}
	void extract_animation_data(aiScene const *const scene, model_anim_t &model_anim) {
		assert(scene && scene->mRootNode);
		model_anim.root_transform = cur_xf.create_xform_matrix();
		model_anim.global_inverse_transform = aiMatrix4x4_to_xform_matrix(scene->mRootNode->mTransformation).inverse();
		model_anim.animations.resize(scene->mNumAnimations);

		for (unsigned a = 0; a < scene->mNumAnimations; ++a) {
			if (scene->mAnimations[a]->mTicksPerSecond) {model_anim.animations[a].ticks_per_sec = scene->mAnimations[a]->mTicksPerSecond;} // defaults to 25
			model_anim.animations[a].duration = scene->mAnimations[a]->mDuration;
		}
		extract_animation_data_recur(scene, scene->mRootNode, model_anim);
	}

	void parse_single_bone(int bone_index, aiBone const *const pBone, mesh_bone_data_t &bone_data, model_anim_t &model_anim, unsigned first_vertex_offset) {
		unsigned const bone_id(model_anim.get_bone_id(pBone->mName.C_Str()));

		if (bone_id == model_anim.bone_transforms.size()) { // maybe add a new bone
			model_anim.bone_transforms.push_back(xform_matrix());
			model_anim.bone_offset_matrices.push_back(aiMatrix4x4_to_xform_matrix(pBone->mOffsetMatrix));
		}
		for (unsigned i = 0; i < pBone->mNumWeights; i++) {
			aiVertexWeight const &vw(pBone->mWeights[i]);
			unsigned const vertex_id(first_vertex_offset + vw.mVertexId);
			assert(vertex_id < bone_data.vertex_to_bones.size());
			bone_data.vertex_to_bones[vertex_id].add(bone_id, vw.mWeight, had_vertex_error);
		}
	}
	void parse_mesh_bones(aiMesh const *const mesh, mesh_bone_data_t &bone_data, model_anim_t &model_anim, unsigned first_vertex_offset) {
		for (unsigned int i = 0; i < mesh->mNumBones; i++) {parse_single_bone(i, mesh->mBones[i], bone_data, model_anim, first_vertex_offset);}
	}
	void process_mesh(aiMesh const *const mesh, aiScene const *const scene, model_anim_t &model_anim) {
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

			if (!load_animations) { // not legal to apply model transform here; must be applied after bone transforms
				cur_xf.xform_pos   (v.v);
				cur_xf.xform_pos_rm(v.n);
			}
			if (mesh->mTextureCoords != nullptr && mesh->mTextureCoords[0] != nullptr) { // TCs are optional and default to (0,0); we only use the first of 8
				v.t[0] = mesh->mTextureCoords[0][i].x; 
				v.t[1] = mesh->mTextureCoords[0][i].y;
			}
			point bcube_pt(v.v);
			if (load_animations) {cur_xf.xform_pos(bcube_pt);} // if we didn't transform the point above, transform it now to compute a (hopefully more accurate) bcube
			if (i == 0) {mesh_bcube.set_from_point(bcube_pt);} else {mesh_bcube.union_with_pt(bcube_pt);}
		} // for i
		assert(mesh->mFaces != nullptr);
		assert(mesh->mNumFaces > 0); // if there were verts, there must be faces

		for (unsigned i = 0; i < mesh->mNumFaces; i++) { // process faces/indices
			aiFace const& face(mesh->mFaces[i]);
			assert(face.mNumIndices == 3); // must be triangles
			for (unsigned j = 0; j < face.mNumIndices; j++) {indices.push_back(face.mIndices[j]);}
		}
		if (!mesh_bcube.is_all_zeros()) {
			if (load_animations) {} // TODO: need to apply model transform to mesh_bcube
			model.union_bcube_with(mesh_bcube);
		}
		//if (mesh->mMaterialIndex >= 0) {} // according to the tutorial, this check should be done; but mMaterialIndex is unsigned, so it can't fail?
		material_t &mat(model.get_material(mesh->mMaterialIndex, 1)); // alloc_if_needed=1
		bool const is_new_mat(mat.empty());
		unsigned const first_vertex_offset(mat.add_triangles(verts, indices, 1)); // add_new_block=1; should return 0
		//cout << TXT(mesh->mName.C_Str()) << TXT(mesh->mNumVertices) << TXT(mesh->mNumFaces) << TXT(mesh->mNumBones) << endl;
		
		if (load_animations && mesh->HasBones()) { // handle bones
			mesh_bone_data_t &bone_data(mat.get_bone_data_for_last_added_tri_mesh());
			bone_data.vertex_to_bones.resize(first_vertex_offset + mesh->mNumVertices);
			parse_mesh_bones(mesh, bone_data, model_anim, first_vertex_offset);
			for (unsigned i = first_vertex_offset; i < bone_data.vertex_to_bones.size(); ++i) {bone_data.vertex_to_bones[i].normalize();} // normalize weights to 1.0
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
			// Note: The version of assimp I have installed in Ubuntu doesn't have AI_MATKEY_TRANSMISSION_FACTOR
			aiGetMaterialFloatArray(material, AI_MATKEY_TRANSPARENCYFACTOR, &mat.tr, &max4);
			//if (aiGetMaterialIntegerArray(mtl, AI_MATKEY_ENABLE_WIREFRAME, &wireframe, &max) == AI_SUCCESS) {}
			//if (aiGetMaterialIntegerArray(mtl, AI_MATKEY_TWOSIDED,         &two_sided, &max) == AI_SUCCESS) {}
			// illum? tf?
		}
	}  
	void process_node_recur(aiNode const *const node, aiScene const *const scene, model_anim_t &model_anim) {
		assert(node != nullptr);
		//print_assimp_matrix(node->mTransformation);
		// process all the node's meshes (if any), in tree order rather than simply iterating over mMeshes
		for (unsigned i = 0; i < node->mNumMeshes; i++) {process_mesh(scene->mMeshes[node->mMeshes[i]], scene, model_anim);}
		// then do the same for each of its children
		for (unsigned i = 0; i < node->mNumChildren; i++) {process_node_recur(node->mChildren[i], scene, model_anim);}
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
		unsigned flags(aiProcess_Triangulate | aiProcess_SortByPType | aiProcess_JoinIdenticalVertices |
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
		process_node_recur(scene->mRootNode, scene, model.model_anim_data);
		if (load_animations) {extract_animation_data(scene, model.model_anim_data);}
		model.finalize(); // optimize vertices, remove excess capacity, compute bounding sphere, subdivide, compute LOD blocks
		model.load_all_used_tids();
		if (verbose) {cout << "bcube: " << model.get_bcube().str() << endl << "model stats: "; model.show_stats();}
		return 1;
	}
};

bool read_assimp_model(string const &filename, model3d &model, geom_xform_t const &xf, int recalc_normals, bool verbose) {
	//timer_t timer("Read AssImp Model");
	bool const load_animations = 1;
	file_reader_assimp reader(model, load_animations);
	return reader.read(filename, xf, recalc_normals, verbose);
}

#else // ENABLE_ASSIMP

bool read_assimp_model(string const &filename, model3d &model, geom_xform_t const &xf, int recalc_normals, bool verbose) {
	cerr << "Error: Assimp model import has not been enabled at compile time" << endl;
	return 0;
}

#endif
