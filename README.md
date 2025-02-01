![REPO SIZE](https://img.shields.io/github/repo-size/fegennari/3DWorld.svg) 
![CODE SIZE](https://img.shields.io/github/languages/code-size/fegennari/3DWorld.svg) 
![License](https://img.shields.io/github/license/fegennari/3DWorld.svg)

3DWorld Blog: https://3dworldgen.blogspot.com

3DWorld is a cross-platform OpenGL-based 3D Game Engine that I've been working on since I took the CS184 computer graphics course at UC Berkeley in 2001.
It has the following features:
* 3D graphics functions, classes, and wrappers around OpenGL
* Shader generator/processor with hot reload
* Procedural content generation for terrain, vegetation, cities, building interiors and exteriors, etc.
* Procedural universe generator with galaxies, stars, planets, moons, etc.
* Procedural voxel 3D terrain generation with realtime user editing
* Terrain generator including various noise functions, erosion, realtime user editing, heightmap read/write
* Procedural building (interior and exterior), road, and city generation
* Physics simulation for primitive object types and others (> 10K dynamic objects)
* Realtime day/night cycle with weather (rain, snow, hail, wind, lightning)
* Physically based materials with reflection and refraction
* Dynamic shadows, ambient occlusion, up to 1024 dynamic light sources, postprocessing effects
* Skeletal animation and procedural animation
* Built-in first person shooter game "smiley killer"
* Build-in spaceship + planet colonization game
* Building open world item collection and zombie gameplay mode
* Computer AI for players in the FPS game, building gameplay, and ships in the universe game
* Importer for Lightwave object file, 3DS formats, and Assimp for other file formats
* Reading support for textures: JPEG, PNG, BMP, TIFF, TGA, DDS
* Optimized for fast load and realtime rendering of large models (> 1GB of vertex/texture data)

I converted the project from svn to git at commit 6607.
Most of the code is written in C++, with GLSL for shaders.
This is intended to be a cross-platform project.
Microsoft Visual Studio 2022 project files are included.
A linux/gcc makefile is also included, but is more experimental. See README.linux for more details.
The project should build under gcc on linux with some work, but it's been a while since I tried this.
I have an old makefile that is out of date, but may not take too much work to fixup and make it usable.

Be warned, this is a large repository, currently about 1GB.
I've included source code, config files, textures, sounds, small models, lighting files, scene data, heightmaps, and project files.
This repo does not contain the large model files used in some scenes, you'll have to download these separately.
This means that some of the scene config files won't work because they can't find their referenced data.
The current list of dependencies is:
* OpenGL 4.5 (Should come with Windows 8/10/11 latest graphics drivers)
* OpenAL 1.1 (optional) (System Install: https://www.openal.org/downloads/ or you can try the newer openal-soft: https://github.com/kcat/openal-soft)
* freeglut-2.8.1 (Current 3.0 version probably works: https://sourceforge.net/projects/freeglut/)
* freealut-1.1.0 (optional) (One version is here: https://github.com/vancegroup/freealut)
* zlib-1.2.11 (You can download a newer version from here: https://zlib.net/)
* glew-2.0.0 (2.1.0 probably works as well: http://glew.sourceforge.net/)
* gli (Latest version: https://github.com/g-truc/gli / header only, included in dependencies directory)
* glm-0.9.9.0 (Latest version: https://glm.g-truc.net/0.9.9/index.html or https://github.com/g-truc/glm / header only, included in dependencies directory)
* libpng-1.2.20 (optional) (My version is very old; Latest version: https://libpng.sourceforge.io/index.html); Can be replaced with stb_image in most cases, except for map image export.)
* libtiff-4.3.0 (optional) (Latest version: http://www.simplesystems.org/libtiff/)
* Assimp (optional) (See build instructions at https://github.com/assimp/assimp/blob/master/Build.md ; The vcpkg build instructions are probably the easiest.)
* libtarga (source included)
* STB headers: stb_image, stb_image_write, stb_dxt (source included)

I've included stripped down versions of most of these libraries in the dependencies directory.
I removed all large files that aren't required by 3DWorld, in some cases even examples/tests/documentation.
These have been built with MS Visual Studio 2022 Community on Windows 11.
If you want to use these, you'll need to copy the directories to the root directory and rebuild any libraries needed for other versions of Windows or Visual Studio.
If you clone/install vcpkg it should be at the same level as the 3DWorld directory.

Note that many of these dependencies are old and could be replaced with newer libraries. I've been concentrating on adding content and I'm not too interested in this.
Freeglut should probably be replaced with SDL, and the image libraries with STB or DevIL. (STB is used as a fallback but doesn't support all of the images used.)

If you want to build 3DWorld, you can use the projects in the dependencies/ folder, or download and build them yourself and change the project settings to use them.
I currently use the x64 MS Visual Studio 2022 Community build target for 3DWorld, but the win32 build target also works.
The MSVS 2019 project 3DWorld_msvs2019.vcxproj is currently out of date but can possibly be made to work.
It should compile and run in 32-bit mode if you copy the DLLs from the lib64/ folder into the root of the repo and make some other project settings changes.

If you have linux, you can try to build using the provided makefile. The file README.linux should be helpful.
I've gotten 3DWorld to build and mostly run on Ubuntu 18.04 with gcc 7 and Ubuntu 20.04 with gcc 9.

3DWorld takes a config filename on the command line. If not found, it reads defaults.txt and uses any config file(s) listed there.
Some of these congig files include models such as the Sponza Atrium, Stanford Dragon, sportscar, etc.
These files are too large to store in the git repo. I've attempted to have 3DWorld generate nonfatal errors if the models can't be found.
Many of the larger models can be found at the McGuire Computer Graphics Archive:
http://casual-effects.com/data/

I've packaged up the 3D models that are too large for the GitHub repo and put them on Google Drive here (up to v6 now):
https://drive.google.com/file/d/1crN9rqT-LSvYyTZTw5wtkhwsubE651ex/view?usp=sharing
Some of these models are stored in 3DWorld's internal format and should not be reused in other projects. Others come from websites such as Mixamo.
There is also a textures directory with additional textures used with building interiors that can be merged with the project textures directory.

System requirements:
* Windows 8/10/11; Linux when using the makefile with gcc.
* Microsoft Visual Studio 2019 or 2022. The professional or community version is needed for OpenMP support. You can also try to use gcc on linux.
* A relatively new generation of Nvidia or AMD GPU (Runs on my laptop with Intel graphics, but at 12-20 FPS)
* At least 8GB system memory for the larger scenes
* At least 4GB GPU memory for the larger scenes; My GPU has 12GB of memory

Troubleshooting:
It seems like some systems (AMD cards in particular) require an OpenGL core context. This can be selected by adding "use_core_context 1" in the config file.
This can also be enabled in scene_config/config_post.txt, which is a file that applies after reading all other top-level config files.
In some situations (some Nvidia cards), using a core context can be slower, which is why I don't have it enabled by default.

Useful Keys (see readme-keys.txt for more key bindings):
* a,s,d,w: Movement
* q,e: Change weapon (gameplay mode)
* 'space': Jump/fire
* 'esc': Quit
* 'tab': Onscreen menu (navigate with arrow keys and 'X' to switch menus)
* b: Enable objects/physics/AI (for ground mode and universe mode gameplay)
* F1: Switch between ground/universe/tiled terrain modes
* F2: Toggle gameplay mode
* m: Toggle fullscreen mode
* h: Toggle flight move
* v: Change camera mode (crystal ball/orbit vs. first person)
* V: Toggle mouse look
* K: Toggle overhead map mode
* x: Pause
* Mouse Left: Turn/Action
* Mouse Right: Fire

I currently have this repo up for educational purposes under the GPLv3 license.
Some sub-modules are available with other licenses compatible with commercial use in my GitHub account.
It's not meant as a commercial tool and I'm not trying to make money here.
I'm also not looking for others to work on the project at this stage, though I'm accepting feedback, bug reports, and suggestions.
Maybe things will change if I decide to make a real game out of this.
If you would like to use something here for your project, please let me know.

There is no further documentation for 3DWorld.
However, I do have a blog that includes descriptions of the algorithms and lots of screenshots:
https://3dworldgen.blogspot.com

Here are some screenshots linked from my blog:

![alt text](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEh7bt1KabOKD8tUIxaACMEEDsk2oP4105fUaz4RMGqucAxEAkSKGN4-sX0yMv2IXHMUkWNUNspp0KDuGZIduZoBXAo_K_PQizCx_CtfGjj88GLMXw44jOGmPAs2VUwy__8GckC5MICg4dUV0cK-ripMByVsHX80Wjkzb4XRuuugZXi5_nT4Rc5VwTTI6nP7/w640-h360/mall_7_levels.jpg)

This is the interior of a 7 level procedurally generated mall interior. (config_heightmap.txt)

![alt text](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEh7PsshxXhRFU_pkkJn9yv4e4-STzjeO1ZnwlKEXL2xGC-nrf0mb9fL4G2etkCRonHcjb7OQ1IdLrNhA9lDyTmitVcDQZeXVBRTsww2FxhEwbCwGVZ4XdthBSYfZL4EM-SsRx6yfBjYDozqpEaje0vDy5pdryYG2Flw45Dj4sFRPyvejizgTxerSYWJUQJU/w640-h360/mall_food_court.jpg)

Here is a procedurally generated mall food court with tables and chairs, and stores to the sides. (config_heightmap.txt)

![alt text](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEiDYcHMxI8d4EmqEt_R9zwbcmtk0wfCadahSXeP3xZU9SdXGkMqMAcFVxMMlcOvDRh7Lg5KhEjgosmId9HtHdtuErQd0okK7Vb7IU0RHeKpONDTc887BDIgWeuR8E0dHChU-vKQxkdnjWzaRLGSWBjmIa3O7z6PYR9v-j4oNYpuOzD9PJcv4XR-45YJYghK/w640-h360/racks1.jpg)

Here is a retail area of a building with racks of shelves. Buildings can contain tens of thousands of interactive objects the player can take and carry around. (config_heightmap.txt)

![alt text](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEgI1KwsSM00BSkCWUr5OPrXMg5lsjeKk7GOPFdmGkk1NsHlYv4w-StMW8NJKZMHpRt3EoriqCyjD0xRM2E2BqLSBX5zts36VHn8s6TEuQCuhHFq6sFzLyqWSdOxroqB_37utAtB8nN0soiLp-cVFJwnL8VHyBdNNupxdF3HbjmtDdPm4-clTvTquk_JgOYK/w640-h360/sec_monitors_no_base.jpg)

Cameras placed in other rooms show realtime video feeds on security monitors. (config_heightmap.txt)

![alt text](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEgUaakhysxiViFwntixF7D2F_FAOnU6mnm6A6jIkpdVrLmu5aW80dFXJ6pTiWWUcoamsLy44aPl3aZtM5QgwLDQxk2LpoGvy03TMAW0jlOr7bpNsoI-en8TXIJDjxkGhgIxhq1yFvRfq5xeYX1ElEK_kIUMt-Fsw9SaT_ap7tpNVXwrxLpJ5b9nCB-rruWz/w640-h338/zombie_pool_party.jpg)

Zombies (with path finding and other AI behaviors) walking around and falling into an indoor pool. The water water surface has ripples, splashes, reflections, and refractions. (config_heightmap.txt)

![alt text](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEjtYUsFpkEeENoYmiyj0tWBJY-vAOfBhisVfK6GRhA5IZlprB0MgQ9KgyxecHigRF9gm5kiR4ADcH0LvxUdFrOtV1YsUYfeF1eYmWJF76vvUa-_LjyDcQlo80Q9ReNcfeSSLJ-kRMVGILU3zMBBXla_8ph0fcRV8DKJJwU1NdH0s5tgLPLDQ5wFe-JikA/w640-h360/sprinklers_with_cars_and_indir.jpg)

Parking garage with parked cars and sprinkler pipes routed on the ceiling. (config_heightmap.txt)

![alt text](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEgbKNuwYqJV9LSFC7hg18_yb76dG9JecG74WqLQs7ZMiS0sXTyJO6k7RbDpoQWH09-gh3eJrN4jabU64DvspqLeZvT1Uvfp87dbF2CKsK-KRDtsE4skVYHW2dcWC1POaIsy3ANCZj-BONgOFuIJtqCP1GcFCLSjt827q8OUnAnbKu_KTN4_iaBlEGQuWA/s1920/basement_indir_pipes.jpg)

Procedural building basement with people, indirect lighting, and pipes routed along the ceiling. (config_heightmap.txt)

![alt text](https://1.bp.blogspot.com/-ZbJUmGiha84/YTRRfRWTE4I/AAAAAAAADEs/yB9tPZcllnM40FOsh8nkut3HtDjm0MLqQCLcBGAsYHQ/s1920/residential_grid.jpg)

Procedural residential neighborhood with office buildings in the backround. There are cars and people both on the sidewalks and inside buildings. (config_heightmap.txt)

![alt text](https://1.bp.blogspot.com/-GW82PSnZt7s/X7DVx4wb4aI/AAAAAAAACxM/4PSV1e2iI8wzVEXVnA9K3GPrPqfdpKmcwCLcBGAsYHQ/w640-h360/office_libraries.jpg)

Procedural building interior - a library with procedurally generated books. (config_heightmap.txt)

![alt text](https://1.bp.blogspot.com/-x8EY3wBRA2c/X-2bpKYDQ5I/AAAAAAAAC1g/4SAgBAdGIzY8LmDFgCLjyIMsrwsBsbplgCLcBGAsYHQ/w640-h360/solar_panels_with_sides.jpg)

Procedurally generated houses and office buildings with full interiors. (config_heightmap.txt)

![alt text](https://1.bp.blogspot.com/-sGemMLy_VXc/Wx10wAEwSxI/AAAAAAAABj8/FKLHXZ27qBsNPXE7fZU1vXo-N06NjV5EQCLcBGAs/s640/spheres_low_noise.jpg)

Reflective (metal), refractive (glass), emissive, and translucent spheres in night time scene with indirect lighting. Drawn in realtime; spheres can be moved interactively. (config_white_plane.txt)

![alt text](https://2.bp.blogspot.com/-f07b_YCw-7Q/XFd5v1jazTI/AAAAAAAAByk/r1xQ0zgTnmUT-ONsv4y7W1X9LwEtQdJuACLcBGAs/s640/peds_waiting.jpg)

Procedural city with buildings, roads, cars, and pedestrians. (config_heightmap.txt)

![alt text](https://1.bp.blogspot.com/-acI3Ly40-Hk/WzSWMckhOiI/AAAAAAAABlg/KvzdEJ9qEjUJPOF7kYvh1RpELBSnnQXtgCEwYBhgL/s640/bridge_night1.jpg)

Procedural city at night with bridge in the foreground. (config_heightmap.txt)

![alt text](https://1.bp.blogspot.com/-H3QY3vua23s/WovTmQk8I6I/AAAAAAAABRE/KSFYKSdDRAAVPA7NcZEQXCpKJJaUY9hWQCLcBGAs/s640/connected_cities_trees.jpg)

Early procedural city with pine trees.

![alt text](https://1.bp.blogspot.com/-4PWdGiTpsfw/WrdOCFbjxnI/AAAAAAAABWo/O0pT-TfRxMMaCNDnxs0UdSayEpU3y-XWACLcBGAs/s640/city_cars3.jpg)

Early procedural city.

![alt text](https://1.bp.blogspot.com/-orCzK6w5xEM/WjmW-0kXVXI/AAAAAAAABN8/Sa73QhiUnvwC0CqC_ZSAPnT1KX8miu85gCLcBGAs/s640/erosion_from_above.jpg)

Terrain using domain warp noise with hydraulic erosion simulation before city has been placed. (config_heightmap.txt)

![alt text](https://4.bp.blogspot.com/-_yqljuQRYFA/WjN0WlQXxyI/AAAAAAAABMo/VLQbV8HF9nMlh0yoqUS57vxr6uer2RrswCEwYBhgL/s640/river.jpg)

Tiled terrain mode with river, trees, grass, etc. (config_t.txt)

![alt text](https://4.bp.blogspot.com/-ZqmYa0act0w/We1_2z6l1VI/AAAAAAAABKU/uXNawQ9xwnAqD0E8pCdz7MouyXVYEdczgCEwYBhgL/s640/fires3.jpg)

Realtime interactive destructive fire/smoke simulation involving trees, plants, and grass. (config_trees.txt)

![alt text](https://2.bp.blogspot.com/-h9eUV4FiM28/WZKe4pr7YpI/AAAAAAAABGE/pXBPNL0OJi48ErNDS6RH0IprW7V_W5XtACLcBGAs/s640/nebula_rings_asteroids.jpg)

Procedural universe solar system with asteroid belt. (universe/config_universe.txt)

![alt text](https://4.bp.blogspot.com/-5dr80n928lw/WT9t81moCVI/AAAAAAAAA_s/YldLPL_Y__gB4q4QtrVrAflmhl2X17qEgCLcB/s640/sponza2.jpg)

Crytek Sponza atrium scene with dynamic shadow casting light source and indirect lighting. (sponza/config_sponza2.txt)

![alt text](https://4.bp.blogspot.com/-kY8qCSsE0ck/WQbOLeCBVtI/AAAAAAAAA9E/SYUEjT1YGEgllNZiaf-bU3JWg5lta0pNACEw/s640/museums_120.jpg)

10,000 instances of a highly detailed museum model placed in tiled terrain mode (Puget Sound heightmap) and drawn in realtime with shadows and indirect lighting. (config_museum_tt_model.txt)

![alt text](https://2.bp.blogspot.com/-MSjc6z9NjRc/WQYkxjOEOuI/AAAAAAAAA74/mrXX01ljwZMkkx1kwlGP4sztMKv7dMTmwCEw/s640/sponza.jpg)

Crytek Sponza atrium with reflective floors and 200 dynamic point light sources in realtime. (sponza/config_sponza2.txt)

![alt text](https://3.bp.blogspot.com/-Ys37EWGm-PU/WKv1nh2y6tI/AAAAAAAAA5A/X4zcAp2f-Y8UD5vIT7-n7AJMwTUjZQSBACEw/s640/many_objects.jpg)

Many reflective/refractive spheres and cubes with density-based light attenuation in the office building courtyard. (mapx/config_mapx.txt)

![alt text](https://4.bp.blogspot.com/-TktxFf1hZ_o/WI14nt33o_I/AAAAAAAAA3Y/DhkJmBlRzuECI4lRdBGDgrikGlarUA_WQCLcB/s640/snow_mask_hr.jpg)

San Miguel scene with dynamic snowfall and path traced snow coverage map rendered with a custom shader. (config_san_miguel.txt)

![alt text](https://1.bp.blogspot.com/-g8X-SATTjbM/WHKTPJ3EURI/AAAAAAAAA2E/18Nf3oK0ZkUoPh_H6_VBzqssbsjdIhLKwCLcB/s640/bright1.jpg)

San Miguel scene rendered in realtime with path traced precomputed indirect lighting, normal maps, and cube mapped reflective surfaces. (config_san_miguel.txt)

![alt text](https://4.bp.blogspot.com/-uxrpc3f1IEY/WCgqzNFqCvI/AAAAAAAAAzs/6lmvfaEJ_dU5MSPpdZgwVyo2EgmIWH0kgCLcB/s640/all_pine_trees.jpg)

2M procedurally generated + placed pine trees (500K visible) drawn in realtime in tiled terrain mode using instanced billboards. (config_t.txt)

![alt text](https://4.bp.blogspot.com/-yP383fqlaRk/V_8JfWkjyWI/AAAAAAAAAxY/2eH5WnWktgwyFUwhRRYGAnT9trfjkFRswCEw/s640/cubes_10k.jpg)

Interactive stacking of 10,000 dynamic boxes in mapx scene. Stacks can be moved and knocked over, and scene can be saved/reloaded/edited. (mapx/config_mapx.txt)

![alt text](https://3.bp.blogspot.com/-8CIw4xIUMdk/V8KKax0v7GI/AAAAAAAAAvc/_PtkZvOaY5IiptMUGgKa4fVKsJKs1kOWwCLcB/s640/ringed_planet.jpg)

Procedurally generated planet with rings containing small asteroids. (universe/config_universe.txt)

![alt text](https://3.bp.blogspot.com/-jWcp7MUFoeU/V4czYLPtHVI/AAAAAAAAAsw/B8x0U9QaGiEnKNcDDoEQRcy4WYCbIermQCEw/s640/asteroid_belt_inside.jpg)

Asteroid belt with 10,000 dynamic (rotating and orbiting) asteroids, 100,000 point sprite particles, and volume billboard clouds for dust. (universe/config_universe.txt)

![alt text](https://2.bp.blogspot.com/-lEunlK-ZyT8/V1EkIeSttaI/AAAAAAAAArQ/X7G170MT_6IfDd74WSyY_biN_dIQEzCsQCLcB/s640/trees_above.jpg)

Tiled terrain mode, showing various forms of procedurally generated vegetation, including grass, flowers, pine trees, palm trees, and deciduous trees. These are actually polygon models that cast shadows. (config_t.txt)

![alt text](https://3.bp.blogspot.com/-1kFHDCXg65Y/VzgfWd1phkI/AAAAAAAAAqA/7rrUIxsVH-Ea2p_NkbRYuyeUZuuj0cejgCKgB/s640/vfog_smoke_rays.jpg)

Realtime, dynamic light shafts through noise-based volumetric smoke in the Crytek Sponza scene. (sponza/config_sponza2.txt)

![alt text](https://4.bp.blogspot.com/-mk9EA7i_3LU/VuoVYWaz4oI/AAAAAAAAAnw/DHnb2yr_v0scME_D8DRdIsiyl6a8AMOyA/s640/sponza_metal_floor3.jpg)

Reflective, metallic floors in the Crytek Sponza scene, with indirect sun + sky lighting. (sponza/config_sponza2.txt)

![alt text](https://3.bp.blogspot.com/-lCw93zflw0k/Vt6EOqYZUZI/AAAAAAAAAlY/yvfndWb5Okw/s1600/reflective_objects_front.jpg)

Reflective spheres, torus, and cube drawn in realtime in a dynamic scene using cube map reflection textures in the office building lobby. Materials are a mix of metals and dielectrics using physically-based models. (mapx/config_mapx.txt)

![alt text](https://4.bp.blogspot.com/-TueDAN3BGEw/VsbLtNP0dwI/AAAAAAAAAjg/-kBjm4EQDd0/s640/museum2.jpg)

Museum scene with indirect lighting, reflective surfaces, and shadow mapping. This model is procedurally textured. (config_museum.txt)

![alt text](https://1.bp.blogspot.com/-2AzAKVCUhvw/VpGvWG6uQwI/AAAAAAAAAgM/3QLnzeiaeCw/s1600/snow_scene.jpg)

Snowy house scene generated by dropping a billion snow particles and accumlating snow. Snow is precomputed but can be rendered in realtime. (house/config_house_winter.txt)

![alt text](https://1.bp.blogspot.com/-2LlXIzcVDnA/VmUt8R4cwGI/AAAAAAAAAeY/Mx_xy30eVCQ/s640/house_rain.jpg)

Rain simulation and rendering, including collision detection for each raindrop and wet/reflective surfaces. (house/config_house.txt)

![alt text](https://3.bp.blogspot.com/-PaceyCZ6W6M/VcmIQAHy8gI/AAAAAAAAARU/DnDUnSyewVs/s1600/dir_plus_indir.jpg)

Indirect lighting and shadows from overlapping spotlights applied to dynamic objects. (mapx/config_mapx.txt)

![alt text](https://3.bp.blogspot.com/-kzdVvau1pM4/VP6DGtZZg9I/AAAAAAAAALM/xxxIH2H_byM/s1600/waves.jpg)

Water simulation with waves using hardware tessellation in tiled terrain mode. (config_t.txt)

![alt text](https://4.bp.blogspot.com/-iQt6mdroDTY/VMQIqQ2PiJI/AAAAAAAAAJE/0TVtUzsvdqM/s1600/voxel_snow_ao.jpg)

Procedurally generated 3D voxel ice caves with indirect lighting. This terrain can be edited in realtime with brushes and weapon fire. (config_ice_caves.txt)

![alt text](https://3.bp.blogspot.com/-nyVNDhCQKvo/VKnl-TOYsrI/AAAAAAAAAGs/UDBBNEQGRBk/s1600/terran_planet.jpg)

Procedurally generated planet drawn entirely in the fragment shader with no textures. This includes terrain generation (ice, snow, rock, forest, desert, dirt), water, clouds, and atmosphere. (universe/config_universe.txt)
