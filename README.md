3DWorld is an OpenGL-based 3D Game engine that I've been working on since I took a computer graphics course at UC Berkeley in 2001.
It has the following features:
* 3D graphics functions, classes, and wrappers around OpenGL
* Shader generator/processor with hot reload
* Procedural content generation for terrain, vegetation, buildings, etc.
* Procedural universe generator with galaxies, stars, planets, moons, etc.
* Procedural voxel 3D terrain generation with realtime user editing
* Terrain generator including various noise functions, erosion, realtime user editing, heightmap read/write
* Physics simulation for primitive object types and others
* Built-in first person shooter game "smiley killer"
* Build-in spaceship + planet colonization game
* Computer AI for players in the FPS game and ships in the universe game
* Importer for Lightwave object file and 3DS formats
* Reading support for textures: JPEG, PNG, BMP, TIFF, TGA, RAW, DDS
* Optimized for fast load and realtime rendering of large models (> 1GB of vertex/texture data)

I converted the project from svn to git at commit 6607.
Most of the code is written in C++, with GLSL for shaders.
This is intended to be a cross-platform project.
Microsoft Visual Studio 2015 project files are included.
The project should build under gcc on linux with some work, but it's been a while since I tried this.
I have an old makefile that is out of date, but may not take too much work to fixup and make it usable.

Be warned, this is a large repository, currently 770MB.
I've included source code, config files, textures, sounds, small models, lighting files, scene data, heightmaps, and project files.
This repo does not contain the dependencies or large model files, you'll have to download these separately.
This means that some of the scene config files won't work because they can't find their referenced data.
The current list of dependencies is:
* OpenGL 4.4 (Should come with Windows 7/8/10 latest graphics drivers)
* OpenAL 1.1 (System Install: https://www.openal.org/downloads/)
* freeglut-2.8.1 (Current 3.0 version probably works: https://sourceforge.net/projects/freeglut/)
* freealut-1.1.0 (One version is here: https://github.com/vancegroup/freealut)
* zlib-1.2.1 (You can download a newer version from here: https://zlib.net/)
* glew-2.0.0 (2.1.0 probably works as well: http://glew.sourceforge.net/)
* gli-0.5.1.0 (Latest version: https://github.com/g-truc/gli)
* glm-0.9.5.2 (Latest version: https://glm.g-truc.net/0.9.8/index.html My version: https://glm.g-truc.net/0.9.5/index.html)
* libjpeg-9a (My version is old; Latest version: http://www.ijg.org/)
* libpng-1.2.20 (My version is very old; Latest version: https://libpng.sourceforge.io/index.html)
* libtiff-4.0.3 (Latest version: http://www.simplesystems.org/libtiff/)
* libtarga (source included)

Note that many of these dependencies are old and could be replaced with newer libraries. I've been concentrating on adding content and I'm not too interested in this.
Freeglut should probably be replaced with SDL, the last 4 image libraries with DevIL, and maybe assimp can be used for model loading.

If you want to build 3DWorld, you'll need to download and build these dependencies somewhere and change the project settings to use them.
I just copy these into the current directory and have these files ignored by git/svn.
I currently use a 32-bit MS Visual Studio build target for 3DWorld.
It should compile in 64-bit mode, but I couldn't find compatible 64-bit debug libraries for OpenAL,
and a few of the other dependencies didn't build cleanly in 64-bit mode.

3DWorld takes a config filename on the command line. If not found, it reads defaults.txt and uses any config file(s) listed there.
Some of these congig files include models such as the Sponza Atrium, Stanford Dragon, sportscar, etc.
These files are too large to store in the git repo.
Many of the larger models can be found at the McGuire Computer Graphics Archive:
http://casual-effects.com/data/

System requirements:
* Windows 7/8/10 (Runs on Windows 7, but I've only built on 8 and 10). Linux if you create a makefile for gcc.
* Microsoft Visual Studio 2015 (or newer?). The professional version is needed for OpenMP support. You can also try to use gcc.
* A relatively new generation of Nvidia or ATI GPU (Runs on my laptop with Intel graphics, but at 12-20 FPS)
* At least 4GB system memory for the larger scenes
* At least 2GB GPU memory for the larger scenes

I currently have this repo up for educational purposes under the GPL license.
It's not meant as a commercial tool and I'm not trying to make money here.
I'm also not looking for others to work on the project at this early stage, though I'm accepting feedback and suggestions.
Maybe things will change if I decide to make a real game out of this.
If you would like to use something here for your project, please let me know.

There is no further documentation for 3DWorld.
However, I do have a blog that includes descriptions of the algorithms and lots of screenshots:
https://3dworldgen.blogspot.com/
