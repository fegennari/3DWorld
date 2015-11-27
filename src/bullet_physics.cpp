// 3D World - Bullet Physics Library Wrapper
// by Frank Gennari
// 10/13/15

#include "3DWorld.h"
#include "physics_objects.h"

//#define ENABLE_BULLET_PHYSICS

#ifdef ENABLE_BULLET_PHYSICS
#include "btBulletCollisionCommon.h"
#include "btBulletDynamicsCommon.h"
#include "Bullet3AppSupport/CommonRigidBodySetup.h"

struct PhysicsSetup3DWorld : public CommonRigidBodySetup {
	virtual void initPhysics(GraphicsPhysicsBridge& gfxBridge) {
		createEmptyDynamicsWorld();
	}
};

void test() {
	PhysicsSetup3DWorld world;
}

#endif // ENABLE_BULLET_PHYSICS



