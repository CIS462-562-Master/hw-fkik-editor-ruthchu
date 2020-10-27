#include "aActor.h"

#pragma warning(disable : 4018)



/****************************************************************
*
*    	    Actor functions
*
****************************************************************/

AActor::AActor() 
{
	m_pInternalSkeleton = new ASkeleton();
	m_pSkeleton = m_pInternalSkeleton;

	m_BVHController = new BVHController();
	m_BVHController->setActor(this);

	m_IKController = new IKController();
	m_IKController->setActor(this);

	// code to update additional Actor data goes here
	resetGuide();

}

AActor::AActor(const AActor* actor)
{
	*this = *actor;
}

AActor& AActor::operator = (const AActor& actor)
{
	// Performs a deep copy
	if (&actor == this)
	{
		return *this;
	}
	m_pSkeleton = actor.m_pSkeleton;

	// code to update additional Actor data goes here


	return *this;
}

AActor::~AActor()
{
	 delete m_IKController;
	 delete m_BVHController;
	 delete m_pInternalSkeleton;

}

void AActor::clear()
{
	// looks like it is clearing more times than the number of actors.  as a result, m_pSkeleton is not defined for last case.
	m_pSkeleton->clear();  

	// code to update additional Actor data goes here
}

void AActor::update()
{
	if (!m_pSkeleton->getRootNode() )
		 return; // Nothing loaded
	else m_pSkeleton->update();

	// code to update additional Actor data goes here

}

ASkeleton* AActor::getSkeleton()
{
	return m_pSkeleton;
}

void AActor::setSkeleton(ASkeleton* pExternalSkeleton)
{
	m_pSkeleton = pExternalSkeleton;
}

void AActor::resetSkeleton()
{
	m_pSkeleton = m_pInternalSkeleton;
}

BVHController* AActor::getBVHController()
{
	return m_BVHController;
}

IKController* AActor::getIKController()
{
	return m_IKController;
}

void AActor::updateGuideJoint(vec3 guideTargetPos)
{
	if (!m_pSkeleton->getRootNode()) { return; }

	// TODO: 
	// 1.	Set the global position of the guide joint to the global position of the root joint
	// 2.	Set the y component of the guide position to 0
	// 3.	Set the global rotation of the guide joint towards the guideTarget
	vec3 rootWorld = m_Guide.getLocal2Global() * m_pSkeleton->getRootNode()->getGlobalTranslation();
	rootWorld[1] = 0.0;
	m_Guide.setGlobalTranslation(rootWorld);
	vec3 newTarget = (guideTargetPos - m_Guide.getGlobalTranslation());
	newTarget[1] = 0.f;
	newTarget.Normalize();
	
	vec3 xCol = vec3(0.f, 1.f, 0.f).Cross(newTarget);
	vec3 yCol = vec3(0.f, 1.f, 0.f);

	mat3 orient = mat3(xCol, yCol, newTarget).Transpose();
	m_Guide.setGlobalRotation(orient);
	m_pSkeleton->update();
}

void AActor::solveFootIK(float leftHeight, float rightHeight, bool rotateLeft, bool rotateRight, vec3 leftNormal, vec3 rightNormal)
{
	if (!m_pSkeleton->getRootNode()) { return; }
	AJoint* leftFoot = m_pSkeleton->getJointByID(m_IKController->mLfootID);
	AJoint* rightFoot = m_pSkeleton->getJointByID(m_IKController->mRfootID);

	// TODO: 
	// The normal and the height given are in the world space

	// 1.	Update the local translation of the root based on the left height and the right height
	AJoint* root = m_pSkeleton->getRootNode();
	root->setLocalTranslation(root->getLocalTranslation() + vec3(0.f, std::max(leftHeight, rightHeight), 0.f));
	m_pSkeleton->update();

	// 2.	Update the character with Limb-based IK 
	
	// Rotate Foot
	if (rotateLeft)
	{
		// Update the local orientation of the left foot based on the left normal
		// Want to rotate foot so it is flush agains the angle of the ground
		ATarget tar;
		vec3 modFoot = vec3(leftFoot->getGlobalTranslation()[0], leftHeight, leftFoot->getGlobalTranslation()[2]);
		tar.setGlobalTranslation(modFoot);
		m_IKController->IKSolver_Limb(m_IKController->mLfootID, tar);
		vec3 yVec = leftFoot->getLocalRotation().Transpose()[1];

		vec3 yCol = leftNormal.Normalize();
		vec3 zCol = yVec.Cross(yCol).Normalize();
		vec3 xCol = zCol.Cross(yCol).Normalize();
		
		mat3 orient = mat3(xCol, yCol, zCol).Transpose();
		
		leftFoot->setLocalRotation(orient);
		m_pSkeleton->update();
	}
	if (rotateRight)
	{
		// Update the local orientation of the right foot based on the right normal
		ATarget tar;
		vec3 modFoot = vec3(rightFoot->getGlobalTranslation()[0], rightHeight, rightFoot->getGlobalTranslation()[2]);
		tar.setGlobalTranslation(modFoot);
		m_IKController->IKSolver_Limb(m_IKController->mRfootID, tar);
		vec3 yVec = rightFoot->getLocalRotation().Transpose()[1];

		vec3 yCol = rightNormal.Normalize();
		vec3 zCol = yVec.Cross(yCol).Normalize();
		vec3 xCol = zCol.Cross(yCol).Normalize();

		mat3 orient = mat3(xCol, yCol, zCol).Transpose();

		rightFoot->setLocalRotation(orient);
		m_pSkeleton->update();
	}
}
