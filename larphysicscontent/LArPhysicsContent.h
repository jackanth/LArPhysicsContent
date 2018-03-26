/**
 *  @file   larphysicscontent/LArPhysicsContent.h
 *
 *  @brief  Header file detailing content for use with physics analyses.
 *
 *  $Log: $
 */
#ifndef LAR_PHYSICS_CONTENT_H
#define LAR_PHYSICS_CONTENT_H 1

namespace pandora { class Pandora; }

/**
 *  @brief  LArPhysicsContent class
 */
class LArPhysicsContent
{
public:
    /**
     *  @brief  Register all the lar content algorithms and tools with pandora
     *
     *  @param  pandora the pandora instance with which to register content
     */
    static pandora::StatusCode RegisterAlgorithms(const pandora::Pandora &pandora);
};

#endif // #ifndef LAR_PHYSICS_CONTENT_H
