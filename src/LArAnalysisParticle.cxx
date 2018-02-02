/**
 *  @file   LArPhysicsContent/src/LArAnalysisParticle.cc
 *
 *  @brief  Implementation of the lar analysis particle class.
 *
 *  $Log: $
 */

#include "LArAnalysisParticle.h"

using namespace pandora;

namespace lar_physics_content
{

LArAnalysisParticle::LArAnalysisParticle(const LArAnalysisParticleParameters &parameters) :
    ParticleFlowObject(parameters),
    m_type(parameters.m_type),
    m_typeTree(parameters.m_typeTree),
    m_analysisEnergy(parameters.m_analysisEnergy),
    m_energyFromCharge(parameters.m_energyFromCharge),
    m_isVertexFiducial(parameters.m_isVertexFiducial),
    m_fiducialHitFraction(parameters.m_fiducialHitFraction),
    m_vertexPosition(parameters.m_vertexPosition),
    m_directionCosines(parameters.m_directionCosines),
    m_analysisMomentum(parameters.m_analysisMomentum),
    m_numberOf3dHits(parameters.m_numberOf3dHits),
    m_numberOfCollectionPlaneHits(parameters.m_numberOfCollectionPlaneHits),
    m_isShower(parameters.m_isShower),
    m_numberOfDownstreamParticles(parameters.m_numberOfDownstreamParticles),
    m_hasMcInfo(parameters.m_hasMcInfo),
    m_mcType(parameters.m_mcType),
    m_mcTypeTree(parameters.m_mcTypeTree),
    m_mcEnergy(parameters.m_mcEnergy),
    m_mcMomentum(parameters.m_mcMomentum),
    m_mcVertexPosition(parameters.m_mcVertexPosition),
    m_mcDirectionCosines(parameters.m_mcDirectionCosines),
    m_mcIsVertexFiducial(parameters.m_mcIsVertexFiducial),
    m_mcContainmentFraction(parameters.m_mcContainmentFraction),
    m_mcIsShower(parameters.m_mcIsShower),
    m_mcPdgCode(parameters.m_mcPdgCode),
    m_mcHitPurity(parameters.m_mcHitPurity),
    m_mcHitCompleteness(parameters.m_mcHitCompleteness),
    m_pMcMainMCParticle(parameters.m_pMcMainMCParticle)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArAnalysisParticle::Print() const
{
    if (this->m_hasMcInfo)
    {
        std::cout << "LArAnalysisParticle at " << this << ": \n"
              << "    - Type:                      " << this->TypeAsString(this->m_type) << "\n"
              << "    - MC type:                   " << this->TypeAsString(this->m_mcType) << "\n"
              << "    - Type tree:                 " << TEXT_GREEN_BOLD << this->TypeTreeAsString(this->m_typeTree) << TEXT_NORMAL << "\n"
              << "    - MC type tree:              " << TEXT_GREEN_BOLD << this->TypeTreeAsString(this->m_mcTypeTree) << TEXT_NORMAL << "\n"
              << "    - Energy:                    " << TEXT_MAGENTA_BOLD << 1000.f * this->m_analysisEnergy << "MeV\n" << TEXT_NORMAL
              << "    - MC energy:                 " << TEXT_MAGENTA_BOLD << 1000.f * this->m_mcEnergy << "MeV\n" << TEXT_NORMAL
              << "    - Energy from charge:        " << 1000.f * this->m_energyFromCharge << "MeV\n"
              << "    - Is vertex fiducial:        " << TEXT_RED_BOLD << std::boolalpha << this->m_isVertexFiducial << std::noboolalpha
                                                     << TEXT_NORMAL << "\n"
              << "    - MC is vertex fiducial:     " << TEXT_RED_BOLD << std::boolalpha << this->m_mcIsVertexFiducial << std::noboolalpha
                                                     << TEXT_NORMAL << "\n"
              << "    - Fiducial hit fraction:     " << TEXT_RED_BOLD << 100.f * this->m_fiducialHitFraction << TEXT_NORMAL << "%\n"
              << "    - MC containment fraction:   " << TEXT_RED_BOLD << 100.f * this->m_mcContainmentFraction << TEXT_NORMAL << "%\n"
              << "    - Vertex:                    " << "(" << this->m_vertexPosition.GetX() << ", " << this->m_vertexPosition.GetY() 
                                                     << ", " << this->m_vertexPosition.GetZ() << ")\n"
              << "    - MC vertex:                 " << "(" << this->m_mcVertexPosition.GetX() << ", " << this->m_mcVertexPosition.GetY()
                                                     << ", " << this->m_mcVertexPosition.GetZ() << ")\n"
              << "    - Direction cosines:         " << "(" << this->m_directionCosines.GetX() << ", " << this->m_directionCosines.GetY() 
                                                     << ", " << this->m_directionCosines.GetZ() << ")\n"
              << "    - MC direction cosines:      " << "(" << this->m_mcDirectionCosines.GetX() << ", " 
                                                     << this->m_mcDirectionCosines.GetY() << ", " << this->m_mcDirectionCosines.GetZ() 
                                                     << ")\n"
              << "    - Momentum:                  " << "(" << 1000.f * this->m_analysisMomentum.GetX() << ", " 
                                                     << 1000.f * this->m_analysisMomentum.GetY() << ", " 
                                                     << 1000.f * this->m_analysisMomentum.GetZ() << ") MeV/c\n"
              << "    - MC momentum:               " << "(" << 1000.f * this->m_mcMomentum.GetX() << ", " 
                                                     << 1000.f * this->m_mcMomentum.GetY() << ", " << 1000.f * this->m_mcMomentum.GetZ()
                                                     << ") MeV/c\n"
              << "    - Num 3D hits:               " << this->m_numberOf3dHits << "\n"
              << "    - Num collection-plane hits: " << this->m_numberOfCollectionPlaneHits << "\n"
              << "    - Is shower:                 " << std::boolalpha << this->m_isShower << std::noboolalpha << "\n"
              << "    - MC is shower:              " << std::boolalpha << this->m_mcIsShower << std::noboolalpha << "\n"
              << "    - Num downstream particles:  " << this->m_numberOfDownstreamParticles << "\n"
              << "    - Has MC info:               " << std::boolalpha << this->m_hasMcInfo << std::noboolalpha << "\n"
              << "    - MC PDG code:               " << this->m_mcPdgCode << "\n"
              << "    - MC hit purity:             " << 100.f * this->m_mcHitPurity << "%\n"
              << "    - MC hit completeness:       " << 100.f * this->m_mcHitCompleteness << "%\n"
              << "    - MC main MC particle at:    " << this->m_pMcMainMCParticle
              << std::endl;
    }
    
    else
    {
        std::cout << "LArAnalysisParticle at " << this << ": \n"
                  << "    - Type:                      " << this->TypeAsString(this->m_type) << "\n"
                  << "    - Type tree:                 " << TEXT_GREEN_BOLD << this->TypeTreeAsString(this->m_typeTree) << TEXT_NORMAL 
                                                         << "\n"
                  << "    - Energy:                    " << TEXT_MAGENTA_BOLD << 1000.f * this->m_analysisEnergy << "MeV\n" << TEXT_NORMAL
                  << "    - Energy from charge:        " << 1000.f * this->m_energyFromCharge << "MeV\n"
                  << "    - Is vertex fiducial:        " << TEXT_RED_BOLD << std::boolalpha << this->m_isVertexFiducial << std::noboolalpha
                                                         << TEXT_NORMAL << "\n"
                  << "    - Fiducial hit fraction:     " << TEXT_RED_BOLD << 100.f * this->m_fiducialHitFraction << TEXT_NORMAL << "%\n"
                  << "    - Vertex:                    " << "(" << this->m_vertexPosition.GetX() << ", " << this->m_vertexPosition.GetY()
                                                         << ", " << this->m_vertexPosition.GetZ() << ")\n"
                  << "    - Direction cosines:         " << "(" << this->m_directionCosines.GetX() << ", "
                                                         << this->m_directionCosines.GetY() << ", " << this->m_directionCosines.GetZ() 
                                                         << ")\n"
                  << "    - Momentum:                  " << "(" << 1000.f * this->m_analysisMomentum.GetX() << ", " 
                                                         << 1000.f * this->m_analysisMomentum.GetY() << ", " 
                                                         << 1000.f * this->m_analysisMomentum.GetZ() << ") MeV/c\n"
                  << "    - Num 3D hits:               " << this->m_numberOf3dHits << "\n"
                  << "    - Num collection-plane hits: " << this->m_numberOfCollectionPlaneHits << "\n"
                  << "    - Is shower:                 " << std::boolalpha << this->m_isShower << std::noboolalpha << "\n"
                  << "    - Num downstream particles:  " << this->m_numberOfDownstreamParticles << "\n"
                  << "    - Has MC info:               " << std::boolalpha << this->m_hasMcInfo << std::noboolalpha
                  << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::string LArAnalysisParticle::TypeAsString(const TYPE type)
{
    switch (type)
    {
        case TYPE::PION_MUON:  return "PION_MUON";
        case TYPE::PROTON:     return "PROTON";
        case TYPE::SHOWER:     return "SHOWER";
        case TYPE::TRACK:      return "TRACK";
        case TYPE::NEUTRINO:   return "NEUTRINO";
        case TYPE::COSMIC_RAY: return "COSMIC_RAY";
        default:               break;
    }
    
    return "UNKNOWN";
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::string LArAnalysisParticle::TypeTreeAsStringImpl(const TypeTree &typeTree, const bool printTrailingDelimiter)
{
    const std::string delimiter = " - ";
    std::string typeTreeString = TypeAsString(typeTree.Type());
    
    if (!typeTree.Daughters().empty())
    {
        typeTreeString += delimiter;
        typeTreeString += "[ ";
        
        for (auto iter = typeTree.Daughters().begin(); iter != typeTree.Daughters().end(); ++iter)
        {
            const bool isLast = (std::next(iter, 1) == typeTree.Daughters().end());
            typeTreeString += TypeTreeAsStringImpl(*iter, !isLast);
        }
            
        typeTreeString += " ]";
    }
    
    if (printTrailingDelimiter)
        typeTreeString += delimiter;
    
    return typeTreeString;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArAnalysisParticleParameters::LArAnalysisParticleParameters() noexcept :
    object_creation::ParticleFlowObject::Parameters(),
    m_type(LArAnalysisParticle::TYPE::UNKNOWN),
    m_typeTree(),
    m_analysisEnergy(0.f),
    m_energyFromCharge(0.f),
    m_isVertexFiducial(false),
    m_fiducialHitFraction(0.f),
    m_vertexPosition(CartesianVector(0.f, 0.f, 0.f)),
    m_directionCosines(CartesianVector(0.f, 0.f, 0.f)),
    m_analysisMomentum(CartesianVector(0.f, 0.f, 0.f)),
    m_numberOf3dHits(0U),
    m_numberOfCollectionPlaneHits(0U),
    m_isShower(false),
    m_numberOfDownstreamParticles(0U),
    m_hasMcInfo(false),
    m_mcType(LArAnalysisParticle::TYPE::UNKNOWN),
    m_mcTypeTree(),
    m_mcEnergy(0.f),
    m_mcMomentum(CartesianVector(0.f, 0.f, 0.f)),
    m_mcVertexPosition(CartesianVector(0.f, 0.f, 0.f)),
    m_mcDirectionCosines(CartesianVector(0.f, 0.f, 0.f)),
    m_mcIsVertexFiducial(false),
    m_mcContainmentFraction(0.f),
    m_mcIsShower(false),
    m_mcPdgCode(0),
    m_mcHitPurity(0.f),
    m_mcHitCompleteness(0.f),
    m_pMcMainMCParticle(nullptr)
{
}

} // namespace lar_physics_content
