/**
 *  @file   larphysicscontent/LArObjects/LArAnalysisParticle.cc
 *
 *  @brief  Implementation of the lar analysis particle class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArObjects/LArAnalysisParticle.h"

#include "larphysicscontent/LArHelpers/LArAnalysisParticleHelper.h"

using namespace pandora;

namespace lar_physics_content
{
LArAnalysisParticle::LArAnalysisParticle(const LArAnalysisParticleParameters &parameters) :
    ParticleFlowObject(parameters),
    m_type(parameters.m_type),
    m_typeTree(parameters.m_typeTree),
    m_kineticEnergy(parameters.m_kineticEnergy),
    m_isVertexFiducial(parameters.m_isVertexFiducial),
    m_fiducialHitFraction(parameters.m_fiducialHitFraction),
    m_vertexPosition(parameters.m_vertexPosition),
    m_directionCosines(parameters.m_directionCosines),
    m_analysisMomentum(parameters.m_analysisMomentum),
    m_numberOf3dHits(parameters.m_numberOf3dHits),
    m_numberOfCollectionPlaneHits(parameters.m_numberOfCollectionPlaneHits),
    m_isShower(parameters.m_isShower),
    m_numberOfDownstreamParticles(parameters.m_numberOfDownstreamParticles),
    m_kineticEnergyFromRangeFraction(parameters.m_kineticEnergyFromRangeFraction),
    m_kineticEnergyFromCorrectedTrackChargeFraction(parameters.m_kineticEnergyFromCorrectedTrackChargeFraction),
    m_kineticEnergyFromUncorrectedTrackChargeFraction(parameters.m_kineticEnergyFromUncorrectedTrackChargeFraction),
    m_kineticEnergyFromShowerChargeFraction(parameters.m_kineticEnergyFromShowerChargeFraction),
    m_hasMcInfo(parameters.m_hasMcInfo),
    m_mcType(parameters.m_mcType),
    m_mcTypeTree(parameters.m_mcTypeTree),
    m_mcEnergy(parameters.m_mcEnergy),
    m_mcKineticEnergy(parameters.m_mcKineticEnergy),
    m_mcMass(parameters.m_mcMass),
    m_mcMomentum(parameters.m_mcMomentum),
    m_mcVertexPosition(parameters.m_mcVertexPosition),
    m_mcDirectionCosines(parameters.m_mcDirectionCosines),
    m_mcIsVertexFiducial(parameters.m_mcIsVertexFiducial),
    m_mcContainmentFraction(parameters.m_mcContainmentFraction),
    m_mcIsShower(parameters.m_mcIsShower),
    m_mcPdgCode(parameters.m_mcPdgCode),
    m_mcHitPurity(parameters.m_mcHitPurity),
    m_mcHitCompleteness(parameters.m_mcHitCompleteness),
    m_mcCollectionPlaneHitPurity(parameters.m_mcCollectionPlaneHitPurity),
    m_mcCollectionPlaneHitCompleteness(parameters.m_mcCollectionPlaneHitCompleteness),
    m_pMcMainMCParticle(parameters.m_pMcMainMCParticle)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArAnalysisParticle::Print() const
{
    if (m_hasMcInfo)
    {
        std::cout << "LArAnalysisParticle at " << this << ": \n"
                  << "    - Type:                      " << LArAnalysisParticle::TypeAsString(m_type) << '\n'
                  << "    - MC type:                   " << LArAnalysisParticle::TypeAsString(m_mcType) << '\n'
                  << "    - Type tree:                 " << LArAnalysisParticleHelper::TypeTreeAsString(m_typeTree) << '\n'
                  << "    - MC type tree:              " << LArAnalysisParticleHelper::TypeTreeAsString(m_mcTypeTree) << '\n'
                  << "    - Kinetic energy:            " << 1000.f * m_kineticEnergy << "MeV\n"
                  << "    - MC kinetic energy:         " << 1000.f * m_mcKineticEnergy << "MeV\n"
                  << "    - MC energy:                 " << 1000.f * m_mcEnergy << "MeV\n"
                  << "    - MC mass:                   " << 1000.f * m_mcMass << "MeV/c^2\n"
                  << "    - KE from range:             " << 100.f * m_kineticEnergyFromRangeFraction << "%\n"
                  << "    - KE from corr track charge: " << 100.f * m_kineticEnergyFromCorrectedTrackChargeFraction << "%\n"
                  << "    - KE from track charge:      " << 100.f * m_kineticEnergyFromUncorrectedTrackChargeFraction << "%\n"
                  << "    - KE from shower charge:     " << 100.f * m_kineticEnergyFromShowerChargeFraction << "%\n"
                  << "    - Is vertex fiducial:        " << std::boolalpha << m_isVertexFiducial << std::noboolalpha << '\n'
                  << "    - MC is vertex fiducial:     " << std::boolalpha << m_mcIsVertexFiducial << std::noboolalpha << '\n'
                  << "    - Fiducial hit fraction:     " << 100.f * m_fiducialHitFraction << "%\n"
                  << "    - MC containment fraction:   " << 100.f * m_mcContainmentFraction << "%\n"
                  << "    - Vertex:                    "
                  << "(" << m_vertexPosition.GetX() << ", " << m_vertexPosition.GetY() << ", " << m_vertexPosition.GetZ() << ") cm\n"
                  << "    - MC vertex:                 "
                  << "(" << m_mcVertexPosition.GetX() << ", " << m_mcVertexPosition.GetY() << ", " << m_mcVertexPosition.GetZ() << ") cm\n"
                  << "    - Direction cosines:         "
                  << "(" << m_directionCosines.GetX() << ", " << m_directionCosines.GetY() << ", " << m_directionCosines.GetZ() << ")\n"
                  << "    - MC direction cosines:      "
                  << "(" << m_mcDirectionCosines.GetX() << ", " << m_mcDirectionCosines.GetY() << ", " << m_mcDirectionCosines.GetZ() << ")\n"
                  << "    - Momentum:                  "
                  << "(" << 1000.f * m_analysisMomentum.GetX() << ", " << 1000.f * m_analysisMomentum.GetY() << ", "
                  << 1000.f * m_analysisMomentum.GetZ() << ") MeV/c\n"
                  << "    - MC momentum:               "
                  << "(" << 1000.f * m_mcMomentum.GetX() << ", " << 1000.f * m_mcMomentum.GetY() << ", " << 1000.f * m_mcMomentum.GetZ() << ") MeV/c\n"
                  << "    - Num 3D hits:               " << m_numberOf3dHits << '\n'
                  << "    - Num collection-plane hits: " << m_numberOfCollectionPlaneHits << '\n'
                  << "    - Is shower:                 " << std::boolalpha << m_isShower << std::noboolalpha << '\n'
                  << "    - MC is shower:              " << std::boolalpha << m_mcIsShower << std::noboolalpha << '\n'
                  << "    - Num downstream particles:  " << m_numberOfDownstreamParticles << '\n'
                  << "    - Has MC info:               " << std::boolalpha << m_hasMcInfo << std::noboolalpha << '\n'
                  << "    - MC PDG code:               " << m_mcPdgCode << '\n'
                  << "    - MC hit purity:             " << 100.f * m_mcHitPurity << "%\n"
                  << "    - MC hit completeness:       " << 100.f * m_mcHitCompleteness << "%\n"
                  << "    - MC hit purity (W):         " << 100.f * m_mcCollectionPlaneHitPurity << "%\n"
                  << "    - MC hit completeness (W):   " << 100.f * m_mcCollectionPlaneHitCompleteness << "%\n"
                  << "    - MC main MC particle at:    " << m_pMcMainMCParticle << std::endl;
    }

    else
    {
        std::cout << "LArAnalysisParticle at " << this << ": \n"
                  << "    - Type:                      " << LArAnalysisParticle::TypeAsString(m_type) << '\n'
                  << "    - Type tree:                 " << LArAnalysisParticleHelper::TypeTreeAsString(m_typeTree) << "\n"
                  << "    - Kinetic energy:            " << 1000.f * m_kineticEnergy << "MeV\n"
                  << "    - KE from range:             " << 100.f * m_kineticEnergyFromRangeFraction << "%\n"
                  << "    - KE from corr track charge: " << 100.f * m_kineticEnergyFromCorrectedTrackChargeFraction << "%\n"
                  << "    - KE from track charge:      " << 100.f * m_kineticEnergyFromUncorrectedTrackChargeFraction << "%\n"
                  << "    - KE from shower charge:     " << 100.f * m_kineticEnergyFromShowerChargeFraction << "%\n"
                  << "    - Is vertex fiducial:        " << std::boolalpha << m_isVertexFiducial << std::noboolalpha << '\n'
                  << "    - Fiducial hit fraction:     " << 100.f * m_fiducialHitFraction << "%\n"
                  << "    - Vertex:                    "
                  << "(" << m_vertexPosition.GetX() << ", " << m_vertexPosition.GetY() << ", " << m_vertexPosition.GetZ() << ") cm\n"
                  << "    - Direction cosines:         "
                  << "(" << m_directionCosines.GetX() << ", " << m_directionCosines.GetY() << ", " << m_directionCosines.GetZ() << ")\n"
                  << "    - Momentum:                  "
                  << "(" << 1000.f * m_analysisMomentum.GetX() << ", " << 1000.f * m_analysisMomentum.GetY() << ", "
                  << 1000.f * m_analysisMomentum.GetZ() << ") MeV/c\n"
                  << "    - Num 3D hits:               " << m_numberOf3dHits << '\n'
                  << "    - Num collection-plane hits: " << m_numberOfCollectionPlaneHits << '\n'
                  << "    - Is shower:                 " << std::boolalpha << m_isShower << std::noboolalpha << '\n'
                  << "    - Num downstream particles:  " << m_numberOfDownstreamParticles << '\n'
                  << "    - Has MC info:               " << std::boolalpha << m_hasMcInfo << std::noboolalpha << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::string LArAnalysisParticle::TypeAsString(const LArAnalysisParticle::TYPE type)
{
    switch (type)
    {
        case TYPE::PION_MUON:
            return "PION_MUON";
        case TYPE::PROTON:
            return "PROTON";
        case TYPE::SHOWER:
            return "SHOWER";
        case TYPE::TRACK:
            return "TRACK";
        case TYPE::NEUTRINO:
            return "NEUTRINO";
        case TYPE::COSMIC_RAY:
            return "COSMIC_RAY";
        default:
            break;
    }

    return "UNKNOWN";
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArAnalysisParticleParameters::LArAnalysisParticleParameters() noexcept :
    object_creation::ParticleFlowObject::Parameters(),
    m_type(LArAnalysisParticle::TYPE::UNKNOWN),
    m_typeTree(),
    m_kineticEnergy(0.f),
    m_isVertexFiducial(false),
    m_fiducialHitFraction(0.f),
    m_vertexPosition(CartesianVector(0.f, 0.f, 0.f)),
    m_directionCosines(CartesianVector(0.f, 0.f, 0.f)),
    m_analysisMomentum(CartesianVector(0.f, 0.f, 0.f)),
    m_numberOf3dHits(0U),
    m_numberOfCollectionPlaneHits(0U),
    m_isShower(false),
    m_numberOfDownstreamParticles(0U),
    m_kineticEnergyFromRangeFraction(0.f),
    m_kineticEnergyFromCorrectedTrackChargeFraction(0.f),
    m_kineticEnergyFromUncorrectedTrackChargeFraction(0.f),
    m_kineticEnergyFromShowerChargeFraction(0.f),
    m_hasMcInfo(false),
    m_mcType(LArAnalysisParticle::TYPE::UNKNOWN),
    m_mcTypeTree(),
    m_mcEnergy(0.f),
    m_mcKineticEnergy(0.f),
    m_mcMass(0.f),
    m_mcMomentum(CartesianVector(0.f, 0.f, 0.f)),
    m_mcVertexPosition(CartesianVector(0.f, 0.f, 0.f)),
    m_mcDirectionCosines(CartesianVector(0.f, 0.f, 0.f)),
    m_mcIsVertexFiducial(false),
    m_mcContainmentFraction(0.f),
    m_mcIsShower(false),
    m_mcPdgCode(0),
    m_mcHitPurity(0.f),
    m_mcHitCompleteness(0.f),
    m_mcCollectionPlaneHitPurity(0.f),
    m_mcCollectionPlaneHitCompleteness(0.f),
    m_pMcMainMCParticle(nullptr)
{
}

} // namespace lar_physics_content
