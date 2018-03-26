/**
 *  @file   larphysicscontent/WriteAnalysisParticlesAlgorithm.h
 *
 *  @brief  Header file for the write AnalysisParticles algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_WRITE_ANALYSIS_PARTICLES_ALGORITHM_H
#define LAR_WRITE_ANALYSIS_PARTICLES_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larphysicscontent/LArAnalysisParticle.h"
#include "larphysicscontent/LArAnalysisParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

/**
 *  @brief  Macro for performing operations on the tree member scalars: (variable name, member variable, type, default value, units)
 *
 *  @param  d the macro
 */
#define TREE_SCALAR_MEMBERS(d)                                                                                                       \
    d("nu_WasReconstructed",                            m_nu_WasReconstructed,                            Bool_t,    false, "")      \
    d("nu_IsVertexFiducial",                            m_nu_IsVertexFiducial,                            Bool_t,    false, "")      \
    d("nu_IsContained",                                 m_nu_IsContained,                                 Bool_t,    false, "")      \
    d("nu_FiducialHitFraction",                         m_nu_FiducialHitFraction,                         Float_t,   0.f,   "")      \
    d("nu_HasMcInfo",                                   m_nu_HasMcInfo,                                   Bool_t,    false, "")      \
    d("nu_VisibleEnergy",                               m_nu_VisibleEnergy,                               Float_t,   0.f,   "GeV")   \
    d("nu_VisibleLongitudinalEnergy",                   m_nu_VisibleLongitudinalEnergy,                   Float_t,   0.f,   "GeV")   \
    d("nu_VisibleTransverseEnergy",                     m_nu_VisibleTransverseEnergy,                     Float_t,   0.f,   "GeV")   \
    d("nu_VisibleEnergyFracFromRange",                  m_nu_VisibleEnergyFracFromRange,                  Float_t,   0.f,   "")      \
    d("nu_VisibleEnergyFracFromCorrectedTrackCharge",   m_nu_VisibleEnergyFracFromCorrectedTrackCharge,   Float_t,   0.f,   "")      \
    d("nu_VisibleEnergyFracFromUncorrectedTrackCharge", m_nu_VisibleEnergyFracFromUncorrectedTrackCharge, Float_t,   0.f,   "")      \
    d("nu_VisibleEnergyFracFromShowerCharge",           m_nu_VisibleEnergyFracFromShowerCharge,           Float_t,   0.f,   "")      \
    d("nu_VertexX",                                     m_nu_VertexX,                                     Float_t,   0.f,   "cm")    \
    d("nu_VertexY",                                     m_nu_VertexY,                                     Float_t,   0.f,   "cm")    \
    d("nu_VertexZ",                                     m_nu_VertexZ,                                     Float_t,   0.f,   "cm")    \
    d("nu_DirectionCosineX",                            m_nu_DirectionCosineX,                            Float_t,   0.f,   "")      \
    d("nu_DirectionCosineY",                            m_nu_DirectionCosineY,                            Float_t,   0.f,   "")      \
    d("nu_DirectionCosineZ",                            m_nu_DirectionCosineZ,                            Float_t,   0.f,   "")      \
    d("nu_Momentum",                                    m_nu_Momentum,                                    Float_t,   0.f,   "GeV/c") \
    d("nu_MomentumX",                                   m_nu_MomentumX,                                   Float_t,   0.f,   "GeV/c") \
    d("nu_MomentumY",                                   m_nu_MomentumY,                                   Float_t,   0.f,   "GeV/c") \
    d("nu_MomentumZ",                                   m_nu_MomentumZ,                                   Float_t,   0.f,   "GeV/c") \
    d("nu_TypeTree",                                    m_nu_TypeTree,                                    TString,   "",    "")      \
    d("nu_NumberOf3dHits",                              m_nu_NumberOf3dHits,                              UInt_t,    0U,    "")      \
    d("nu_NumberOfCollectionPlaneHits",                 m_nu_NumberOfCollectionPlaneHits,                 UInt_t,    0U,    "")      \
    d("nu_NumberOfRecoParticles",                       m_nu_NumberOfRecoParticles,                       UInt_t,    0U,    "")      \
    d("nu_NumberOfRecoTracks",                          m_nu_NumberOfRecoTracks,                          UInt_t,    0U,    "")      \
    d("nu_NumberOfRecoShowers",                         m_nu_NumberOfRecoShowers,                         UInt_t,    0UL,   "")      \
    d("nu_mc_McParticleUid",                            m_nu_mc_McParticleUid,                            ULong64_t, 0ULL,  "")      \
    d("nu_mc_Energy",                                   m_nu_mc_Energy,                                   Float_t,   0.f,   "GeV")   \
    d("nu_mc_LongitudinalEnergy",                       m_nu_mc_LongitudinalEnergy,                       Float_t,   0.f,   "GeV")   \
    d("nu_mc_TransverseEnergy",                         m_nu_mc_TransverseEnergy,                         Float_t,   0.f,   "GeV")   \
    d("nu_mc_VisibleEnergy",                            m_nu_mc_VisibleEnergy,                            Float_t,   0.f,   "GeV")   \
    d("nu_mc_VisibleLongitudinalEnergy",                m_nu_mc_VisibleLongitudinalEnergy,                Float_t,   0.f,   "GeV")   \
    d("nu_mc_VisibleTransverseEnergy",                  m_nu_mc_VisibleTransverseEnergy,                  Float_t,   0.f,   "GeV")   \
    d("nu_mc_VertexX",                                  m_nu_mc_VertexX,                                  Float_t,   0.f,   "cm")    \
    d("nu_mc_VertexY",                                  m_nu_mc_VertexY,                                  Float_t,   0.f,   "cm")    \
    d("nu_mc_VertexZ",                                  m_nu_mc_VertexZ,                                  Float_t,   0.f,   "cm")    \
    d("nu_mc_DirectionCosineX",                         m_nu_mc_DirectionCosineX,                         Float_t,   0.f,   "")      \
    d("nu_mc_DirectionCosineY",                         m_nu_mc_DirectionCosineY,                         Float_t,   0.f,   "")      \
    d("nu_mc_DirectionCosineZ",                         m_nu_mc_DirectionCosineZ,                         Float_t,   0.f,   "")      \
    d("nu_mc_Momentum",                                 m_nu_mc_Momentum,                                 Float_t,   0.f,   "GeV/c") \
    d("nu_mc_MomentumX",                                m_nu_mc_MomentumX,                                Float_t,   0.f,   "GeV/c") \
    d("nu_mc_MomentumY",                                m_nu_mc_MomentumY,                                Float_t,   0.f,   "GeV/c") \
    d("nu_mc_MomentumZ",                                m_nu_mc_MomentumZ,                                Float_t,   0.f,   "GeV/c") \
    d("nu_mc_IsVertexFiducial",                         m_nu_mc_IsVertexFiducial,                         Bool_t,    false, "")      \
    d("nu_mc_IsContained",                              m_nu_mc_IsContained,                              Bool_t,    false, "")      \
    d("nu_mc_ContainmentFraction",                      m_nu_mc_ContainmentFraction,                      Float_t,   0.f,   "")      \
    d("nu_mc_TypeTree",                                 m_nu_mc_TypeTree,                                 TString,   "",    "")      \
    d("nu_mc_InteractionType",                          m_nu_mc_InteractionType,                          TString,   "",    "")      \
    d("nu_mc_IsChargedCurrent",                         m_nu_mc_IsChargedCurrent,                         Bool_t,    false, "")      \
    d("nu_mc_VisibleEnergyFraction",                    m_nu_mc_VisibleEnergyFraction,                    Float_t,   0.f,   "")      \
    d("nu_mc_PdgCode",                                  m_nu_mc_PdgCode,                                  Int_t,     0,     "")      \
    d("primary_Number",                                 m_primary_Number,                                 UInt_t,    0U,    "")      \
    d("cr_Number",                                      m_cr_Number,                                      UInt_t,    0U,    "")

/**
 *  @brief  Macro for performing operations on the primary non-MC tree member vectors: (variable name, member variable, type, default value, units, size)
 *
 *  @param  d the macro
 */
#define TREE_VECTOR_MEMBERS_PRIMARY(d)                                                                                                                                    \
    d("primary_WasReconstructed",                            m_primary_WasReconstructed,                            RBoolVector,      false, "",        m_primary_Number) \
    d("primary_IsVertexFiducial",                            m_primary_IsVertexFiducial,                            RBoolVector,      false, "",        m_primary_Number) \
    d("primary_IsContained",                                 m_primary_IsContained,                                 RBoolVector,      false, "",        m_primary_Number) \
    d("primary_FiducialHitFraction",                         m_primary_FiducialHitFraction,                         RFloatVector,     0.f,   "",        m_primary_Number) \
    d("primary_HasMcInfo",                                   m_primary_HasMcInfo,                                   RBoolVector,      false, "",        m_primary_Number) \
    d("primary_KineticEnergy",                               m_primary_KineticEnergy,                               RFloatVector,     0.f,   "GeV",     m_primary_Number) \
    d("primary_KineticEnergyFracFromRange",                  m_primary_KineticEnergyFracFromRange,                  RFloatVector,     0.f,   "",        m_primary_Number) \
    d("primary_KineticEnergyFracFromCorrectedTrackCharge",   m_primary_KineticEnergyFracFromCorrectedTrackCharge,   RFloatVector,     0.f,   "",        m_primary_Number) \
    d("primary_KineticEnergyFracFromUncorrectedTrackCharge", m_primary_KineticEnergyFracFromUncorrectedTrackCharge, RFloatVector,     0.f,   "",        m_primary_Number) \
    d("primary_KineticEnergyFracFromShowerCharge",           m_primary_KineticEnergyFracFromShowerCharge,           RFloatVector,     0.f,   "",        m_primary_Number) \
    d("primary_VertexX",                                     m_primary_VertexX,                                     RFloatVector,     0.f,   "cm",      m_primary_Number) \
    d("primary_VertexY",                                     m_primary_VertexY,                                     RFloatVector,     0.f,   "cm",      m_primary_Number) \
    d("primary_VertexZ",                                     m_primary_VertexZ,                                     RFloatVector,     0.f,   "cm",      m_primary_Number) \
    d("primary_DirectionCosineX",                            m_primary_DirectionCosineX,                            RFloatVector,     0.f,   "",        m_primary_Number) \
    d("primary_DirectionCosineY",                            m_primary_DirectionCosineY,                            RFloatVector,     0.f,   "",        m_primary_Number) \
    d("primary_DirectionCosineZ",                            m_primary_DirectionCosineZ,                            RFloatVector,     0.f,   "",        m_primary_Number) \
    d("primary_Momentum",                                    m_primary_Momentum,                                    RFloatVector,     0.f,   "GeV/c",   m_primary_Number) \
    d("primary_MomentumX",                                   m_primary_MomentumX,                                   RFloatVector,     0.f,   "GeV/c",   m_primary_Number) \
    d("primary_MomentumY",                                   m_primary_MomentumY,                                   RFloatVector,     0.f,   "GeV/c",   m_primary_Number) \
    d("primary_MomentumZ",                                   m_primary_MomentumZ,                                   RFloatVector,     0.f,   "GeV/c",   m_primary_Number) \
    d("primary_TypeTree",                                    m_primary_TypeTree,                                    RTextVector,      "",    "",        m_primary_Number) \
    d("primary_IsShower",                                    m_primary_IsShower,                                    RBoolVector,      false, "",        m_primary_Number) \
    d("primary_IsTrack",                                     m_primary_IsTrack,                                     RBoolVector,      false, "",        m_primary_Number) \
    d("primary_IsProton",                                    m_primary_IsProton,                                    RBoolVector,      false, "",        m_primary_Number) \
    d("primary_IsPionOrMuon",                                m_primary_IsPionOrMuon,                                RBoolVector,      false, "",        m_primary_Number) \
    d("primary_NumberOf3dHits",                              m_primary_NumberOf3dHits,                              RUnsignedVector,  0U,    "",        m_primary_Number) \
    d("primary_NumberOfCollectionPlaneHits",                 m_primary_NumberOfCollectionPlaneHits,                 RUnsignedVector,  0U,    "",        m_primary_Number) \
    d("primary_NumberOfDownstreamParticles",                 m_primary_NumberOfDownstreamParticles,                 RUnsignedVector,  0U,    "",        m_primary_Number)

/**
 *  @brief  Macro for performing operations on the primary MC tree member vectors: (variable name, member variable, type, default value, units, size)
 *
 *  @param  d the macro
 */
#define TREE_VECTOR_MEMBERS_PRIMARY_MC(d)                                                                                                                                 \
    d("primary_mc_McParticleUid",                            m_primary_mc_McParticleUid,                            RUInt64Vector,    0LL,   "",        m_primary_Number) \
    d("primary_mc_IsParticleSplitByReco",                    m_primary_mc_IsParticleSplitByReco,                    RBoolVector,      false, "",        m_primary_Number) \
    d("primary_mc_Energy",                                   m_primary_mc_Energy,                                   RFloatVector,     0.f,   "GeV",     m_primary_Number) \
    d("primary_mc_KineticEnergy",                            m_primary_mc_KineticEnergy,                            RFloatVector,     0.f,   "GeV",     m_primary_Number) \
    d("primary_mc_Mass",                                     m_primary_mc_Mass,                                     RFloatVector,     0.f,   "GeV/c^2", m_primary_Number) \
    d("primary_mc_VertexX",                                  m_primary_mc_VertexX,                                  RFloatVector,     0.f,   "cm",      m_primary_Number) \
    d("primary_mc_VertexY",                                  m_primary_mc_VertexY,                                  RFloatVector,     0.f,   "cm",      m_primary_Number) \
    d("primary_mc_VertexZ",                                  m_primary_mc_VertexZ,                                  RFloatVector,     0.f,   "cm",      m_primary_Number) \
    d("primary_mc_DirectionCosineX",                         m_primary_mc_DirectionCosineX,                         RFloatVector,     0.f,   "",        m_primary_Number) \
    d("primary_mc_DirectionCosineY",                         m_primary_mc_DirectionCosineY,                         RFloatVector,     0.f,   "",        m_primary_Number) \
    d("primary_mc_DirectionCosineZ",                         m_primary_mc_DirectionCosineZ,                         RFloatVector,     0.f,   "",        m_primary_Number) \
    d("primary_mc_Momentum",                                 m_primary_mc_Momentum,                                 RFloatVector,     0.f,   "GeV/c",   m_primary_Number) \
    d("primary_mc_MomentumX",                                m_primary_mc_MomentumX,                                RFloatVector,     0.f,   "GeV/c",   m_primary_Number) \
    d("primary_mc_MomentumY",                                m_primary_mc_MomentumY,                                RFloatVector,     0.f,   "GeV/c",   m_primary_Number) \
    d("primary_mc_MomentumZ",                                m_primary_mc_MomentumZ,                                RFloatVector,     0.f,   "GeV/c",   m_primary_Number) \
    d("primary_mc_IsVertexFiducial",                         m_primary_mc_IsVertexFiducial,                         RBoolVector,      false, "",        m_primary_Number) \
    d("primary_mc_IsContained",                              m_primary_mc_IsContained,                              RBoolVector,      false, "",        m_primary_Number) \
    d("primary_mc_ContainmentFraction",                      m_primary_mc_ContainmentFraction,                      RFloatVector,     0.f,   "",        m_primary_Number) \
    d("primary_mc_TypeTree",                                 m_primary_mc_TypeTree,                                 RTextVector,      "",    "",        m_primary_Number) \
    d("primary_mc_IsShower",                                 m_primary_mc_IsShower,                                 RBoolVector,      false, "",        m_primary_Number) \
    d("primary_mc_IsTrack",                                  m_primary_mc_IsTrack,                                  RBoolVector,      false, "",        m_primary_Number) \
    d("primary_mc_IsProton",                                 m_primary_mc_IsProton,                                 RBoolVector,      false, "",        m_primary_Number) \
    d("primary_mc_IsPionOrMuon",                             m_primary_mc_IsPionOrMuon,                             RBoolVector,      false, "",        m_primary_Number) \
    d("primary_mc_IsCosmicRay",                              m_primary_mc_IsCosmicRay,                              RBoolVector,      false, "",        m_primary_Number) \
    d("primary_mc_PdgCode",                                  m_primary_mc_PdgCode,                                  RIntVector,       0,     "",        m_primary_Number)

/**
 *  @brief  Macro for performing operations on the CR non-MC tree member vectors: (variable name, member variable, type, default value, units, size)
 *
 *  @param  d the macro
 */
#define TREE_VECTOR_MEMBERS_CR(d)                                                                                                                                    \
    d("cr_WasReconstructed",                                 m_cr_WasReconstructed,                                 RBoolVector,      false, "",        m_cr_Number) \
    d("cr_IsVertexFiducial",                                 m_cr_IsVertexFiducial,                                 RBoolVector,      false, "",        m_cr_Number) \
    d("cr_IsContained",                                      m_cr_IsContained,                                      RBoolVector,      false, "",        m_cr_Number) \
    d("cr_FiducialHitFraction",                              m_cr_FiducialHitFraction,                              RFloatVector,     0.f,   "",        m_cr_Number) \
    d("cr_HasMcInfo",                                        m_cr_HasMcInfo,                                        RBoolVector,      false, "",        m_cr_Number) \
    d("cr_KineticEnergy",                                    m_cr_KineticEnergy,                                    RFloatVector,     0.f,   "GeV",     m_cr_Number) \
    d("cr_KineticEnergyFracFromRange",                       m_cr_KineticEnergyFracFromRange,                       RFloatVector,     0.f,   "",        m_cr_Number) \
    d("cr_KineticEnergyFracFromCorrectedTrackCharge",        m_cr_KineticEnergyFracFromCorrectedTrackCharge,        RFloatVector,     0.f,   "",        m_cr_Number) \
    d("cr_KineticEnergyFracFromUncorrectedTrackCharge",      m_cr_KineticEnergyFracFromUncorrectedTrackCharge,      RFloatVector,     0.f,   "",        m_cr_Number) \
    d("cr_KineticEnergyFracFromShowerCharge",                m_cr_KineticEnergyFracFromShowerCharge,                RFloatVector,     0.f,   "",        m_cr_Number) \
    d("cr_VertexX",                                          m_cr_VertexX,                                          RFloatVector,     0.f,   "cm",      m_cr_Number) \
    d("cr_VertexY",                                          m_cr_VertexY,                                          RFloatVector,     0.f,   "cm",      m_cr_Number) \
    d("cr_VertexZ",                                          m_cr_VertexZ,                                          RFloatVector,     0.f,   "cm",      m_cr_Number) \
    d("cr_DirectionCosineX",                                 m_cr_DirectionCosineX,                                 RFloatVector,     0.f,   "",        m_cr_Number) \
    d("cr_DirectionCosineY",                                 m_cr_DirectionCosineY,                                 RFloatVector,     0.f,   "",        m_cr_Number) \
    d("cr_DirectionCosineZ",                                 m_cr_DirectionCosineZ,                                 RFloatVector,     0.f,   "",        m_cr_Number) \
    d("cr_Momentum",                                         m_cr_Momentum,                                         RFloatVector,     0.f,   "GeV/c",   m_cr_Number) \
    d("cr_MomentumX",                                        m_cr_MomentumX,                                        RFloatVector,     0.f,   "GeV/c",   m_cr_Number) \
    d("cr_MomentumY",                                        m_cr_MomentumY,                                        RFloatVector,     0.f,   "GeV/c",   m_cr_Number) \
    d("cr_MomentumZ",                                        m_cr_MomentumZ,                                        RFloatVector,     0.f,   "GeV/c",   m_cr_Number) \
    d("cr_TypeTree",                                         m_cr_TypeTree,                                         RTextVector,      "",    "",        m_cr_Number) \
    d("cr_NumberOf3dHits",                                   m_cr_NumberOf3dHits,                                   RUnsignedVector,  0U,    "",        m_cr_Number) \
    d("cr_NumberOfCollectionPlaneHits",                      m_cr_NumberOfCollectionPlaneHits,                      RUnsignedVector,  0U,    "",        m_cr_Number) \
    d("cr_NumberOfDownstreamParticles",                      m_cr_NumberOfDownstreamParticles,                      RUnsignedVector,  0U,    "",        m_cr_Number)

/**
 *  @brief  Macro for performing operations on the CR MC tree member vectors: (variable name, member variable, type, default value, units, size)
 *
 *  @param  d the macro
 */
#define TREE_VECTOR_MEMBERS_CR_MC(d)                                                                                                                                 \
    d("cr_mc_McParticleUid",                                 m_cr_mc_McParticleUid,                                 RUInt64Vector,    0ULL,  "",        m_cr_Number) \
    d("cr_mc_IsParticleSplitByReco",                         m_cr_mc_IsParticleSplitByReco,                         RBoolVector,      false, "",        m_cr_Number) \
    d("cr_mc_Energy",                                        m_cr_mc_Energy,                                        RFloatVector,     0.f,   "GeV",     m_cr_Number) \
    d("cr_mc_KineticEnergy",                                 m_cr_mc_KineticEnergy,                                 RFloatVector,     0.f,   "GeV",     m_cr_Number) \
    d("cr_mc_Mass",                                          m_cr_mc_Mass,                                          RFloatVector,     0.f,   "GeV/c^2", m_cr_Number) \
    d("cr_mc_VertexX",                                       m_cr_mc_VertexX,                                       RFloatVector,     0.f,   "cm",      m_cr_Number) \
    d("cr_mc_VertexY",                                       m_cr_mc_VertexY,                                       RFloatVector,     0.f,   "cm",      m_cr_Number) \
    d("cr_mc_VertexZ",                                       m_cr_mc_VertexZ,                                       RFloatVector,     0.f,   "cm",      m_cr_Number) \
    d("cr_mc_DirectionCosineX",                              m_cr_mc_DirectionCosineX,                              RFloatVector,     0.f,   "",        m_cr_Number) \
    d("cr_mc_DirectionCosineY",                              m_cr_mc_DirectionCosineY,                              RFloatVector,     0.f,   "",        m_cr_Number) \
    d("cr_mc_DirectionCosineZ",                              m_cr_mc_DirectionCosineZ,                              RFloatVector,     0.f,   "",        m_cr_Number) \
    d("cr_mc_Momentum",                                      m_cr_mc_Momentum,                                      RFloatVector,     0.f,   "GeV/c",   m_cr_Number) \
    d("cr_mc_MomentumX",                                     m_cr_mc_MomentumX,                                     RFloatVector,     0.f,   "GeV/c",   m_cr_Number) \
    d("cr_mc_MomentumY",                                     m_cr_mc_MomentumY,                                     RFloatVector,     0.f,   "GeV/c",   m_cr_Number) \
    d("cr_mc_MomentumZ",                                     m_cr_mc_MomentumZ,                                     RFloatVector,     0.f,   "GeV/c",   m_cr_Number) \
    d("cr_mc_IsVertexFiducial",                              m_cr_mc_IsVertexFiducial,                              RBoolVector,      false, "",        m_cr_Number) \
    d("cr_mc_IsContained",                                   m_cr_mc_IsContained,                                   RBoolVector,      false, "",        m_cr_Number) \
    d("cr_mc_ContainmentFraction",                           m_cr_mc_ContainmentFraction,                           RFloatVector,     0.f,   "",        m_cr_Number) \
    d("cr_mc_TypeTree",                                      m_cr_mc_TypeTree,                                      RTextVector,      "",    "",        m_cr_Number) \
    d("cr_mc_IsShower",                                      m_cr_mc_IsShower,                                      RBoolVector,      false, "",        m_cr_Number) \
    d("cr_mc_IsTrack",                                       m_cr_mc_IsTrack,                                       RBoolVector,      false, "",        m_cr_Number) \
    d("cr_mc_IsProton",                                      m_cr_mc_IsProton,                                      RBoolVector,      false, "",        m_cr_Number) \
    d("cr_mc_IsPionOrMuon",                                  m_cr_mc_IsPionOrMuon,                                  RBoolVector,      false, "",        m_cr_Number) \
    d("cr_mc_IsCosmicRay",                                   m_cr_mc_IsCosmicRay,                                   RBoolVector,      false, "",        m_cr_Number) \
    d("cr_mc_PdgCode",                                       m_cr_mc_PdgCode,                                       RIntVector,       0,     "",        m_cr_Number)

/**
 *  @brief  Declare a set of scalar members
 */
#define DECLARE_SCALAR_MEMBER(name, memberVariable, memberType, defaultValue, units) \
    memberType memberVariable;

/**
 *  @brief  Declare a set of vector members
 */
#define DECLARE_VECTOR_MEMBER(name, memberVariable, memberType, defaultValue, units, sizeMember) \
    memberType memberVariable;

/**
 *  @brief  Initialize a set of scalar members
 */
#define INITIALIZE_SCALAR_MEMBER(name, memberVariable, memberType, defaultValue, units) \
    memberVariable(defaultValue),

/**
 *  @brief  Initialize a set of vector members
 */
#define INITIALIZE_VECTOR_MEMBER(name, memberVariable, memberType, defaultValue, units, sizeMember) \
    memberVariable(),

/**
 *  @brief  Set a set of scalar member TTree branches
 */
#define SET_SCALAR_MEMBER_BRANCH(name, memberVariable, memberType, defaultValue, units) \
    m_pOutputTree->Branch(name, &m_treeParameters.memberVariable);

/**
 *  @brief  Set a set of vector member TTree branches
 */
#define SET_VECTOR_MEMBER_BRANCH(name, memberVariable, memberType, defaultValue, units, sizeMember) \
    m_pOutputTree->Branch(name, &m_treeParameters.memberVariable);

/**
 *  @brief  Check the sizes of a set of vector members
 */
#define CHECK_VECTOR_MEMBER_SIZE(name, memberVariable, memberType, defaultValue, units, sizeMember) \
    this->CheckSize(name, m_treeParameters.memberVariable, m_treeParameters.sizeMember);

/**
 *  @brief  Push the default values to a set of vector members
 */
#define VECTOR_MEMBER_PUSH_DEFAULT(name, memberVariable, memberType, defaultValue, units, sizeMember) \
    m_treeParameters.memberVariable.push_back(defaultValue);

/**
 *  @brief  Macro for pushing a new value to a tree member vector
 */
#define PUSH_TREE_RECORD(treeMember, value) \
    m_treeParameters.treeMember.push_back(value)

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{

/**
 *  @brief  TreeParameters class
 */
class TreeParameters
{
public:
    using RFloatVector    = std::vector<Float_t>;   ///< Alias for a vector of ROOT floats
    using RIntVector      = std::vector<Int_t>;     ///< Alias for a vector of ROOT ints
    using RBoolVector     = std::vector<Bool_t>;    ///< Alias for a vector of ROOT bools
    using RUnsignedVector = std::vector<UInt_t>;    ///< Alias for a vector of ROOT unsigned ints
    using RUInt64Vector   = std::vector<ULong64_t>; ///< Alias for a vector of ROOT unsigned 64-bit ints
    using RTextVector     = std::vector<TString>;   ///< Alias for a vector of ROOT strings

    /**
     *  @brief  Default constructor
     */
    TreeParameters() noexcept;

    TREE_SCALAR_MEMBERS(DECLARE_SCALAR_MEMBER)
    TREE_VECTOR_MEMBERS_PRIMARY(DECLARE_VECTOR_MEMBER)
    TREE_VECTOR_MEMBERS_PRIMARY_MC(DECLARE_VECTOR_MEMBER)
    TREE_VECTOR_MEMBERS_CR(DECLARE_VECTOR_MEMBER)
    TREE_VECTOR_MEMBERS_CR_MC(DECLARE_VECTOR_MEMBER)

    bool  m_dummy; ///< A dummy variable for the preprocessor trick.
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  WriteAnalysisParticlesAlgorithm class
 *
 */
class WriteAnalysisParticlesAlgorithm : public Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    WriteAnalysisParticlesAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~WriteAnalysisParticlesAlgorithm();

protected:
    StatusCode ReadSettings(const TiXmlHandle xmlHandle);
    StatusCode Run();

private:
    using AnalysisParticleList = std::list<const LArAnalysisParticle *>; ///< Alias for a list of AnalysisParticles
    using MCPrimaryMap = std::unordered_multimap<const MCParticle *, const LArAnalysisParticle *>; ///< Alias for a map from MC primaries to AnalysisParticles

    /**
     *  @brief  Get the map from the MC primaries to their analysis particles
     *
     *  @param  pfoList the list of PFOs
     *
     *  @return the map
     */
    MCPrimaryMap GetMainMcParticleMap(const PfoList &pfoList) const;

    /**
     *  @brief  Write the parameters for any analysis particle
     *
     *  @param  pAnalysisParticle address of the analysis particle
     *  @param  mainMcParticleMap the map from MC primaries to their analysis particles
     *  @param  pMCParticleList address of the MC particle list
     *
     *  @return success
     */
    bool ProcessAnalysisParticle(const LArAnalysisParticle *const pAnalysisParticle, const MCPrimaryMap &mainMcParticleMap,
        const MCParticleList *const pMCParticleList) const;

    /**
     *  @brief  Record MC information for all the unreconstructed particles
     *
     *  @param  pMCParticleList address of the MC particle list
     *  @param  mainMcParticleMap the map from MC primaries to their analysis particles
     */
    void RecordUnreconstructedParticles(const MCParticleList *const pMCParticleList, const MCPrimaryMap &mainMcParticleMap) const;

    /**
     *  @brief  Record MC information for a given unreconstructed particle
     *
     *  @param  pMCParticleList address of the MC particle list
     *  @param  pMCPrimary address of the unreconstructed MC primary
     */
    void RecordUnreconstructedParticle(const MCParticleList *const pMCParticleList, const MCParticle *const pMCPrimary) const;

    /**
     *  @brief  Get all the MC primaries
     *
     *  @param  pMCParticleList address of the MC particle list
     *
     *  @return the MC primaries
     */
    MCParticleSet GetAllMcPrimaries(const MCParticleList *const pMCParticleList) const;

    /**
     *  @brief  Populate the tree parameters with neutrino information
     *
     *  @param  neutrinoAnalysisParticle the neutrino analysis particle
     *  @param  pMCParticleList address of the MC particle list
     */
    void PopulateNeutrinoParameters(const LArAnalysisParticle &neutrinoAnalysisParticle, const MCParticleList *const pMCParticleList) const;

    /**
     *  @brief  Recurse through the LArAnalysisParticle hierarchy and count the numbers of tracks and showers
     *
     *  @param  currentAnalysisParticle the current analysis particle
     *  @param  numberOfRecoTracks the number of reco tracks (to populate)
     *  @param  numberOfRecoShowers the number of reco showers (to populate)
     */
    void CountRecoTracksAndShowers(const LArAnalysisParticle &currentAnalysisParticle, unsigned int &numberOfRecoTracks,
        unsigned int &numberOfRecoShowers) const;

    /**
     *  @brief  Populate the tree parameters with neutrino MC information
     *
     *  @param  pfoMcInfo the PFO MC info object
     *  @param  pMCParticleList address of the MC particle list
     */
    void PopulateNeutrinoMcParameters(const LArAnalysisParticleHelper::PfoMcInfo &pfoMcInfo, const MCParticleList *const pMCParticleList) const;

    /**
     *  @brief  Caculate the neutrino MC visible energy and momentum
     *
     *  @param  pMCParticleList address of the MC particle list
     *  @param  visibleEnergy the visible energy (to populate)
     *  @param  visibleMomentum the visible momentum (to populate)
     */
    void CalculateNeutrinoMcVisibleMomentum(const MCParticleList *const pMCParticleList, float &visibleEnergy, CartesianVector &visibleMomentum) const;

    /**
     *  @brief  Get all the MC primary neutrino daughters
     *
     *  @param  pMCParticleList address of the MC particle list
     *
     *  @return the primary neutrino daughters
     */
    MCParticleSet GetAllMcPrimaryDaughters(const MCParticleList *const pMCParticleList) const;

    /**
     *  @brief  Get the interaction type for the event
     *
     *  @param  pMCParticleList address of the MC particle list
     *
     *  @return the interaction type
     */
    LArInteractionTypeHelper::InteractionType GetInteractionType(const MCParticleList *const pMCParticleList) const;

    /**
     *  @brief  Add a primary daughter record to the tree parameters
     *
     *  @param  primaryAnalysisParticle the primary daughter analysis particle
     *  @param  coveredMCPrimaries the list of MC primaries that have been covered so far
     */
    void AddPrimaryDaughterRecord(const LArAnalysisParticle &primaryAnalysisParticle, const MCPrimaryMap &coveredMCPrimaries) const;

    /**
     *  @brief  Add a primary daughter MC record
     *
     *  @param  pfoMcInfo the PFO MC info object
     *  @param  particleSplitByReco whether the particle has been split by the reco
     */
    void AddPrimaryDaughterMcRecord(const LArAnalysisParticleHelper::PfoMcInfo &pfoMcInfo, const bool particleSplitByReco) const;

    /**
     *  @brief  Add a blank primary daughter MC record
     */
    void AddBlankPrimaryDaughterMcRecord() const;

    /**
     *  @brief  Add an MC-only primary daughter record
     *
     *  @param  pfoMcInfo the PFO MC info object
     */
    void AddMcOnlyPrimaryDaughterRecord(const LArAnalysisParticleHelper::PfoMcInfo &pfoMcInfo) const;

    /**
     *  @brief  Add a cosmic ray record to the tree parameters
     *
     *  @param  cosmicRayAnalysisParticle the cosmic ray analysis particle
     *  @param  coveredMCPrimaries the list of MC primaries that have been covered so far
     */
    void AddCosmicRayRecord(const LArAnalysisParticle &cosmicRayAnalysisParticle, const MCPrimaryMap &coveredMCPrimaries) const;

    /**
     *  @brief  Add a cosmic ray MC record
     *
     *  @param  pfoMcInfo the PFO MC info object
     *  @param  particleSplitByReco whether the particle has been split by the reco
     */
    void AddCosmicRayMcRecord(const LArAnalysisParticleHelper::PfoMcInfo &pfoMcInfo, const bool particleSplitByReco) const;

    /**
     *  @brief  Add a blank cosmic ray MC record
     */
    void AddBlankCosmicRayMcRecord() const;

    /**
     *  @brief  Add an MC-only cosmic ray record
     *
     *  @param  pfoMcInfo the PFO MC info object
     */
    void AddMcOnlyCosmicRayRecord(const LArAnalysisParticleHelper::PfoMcInfo &pfoMcInfo) const;

    /**
     *  @brief  Check that the sizes of the tree vectors are consistent
     */
    void CheckTreeVectorSizes() const;

    /**
     *  @brief  Check the size of a given vector, throwing an error if it does not match
     *
     *  @param  variableName the variable name (for printing an error message)
     *  @param  vector the vector
     *  @param  desiredSize the desired vector size
     */
    template <typename T>
    void CheckSize(const std::string &variableName, const std::vector<T> &vector, const std::size_t desiredSize) const;

    /**
     *  @brief  Print the tree parameters
     */
    void PrintTree() const;

    /**
     *  @brief  Find out whether a given interaction type is charged-current
     *
     *  @param  interactionType the interaction type enum
     *
     *  @return whether the interaction type is charged-current
     */
    bool IsChargedCurrent(const LArInteractionTypeHelper::InteractionType interactionType) const;

    /**
     *  @brief  Initialize the analysis TTree
     */
    void InitializeTree();

    std::string               m_pfoListName;                        ///< The neutrino PFO list name
    std::string               m_outputFile;                         ///< The output file path
    std::string               m_treeName;                           ///< The name of the TTree
    std::string               m_treeTitle;                          ///< The title of the TTree
    TFile                    *m_pOutputTFile;                       ///< The ROOT TFile associated with the tree
    TTree                    *m_pOutputTree;                        ///< The ROOT TTree to which to write the data
    bool                      m_verbose;                            ///< Whether to print some AnalysisParticle information to screen
    mutable TreeParameters    m_treeParameters;                     ///< The tree parameters
    std::string               m_mcParticleListName;                 ///< The name of the MC particle list
    float                     m_fiducialCutLowXMargin;              ///< The low-x fiducial volume margin
    float                     m_fiducialCutHighXMargin;             ///< The high-x fiducial volume margin
    float                     m_fiducialCutLowYMargin;              ///< The low-y fiducial volume margin
    float                     m_fiducialCutHighYMargin;             ///< The high-y fiducial volume margin
    float                     m_fiducialCutLowZMargin;              ///< The low-z fiducial volume margin
    float                     m_fiducialCutHighZMargin;             ///< The high-z fiducial volume margin
    CartesianVector           m_minCoordinates;                     ///< The set of detector minimum coordinates
    CartesianVector           m_maxCoordinates;                     ///< The set of detector maximum coordinates
    float                     m_mcContainmentFractionLowerBound;    ///< The lower containment fraction bound for MC containment
    float                     m_fiducialHitFractionLowerBound;      ///< The lower fiducial hit fraction bound for containment
    float                     m_mcOnlyParticleContainmentCut;       ///< The lower containment fraction bound for including MC-only particles
    float                     m_mcOnlyParticleEnergyCut;            ///< The lower energy bound for including MC-only particles
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline void WriteAnalysisParticlesAlgorithm::CheckSize(const std::string &variableName, const std::vector<T> &vector,
    const std::size_t desiredSize) const
{
    if (vector.size() != desiredSize)
    {
        std::cout << "WriteAnalysisParticlesAlgorithm: tree member variable '" <<  variableName << "' was of the wrong size" << std::endl;
        throw STATUS_CODE_FAILURE;
    }
}

} // namespace lar_physics_content

#endif // #ifndef LAR_WRITE_ANALYSIS_PARTICLES_ALGORITHM_H
