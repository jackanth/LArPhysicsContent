/**
 *  @file   LArPhysicsContent/include/WriteAnalysisParticlesAlgorithm.h
 *
 *  @brief  Header file for the write AnalysisParticles algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_WRITE_ANALYSIS_PARTICLES_ALGORITHM_H
#define LAR_WRITE_ANALYSIS_PARTICLES_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArAnalysisParticle.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

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
     *  @brief  Constructor
     */
    TreeParameters() noexcept;

    Bool_t             m_nu_WasReconstructed;                                 ///< Whether the neutrino has been reconstructed
    Bool_t             m_nu_IsVertexFiducial;                                 ///< Whether the neutrino vertex is fiducial
    Bool_t             m_nu_IsContained;                                      ///< Whether the neutrino looks contained
    Float_t            m_nu_FiducialHitFraction;                              ///< The fraction of the neutrino's hits that are fiducial
    Bool_t             m_nu_HasMcInfo;                                        ///< Whether the neutrino has MC information
    Float_t            m_nu_VisibleEnergy;                                    ///< The neutrino's reconstructed energy in GeV
    Float_t            m_nu_VisibleLongitudinalEnergy;                        ///< The neutrino's reconstructed transverse energy in GeV
    Float_t            m_nu_VisibleTransverseEnergy;                          ///< The neutrino's reconstructed longitudinal energy in GeV
    Float_t            m_nu_VisibleEnergyFracFromRange;                       ///< The fraction of the reconstructed neutrino energy calculated from range
    Float_t            m_nu_VisibleEnergyFracFromCorrectedTrackCharge;        ///< The fraction of the reconstructed neutrino energy calculated from recombination-corrected track charge
    Float_t            m_nu_VisibleEnergyFracFromUncorrectedTrackCharge;      ///< The fraction of the reconstructed neutrino energy calculated from uncorrected track charge
    Float_t            m_nu_VisibleEnergyFracFromShowerCharge;                ///< The fraction of the reconstructed neutrino energy calculated from shower charge
    Float_t            m_nu_VertexX;                                          ///< The x-component of the neutrino vertex position in cm
    Float_t            m_nu_VertexY;                                          ///< The y-component of the neutrino vertex position in cm
    Float_t            m_nu_VertexZ;                                          ///< The z-component of the neutrino vertex position in cm
    Float_t            m_nu_DirectionCosineX;                                 ///< The direction cosine of the neutrino in the x-direction
    Float_t            m_nu_DirectionCosineY;                                 ///< The direction cosine of the neutrino in the y-direction
    Float_t            m_nu_DirectionCosineZ;                                 ///< The direction cosine of the neutrino in the z-direction
    TString            m_nu_TypeTree;                                         ///< The neutrino type tree
    UInt_t             m_nu_NumberOf3dHits;                                   ///< The number of 3D hits in the neutrino
    UInt_t             m_nu_NumberOfCollectionPlaneHits;                      ///< The number of collection-plane hits in the neutrino
    UInt_t             m_nu_NumberOfDownstreamParticles;                      ///< The number of particles downstream of the neutrino
    ULong64_t          m_nu_mc_McParticleUid;                                 ///< The UID of the MC particle
    Float_t            m_nu_mc_Energy;                                        ///< The MC energy of the neutrino in GeV
    Float_t            m_nu_mc_LongitudinalEnergy;                            ///< The MC longitudinal energy of the neutrino in GeV
    Float_t            m_nu_mc_TransverseEnergy;                              ///< The MC transverse energy of the neutrino in GeV
    Float_t            m_nu_mc_VisibleEnergy;                                 ///< The MC visible energy of the neutrino in GeV
    Float_t            m_nu_mc_VisibleLongitudinalEnergy;                     ///< The MC visible longitudinal energy of the neutrino in GeV
    Float_t            m_nu_mc_VisibleTransverseEnergy;                       ///< The MC visible transverse energy of the neutrino in GeV
    Float_t            m_nu_mc_VertexX;                                       ///< The MC x-component of the neutrino vertex position in cm
    Float_t            m_nu_mc_VertexY;                                       ///< The MC y-component of the neutrino vertex position in cm
    Float_t            m_nu_mc_VertexZ;                                       ///< The MC x-component of the neutrino vertex position in cm
    Float_t            m_nu_mc_DirectionCosineX;                              ///< The MC direction cosine of the neutrino in the x-direction
    Float_t            m_nu_mc_DirectionCosineY;                              ///< The MC direction cosine of the neutrino in the y-direction
    Float_t            m_nu_mc_DirectionCosineZ;                              ///< The MC direction cosine of the neutrino in the z-direction
    Float_t            m_nu_mc_Momentum;                                      ///< The MC momentum of the neutrino in GeV/c
    Float_t            m_nu_mc_MomentumX;                                     ///< The MC momentum of the neutrino in the x-direction in GeV/c
    Float_t            m_nu_mc_MomentumY;                                     ///< The MC momentum of the neutrino in the y-direction in GeV/c
    Float_t            m_nu_mc_MomentumZ;                                     ///< The MC momentum of the neutrino in the z-direction in GeV/c
    Bool_t             m_nu_mc_IsVertexFiducial;                              ///< Whether the neutrino vertex is fiducial (MC quantity)
    Bool_t             m_nu_mc_IsContained;                                   ///< Whether the neutrino is contained (MC quantity)
    Float_t            m_nu_mc_ContainmentFraction;                           ///< The fraction of the neutrino that is contained (MC quantity)
    TString            m_nu_mc_TypeTree;                                      ///< The neutrino's MC type tree
    TString            m_nu_mc_InteractionType;                               ///< The neutrino interaction type (MC quantity)
    Bool_t             m_nu_mc_IsChargedCurrent;                              ///< Whether the interaction is charged-current (MC quantity)
    Float_t            m_nu_mc_VisibleEnergyFraction;                         ///< The fraction of the neutrino's energy that is visible (MC quantity)
    Int_t              m_nu_mc_PdgCode;                                       ///< The neutrino's PDG code (MC quantity)
    Float_t            m_nu_mc_HitPurity;                                     ///< The neutrino's hit number purity (MC quantity)
    Float_t            m_nu_mc_HitCompleteness;                               ///< The neutrino's hit number completeness (MC quantity)
    Float_t            m_nu_mc_CollectionPlaneHitPurity;                      ///< The neutrino's hit number purity in the collection plane (MC quantity)
    Float_t            m_nu_mc_CollectionPlaneHitCompleteness;                ///< The neutrino's hit number completeness in the collection plane (MC quantity)
                                                                      
    UInt_t             m_primary_Number;                                      ///< The number of primary daughters
    RBoolVector        m_primary_WasReconstructed;                            ///< Whether each primary daughter has been reconstructed
    RBoolVector        m_primary_IsVertexFiducial;                            ///< Whether each primary daughter's vertex is fiducial
    RBoolVector        m_primary_IsContained;                                 ///< Whether the primary looks contained
    RFloatVector       m_primary_FiducialHitFraction;                         ///< The fraction of the primary's hits that are fiducial
    RBoolVector        m_primary_HasMcInfo;                                   ///< Whether each primary daughter has MC information
    RFloatVector       m_primary_KineticEnergy;                               ///< The reconstructed kinetic energy of each primary daughter in GeV
    RFloatVector       m_primary_KineticEnergyFracFromRange;                  ///< The fraction of each primary's reconstructed kinetic energy calculated from range
    RFloatVector       m_primary_KineticEnergyFracFromCorrectedTrackCharge;   ///< The fraction of each primary's reconstructed kinetic energy calculated from recombination-corrected track charge
    RFloatVector       m_primary_KineticEnergyFracFromUncorrectedTrackCharge; ///< The fraction of each primary's reconstructed kinetic energy calculated from uncorrected track charge
    RFloatVector       m_primary_KineticEnergyFracFromShowerCharge;           ///< The fraction of each primary's reconstructed kinetic energy calculated from shower charge
    RFloatVector       m_primary_VertexX;                                     ///< The x-component of each primary daughter's vertex position in cm
    RFloatVector       m_primary_VertexY;                                     ///< The y-component of each primary daughter's vertex position in cm
    RFloatVector       m_primary_VertexZ;                                     ///< The z-component of each primary daughter's vertex position in cm
    RFloatVector       m_primary_DirectionCosineX;                            ///< The direction cosine of each primary daughter in the x-direction
    RFloatVector       m_primary_DirectionCosineY;                            ///< The direction cosine of each primary daughter in the y-direction
    RFloatVector       m_primary_DirectionCosineZ;                            ///< The direction cosine of each primary daughter in the z-direction
    RTextVector        m_primary_TypeTree;                                    ///< The type tree of each primary daughter
    RBoolVector        m_primary_IsShower;                                    ///< Whether each primary daughter is a shower
    RBoolVector        m_primary_IsTrack;                                     ///< Whether each primary daughter is a track
    RBoolVector        m_primary_IsProton;                                    ///< Whether each primary daughter is a proton
    RBoolVector        m_primary_IsPionOrMuon;                                ///< Whether each primary daughter is a pion or muon
    RUnsignedVector    m_primary_NumberOf3dHits;                              ///< The number of 3D hits in each primary daughter
    RUnsignedVector    m_primary_NumberOfCollectionPlaneHits;                 ///< The number of collection-plane hits in each primary daughter
    RUnsignedVector    m_primary_NumberOfDownstreamParticles;                 ///< The number of particles downstream of each primary daughter
    RUInt64Vector      m_primary_mc_McParticleUid;                            ///< The UID of the MC particle corresponding to each primary
    RBoolVector        m_primary_mc_IsParticleSplitByReco;                    ///< Whether the primary daughter's MC particle has been split by the reconstruction
    RFloatVector       m_primary_mc_Energy;                                   ///< The MC energy of each primary daughter in GeV
    RFloatVector       m_primary_mc_KineticEnergy;                            ///< The MC kinetic energy of each primary daughter in GeV
    RFloatVector       m_primary_mc_Mass;                                     ///< The MC mass of each primary daughter in GeV/c^2
    RFloatVector       m_primary_mc_VertexX;                                  ///< The MC x-component of each primary daughter's vertex position in cm
    RFloatVector       m_primary_mc_VertexY;                                  ///< The MC y-component of each primary daughter's vertex position in cm
    RFloatVector       m_primary_mc_VertexZ;                                  ///< The MC z-component of each primary daughter's vertex position inc m
    RFloatVector       m_primary_mc_DirectionCosineX;                         ///< The MC direction cosine of each primary daughter in the x-direction
    RFloatVector       m_primary_mc_DirectionCosineY;                         ///< The MC direction cosine of each primary daughter in the y-direction
    RFloatVector       m_primary_mc_DirectionCosineZ;                         ///< The MC direction cosine of each primary daughter in the z-direction
    RFloatVector       m_primary_mc_Momentum;                                 ///< The MC momentum of each primary daughter in GeV/c
    RFloatVector       m_primary_mc_MomentumX;                                ///< The MC momentum of each primary daughter in the x-direction in GeV/c
    RFloatVector       m_primary_mc_MomentumY;                                ///< The MC momentum of each primary daughter in the y-direction in GeV/c
    RFloatVector       m_primary_mc_MomentumZ;                                ///< The MC momentum of each primary daughter in the z-direction in GeV/c
    RBoolVector        m_primary_mc_IsVertexFiducial;                         ///< Whether each primary daughter's vertex is fiducial (MC quantity)
    RBoolVector        m_primary_mc_IsContained;                              ///< Whether each primary daughter is contained (MC quantity)
    RFloatVector       m_primary_mc_ContainmentFraction;                      ///< The fraction of the primary that is contained (MC quantity)
    RTextVector        m_primary_mc_TypeTree;                                 ///< The MC type tree for each primary daughter
    RBoolVector        m_primary_mc_IsShower;                                 ///< Whether each primary daughter is a shower (MC quantity)
    RBoolVector        m_primary_mc_IsTrack;                                  ///< Whether each primary daughter is a track (MC quantity)
    RBoolVector        m_primary_mc_IsProton;                                 ///< Whether each primary daughter is a proton (MC quantity)
    RBoolVector        m_primary_mc_IsPionOrMuon;                             ///< Whether each primary daughter is a pion or muon (MC quantity)
    RBoolVector        m_primary_mc_IsCosmicRay;                              ///< Whether each primary daughter is a cosmic ray (MC quantity)
    RIntVector         m_primary_mc_PdgCode;                                  ///< The primary daughter's PDG code (MC quantity)
    RFloatVector       m_primary_mc_HitPurity;                                ///< The primary daughter's hit number purity (MC quantity)
    RFloatVector       m_primary_mc_HitCompleteness;                          ///< The primary daughter's hit number completeness (MC quantity)
    RFloatVector       m_primary_mc_CollectionPlaneHitPurity;                 ///< The primary daughter's hit number purity in the collection plane (MC quantity)
    RFloatVector       m_primary_mc_CollectionPlaneHitCompleteness;           ///< The primary daughter's hit number completeness in the collection plane (MC quantity)
                                                                  
    UInt_t             m_cr_Number;                                           ///< The number of cosmic rays
    RBoolVector        m_cr_WasReconstructed;                                 ///< Whether each cosmic ray has been reconstructed
    RBoolVector        m_cr_IsVertexFiducial;                                 ///< Whether each cosmic ray's vertex is fiducial
    RBoolVector        m_cr_IsContained;                                      ///< Whether the cosmic ray looks contained
    RFloatVector       m_cr_FiducialHitFraction;                              ///< The fraction of the cosmic ray's hits that are fiducial
    RBoolVector        m_cr_HasMcInfo;                                        ///< Whether each cosmic ray has MC information
    RFloatVector       m_cr_KineticEnergy;                                    ///< The reconstructed kinetic energy of each cosmic ray in GeV
    RFloatVector       m_cr_KineticEnergyFracFromRange;                       ///< The fraction of each cosmic ray's reconstructed kinetic energy calculated from range
    RFloatVector       m_cr_KineticEnergyFracFromCorrectedTrackCharge;        ///< The fraction of each cosmic ray's reconstructed kinetic energy calculated from recombination-corrected track charge
    RFloatVector       m_cr_KineticEnergyFracFromUncorrectedTrackCharge;      ///< The fraction of each cosmic ray's reconstructed kinetic energy calculated from uncorrected track charge
    RFloatVector       m_cr_KineticEnergyFracFromShowerCharge;                ///< The fraction of each cosmic ray's reconstructed kinetic energy calculated from shower charge
    RFloatVector       m_cr_VertexX;                                          ///< The x-component of each cosmic ray's vertex position in cm
    RFloatVector       m_cr_VertexY;                                          ///< The y-component of each cosmic ray's vertex position in cm
    RFloatVector       m_cr_VertexZ;                                          ///< The z-component of each cosmic ray's vertex position in cm
    RFloatVector       m_cr_DirectionCosineX;                                 ///< The direction cosine of each cosmic ray in the x-direction
    RFloatVector       m_cr_DirectionCosineY;                                 ///< The direction cosine of each cosmic ray in the y-direction
    RFloatVector       m_cr_DirectionCosineZ;                                 ///< The direction cosine of each cosmic ray in the z-direction
    RTextVector        m_cr_TypeTree;                                         ///< The type tree of each cosmic ray
    RUnsignedVector    m_cr_NumberOf3dHits;                                   ///< The number of 3D hits in each cosmic ray
    RUnsignedVector    m_cr_NumberOfCollectionPlaneHits;                      ///< The number of collection-plane hits in each cosmic ray
    RUnsignedVector    m_cr_NumberOfDownstreamParticles;                      ///< The number of particles downstream of each cosmic ray
    RUInt64Vector      m_cr_mc_McParticleUid;                                 ///< The UID of the MC particle corresponding to each cosmic ray
    RBoolVector        m_cr_mc_IsParticleSplitByReco;                         ///< Whether the cosmic ray's MC particle has been split by the reconstruction
    RFloatVector       m_cr_mc_Energy;                                        ///< The MC energy of each cosmic ray in GeV
    RFloatVector       m_cr_mc_KineticEnergy;                                 ///< The MC kinetic energy of each cosmic ray in GeV
    RFloatVector       m_cr_mc_Mass;                                          ///< The MC mass of each cosmic ray in GeV/c^2
    RFloatVector       m_cr_mc_VertexX;                                       ///< The MC x-component of each cosmic ray's vertex position in cm
    RFloatVector       m_cr_mc_VertexY;                                       ///< The MC y-component of each cosmic ray's vertex position in cm
    RFloatVector       m_cr_mc_VertexZ;                                       ///< The MC z-component of each cosmic ray's vertex position in cm
    RFloatVector       m_cr_mc_DirectionCosineX;                              ///< The MC direction cosine of each cosmic ray in the x-direction
    RFloatVector       m_cr_mc_DirectionCosineY;                              ///< The MC direction cosine of each cosmic ray in the y-direction
    RFloatVector       m_cr_mc_DirectionCosineZ;                              ///< The MC direction cosine of each cosmic ray in the z-direction
    RFloatVector       m_cr_mc_Momentum;                                      ///< The MC momentum of each cosmic ray in GeV/c
    RFloatVector       m_cr_mc_MomentumX;                                     ///< The MC momentum of each cosmic ray in the x-direction in GeV/c
    RFloatVector       m_cr_mc_MomentumY;                                     ///< The MC momentum of each cosmic ray in the y-direction in GeV/c
    RFloatVector       m_cr_mc_MomentumZ;                                     ///< The MC momentum of each cosmic ray in the z-direction in GeV/c
    RBoolVector        m_cr_mc_IsVertexFiducial;                              ///< Whether each cosmic ray's vertex is fiducial (MC quantity)
    RBoolVector        m_cr_mc_IsContained;                                   ///< Whether each cosmic ray is contained (MC quantity)
    RFloatVector       m_cr_mc_ContainmentFraction;                           ///< The fraction of the CR that is contained (MC quantity)
    RTextVector        m_cr_mc_TypeTree;                                      ///< The MC type tree for each cosmic ray
    RBoolVector        m_cr_mc_IsShower;                                      ///< Whether each cosmic ray is a shower (MC quantity)
    RBoolVector        m_cr_mc_IsTrack;                                       ///< Whether each cosmic ray is a track (MC quantity)
    RBoolVector        m_cr_mc_IsProton;                                      ///< Whether each cosmic ray is a proton (MC quantity)
    RBoolVector        m_cr_mc_IsPionOrMuon;                                  ///< Whether each cosmic ray is a pion or muon (MC quantity)
    RBoolVector        m_cr_mc_IsCosmicRay;                                   ///< Whether each cosmic ray is a cosmic ray (MC quantity)
    RIntVector         m_cr_mc_PdgCode;                                       ///< The cosmic ray's PDG code (MC quantity)
    RFloatVector       m_cr_mc_HitPurity;                                     ///< The cosmic ray's hit number purity (MC quantity)
    RFloatVector       m_cr_mc_HitCompleteness;                               ///< The cosmic ray's hit number completeness (MC quantity)
    RFloatVector       m_cr_mc_CollectionPlaneHitPurity;                      ///< The cosmic ray's hit number purity in the collection plane (MC quantity)
    RFloatVector       m_cr_mc_CollectionPlaneHitCompleteness;                ///< The cosmic ray's hit number completeness in the collection plane (MC quantity)
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
     *  @brief  Constructor
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
     *  @brief  Populate the tree parameters with neutrino information
     * 
     *  @param  neutrinoAnalysisParticle the neutrino analysis particle
     *  @param  pMCParticleList address of the MC particle list
     *  @param  pCaloHitList address of the CaloHitList
     */
    void PopulateNeutrinoParameters(const LArAnalysisParticle &neutrinoAnalysisParticle, const MCParticleList *const pMCParticleList,
        const CaloHitList *const pCaloHitList) const;
    
    /**
     *  @brief  Populate the tree parameters with neutrino MC information
     * 
     *  @param  pMainMcParticle address of the main MC particle
     *  @param  mcEnergy the MC energy
     *  @param  mcVertexPosition the MC vertex position
     *  @param  mcDirectionCosines the MC direction cosines
     *  @param  mcMomentum the MC momentum
     *  @param  mcIsVertexFiducial whether the MC vertex is fiducial
     *  @param  mcContainmentFraction the MC containment fraction
     *  @param  mcPdgCode the MC PDG code
     *  @param  pMCParticleList address of the MC particle list
     *  @param  pCaloHitList address of the CaloHit list
     *  @param  mcHitPurity ahe MC hit purity
     *  @param  mcHitCompleteness the MC hit completeness
     *  @param  mcCollectionPlaneHitPurity the MC collection plane hit purity
     *  @param  mcCollectionPlaneHitCompleteness the MC collection plane hit completeness
     *  @param  mcTypeTree the MC type tree
     */
    void PopulateNeutrinoMcParameters(const MCParticle *const pMainMcParticle, const float mcEnergy,
        const CartesianVector &mcVertexPosition, const CartesianVector &mcDirectionCosines, const CartesianVector &mcMomentum, 
        const bool mcIsVertexFiducial, const float mcContainmentFraction, const int mcPdgCode, const MCParticleList *const pMCParticleList,
        const CaloHitList *const pCaloHitList, const float mcHitPurity, const float mcHitCompleteness, 
        const float mcCollectionPlaneHitPurity, const float mcCollectionPlaneHitCompleteness, 
        const LArAnalysisParticle::TypeTree mcTypeTree) const;
    
    /**
     *  @brief  Get the interaction type for the event
     * 
     *  @param  pMCParticleList address of the MC particle list
     *  @param  pCaloHitList address of the CaloHit list
     * 
     *  @return the interaction type
     */
    LArInteractionTypeHelper::InteractionType GetInteractionType(const MCParticleList *const pMCParticleList,
        const CaloHitList *const pCaloHitList) const;
    
    /**
     *  @brief  Add a primary daughter record to the tree parameters
     * 
     *  @param  primaryAnalysisParticle the primary daughter analysis particle
     *  @param  coveredMCPrimaries the list of MC primaries that have been covered so far
     */
    void AddPrimaryDaughterRecord(const LArAnalysisParticle &primaryAnalysisParticle, const MCPrimaryMap &coveredMCPrimaries) const;
    
    /**
     *  @brief  Add an MC-only primary daughter record
     * 
     *  @param  pMainMcParticle address of the main MC particle
     *  @param  mcEnergy the MC energy
     *  @param  mcKineticEnergy the MC kinetic energy
     *  @param  mcMass the MC mass
     *  @param  mcVertexPosition the MC vertex position
     *  @param  mcDirectionCosines the MC direction cosines
     *  @param  mcMomentum the MC momentum
     *  @param  mcIsVertexFiducial whether the MC vertex is fiducial
     *  @param  mcContainmentFraction the MC containment fraction
     *  @param  mcType the MC type
     *  @param  mcIsShower whether the MC particle is a shower
     *  @param  mcPdgCode the MC PDG code
     *  @param  mcTypeTree the MC type tree
     */
    void AddMcOnlyPrimaryDaughterRecord(const MCParticle *const pMainMcParticle, const float mcEnergy, const float mcKineticEnergy, 
        const float mcMass, const CartesianVector &mcVertexPosition, const CartesianVector &mcDirectionCosines, 
        const CartesianVector &mcMomentum, const bool mcIsVertexFiducial, const float mcContainmentFraction, 
        const LArAnalysisParticle::TYPE mcType, const bool mcIsShower, const int mcPdgCode, 
        const LArAnalysisParticle::TypeTree mcTypeTree) const;
        
    /**
     *  @brief  Add a cosmic ray record to the tree parameters
     * 
     *  @param  cosmicRayAnalysisParticle the cosmic ray analysis particle
     *  @param  coveredMCPrimaries the list of MC primaries that have been covered so far
     */
    void AddCosmicRayRecord(const LArAnalysisParticle &cosmicRayAnalysisParticle, const MCPrimaryMap &coveredMCPrimaries) const;
    
    /**
     *  @brief  Add an MC-only cosmic ray record
     * 
     *  @param pMainMcParticle address of the main MC particle
     *  @param mcEnergy the MC energy
     *  @param mcKineticEnergy the MC kinetic energy
     *  @param mcMass the MC mass
     *  @param mcVertexPosition the MC vertex position
     *  @param mcDirectionCosines the MC direction cosines
     *  @param mcMomentum the MC momentum
     *  @param mcIsVertexFiducial whether the MC vertex is fiducial
     *  @param mcContainmentFraction the MC containment fraction
     *  @param mcType the MC type
     *  @param mcIsShower whether the MC particle is a shower
     *  @param mcPdgCode the MC PDG code
     *  @param mcTypeTree the MC type tree
     */
    void AddMcOnlyCosmicRayRecord(const MCParticle *const pMainMcParticle, const float mcEnergy, const float mcKineticEnergy, 
        const float mcMass, const CartesianVector &mcVertexPosition, const CartesianVector &mcDirectionCosines, 
        const CartesianVector &mcMomentum, const bool mcIsVertexFiducial, const float mcContainmentFraction, 
        const LArAnalysisParticle::TYPE mcType, const bool mcIsShower, const int mcPdgCode, 
        const LArAnalysisParticle::TypeTree mcTypeTree) const;
    
    /**
     *  @brief  Dump the tree parameters
     */
    void DumpTree() const;
        
    /**
     *  @brief  Find out whether a given interaction type is charged-current
     * 
     *  @param  interactionType the interaction type enum
     * 
     *  @return whether the interaction type is charged-current
     */
    bool IsChargedCurrent(const LArInteractionTypeHelper::InteractionType interactionType) const;

    std::string               m_pfoListName;                     ///< The neutrino PFO list name
    std::string               m_outputFile;                      ///< The output file path
    TFile                    *m_pOutputTFile;                    ///< The ROOT TFile associated with the tree
    TTree                    *m_pOutputTree;                     ///< The ROOT TTree to which to write the data
    bool                      m_verbose;                         ///< Whether to print some AnalysisParticle information to screen
    mutable TreeParameters    m_treeParameters;                  ///< The tree parameters
    std::string               m_mcParticleListName;              ///< The name of the MC particle list
    std::string               m_caloHitListName;                 ///< The name of the CaloHit list
    float                     m_fiducialCutLowXMargin;           ///< The low-X fiducial volume margin
    float                     m_fiducialCutHighXMargin;          ///< The high-X fiducial volume margin
    float                     m_fiducialCutLowYMargin;           ///< The low-Y fiducial volume margin
    float                     m_fiducialCutHighYMargin;          ///< The high-Y fiducial volume margin
    float                     m_fiducialCutLowZMargin;           ///< The low-Z fiducial volume margin
    float                     m_fiducialCutHighZMargin;          ///< The high-Z fiducial volume margin
    CartesianVector           m_minCoordinates;                  ///< The set of detector minimum coordinates
    CartesianVector           m_maxCoordinates;                  ///< The set of detector maximum coordinates
    float                     m_mcContainmentFractionLowerBound; ///< The lower containment fraction bound for MC containment
    float                     m_fiducialHitFractionLowerBound;   ///< The lower fiducial hit fraction bound for containment
    float                     m_mcOnlyParticleContainmentCut;    ///< The lower containment fraction bound for including MC-only particles
    float                     m_mcOnlyParticleEnergyCut;         ///< The lower energy bound for including MC-only particles
};

} // namespace lar_physics_content

#endif // #ifndef LAR_WRITE_ANALYSIS_PARTICLES_ALGORITHM_H
