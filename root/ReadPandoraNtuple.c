#include "Common.h"

//------------------------------------------------------------------------------------------------------------------------------------------

#define QUICK_PLOT_1D(variableName, ntupleAddress)                                                        \
{                                                                                                         \
    TString identifier = variableName;                                                                    \
    struct PlotSettings1D plotSettings = g_defaultPlotSettings1D;                                         \
    PlotNtuple1D(ntupleAddress, identifier, identifier, plotSettings)->SaveAs(identifier.Append(".png")); \
}

#define QUICK_PLOT_2D(variableName1, variableName2, ntupleAddress)                                                      \
{                                                                                                                       \
    TString identifier1 = variableName1;                                                                                \
    TString identifier2 = variableName2;                                                                                \
    TString identifier = identifier1 + "_" + identifier2;                                                               \
    struct PlotSettings2D plotSettings = g_defaultPlotSettings2D;                                                       \
    PlotNtuple2D(ntupleAddress, identifier1, identifier2, identifier, plotSettings)->SaveAs(identifier.Append(".png")); \
}
   
#define PLOT_2D(variableName1, variableName2, ntupleAddress, minX, maxX, minY, maxY)                                    \
{                                                                                                                       \
    TString identifier1 = variableName1;                                                                                \
    TString identifier2 = variableName2;                                                                                \
    TString identifier = identifier1 + "_" + identifier2;                                                               \
    struct PlotSettings2D plotSettings = g_defaultPlotSettings2D;                                                       \
    plotSettings.xMin = minX;                                                                                           \
    plotSettings.xMax = maxX;                                                                                           \
    plotSettings.yMin = minY;                                                                                           \
    plotSettings.yMax = maxY;                                                                                           \
    PlotNtuple2D(ntupleAddress, identifier1, identifier2, identifier, plotSettings)->SaveAs(identifier.Append(".png")); \
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ReadPandoraNtuple(const char *const inputFilePath)
{
    typedef std::vector<Bool_t>          BoolVector;
    typedef std::vector<Int_t>           IntVector;
    typedef std::vector<Float_t>         FloatVector;
    typedef std::vector<ULong64_t> UInt64Vector;
    typedef std::vector<TString>   StringVector;
    typedef std::vector<UInt_t>      UnsignedVector;
    
    TTree *const pTree = LoadTreeFromFile(inputFilePath, "PandoraTree");
    
    Bool_t          nu_WasReconstructed                                    = false;
    Bool_t          nu_IsVertexFiducial                                    = false;
    Bool_t          nu_IsContained                                         = false;
    Float_t         nu_FiducialHitFraction                                 = 0.f;
    Bool_t          nu_HasMcInfo                                           = false;
    Float_t         nu_VisibleEnergy                                       = 0.f;
    Float_t         nu_VisibleLongitudinalEnergy                           = 0.f;
    Float_t         nu_VisibleTransverseEnergy                             = 0.f;
    Float_t         nu_VisibleEnergyFracFromRange                          = 0.f;
    Float_t         nu_VisibleEnergyFracFromCorrectedTrackCharge           = 0.f;
    Float_t         nu_VisibleEnergyFracFromUncorrectedTrackCharge         = 0.f;
    Float_t         nu_VisibleEnergyFracFromShowerCharge                   = 0.f;
    Float_t         nu_VertexX                                             = 0.f;
    Float_t         nu_VertexY                                             = 0.f;
    Float_t         nu_VertexZ                                             = 0.f;
    Float_t         nu_DirectionCosineX                                    = 0.f;
    Float_t         nu_DirectionCosineY                                    = 0.f;
    Float_t         nu_DirectionCosineZ                                    = 0.f;
    TString        *p_nu_TypeTree                                          = nullptr;
    UInt_t          nu_NumberOf3dHits                                      = 0U;
    UInt_t          nu_NumberOfCollectionPlaneHits                         = 0U;
    UInt_t          nu_NumberOfDownstreamParticles                         = 0U;
    Long64_t        nu_mc_McParticleUid                                    = 0LL;
    Float_t         nu_mc_Energy                                           = 0.f;
    Float_t         nu_mc_LongitudinalEnergy                               = 0.f;
    Float_t         nu_mc_TransverseEnergy                                 = 0.f;
    Float_t         nu_mc_VisibleEnergy                                    = 0.f;
    Float_t         nu_mc_VisibleLongitudinalEnergy                        = 0.f;
    Float_t         nu_mc_VisibleTransverseEnergy                          = 0.f;
    Float_t         nu_mc_VertexX                                          = 0.f;
    Float_t         nu_mc_VertexY                                          = 0.f;
    Float_t         nu_mc_VertexZ                                          = 0.f;
    Float_t         nu_mc_DirectionCosineX                                 = 0.f;
    Float_t         nu_mc_DirectionCosineY                                 = 0.f;
    Float_t         nu_mc_DirectionCosineZ                                 = 0.f;
    Float_t         nu_mc_Momentum                                         = 0.f;
    Float_t         nu_mc_MomentumX                                        = 0.f;
    Float_t         nu_mc_MomentumY                                        = 0.f;
    Float_t         nu_mc_MomentumZ                                        = 0.f;
    Bool_t          nu_mc_IsVertexFiducial                                 = false;
    Bool_t          nu_mc_IsContained                                      = false;
    Float_t         nu_mc_ContainmentFraction                              = 0.f;
    TString        *p_nu_mc_TypeTree                                       = nullptr;
    TString        *p_nu_mc_InteractionType                                = nullptr;
    Bool_t          nu_mc_IsChargedCurrent                                 = false;
    Float_t         nu_mc_VisibleEnergyFraction                            = 0.f;
    Int_t           nu_mc_PdgCode                                          = 0;
    Float_t         nu_mc_HitPurity                                        = 0.f;
    Float_t         nu_mc_HitCompleteness                                  = 0.f;
    Float_t         nu_mc_CollectionPlaneHitPurity                         = 0.f;
    Float_t         nu_mc_CollectionPlaneHitCompleteness                   = 0.f;
    UInt_t          primary_Number                                         = 0U;
    BoolVector     *p_primary_WasReconstructed                             = nullptr;
    BoolVector     *p_primary_IsVertexFiducial                             = nullptr;
    BoolVector     *p_primary_IsContained                                  = nullptr;
    FloatVector    *p_primary_FiducialHitFraction                          = nullptr;
    BoolVector     *p_primary_HasMcInfo                                    = nullptr;
    FloatVector    *p_primary_KineticEnergy                                = nullptr;
    FloatVector    *p_primary_KineticEnergyFracFromRange                   = nullptr;
    FloatVector    *p_primary_KineticEnergyFracFromCorrectedTrackCharge    = nullptr;
    FloatVector    *p_primary_KineticEnergyFracFromUncorrectedTrackCharge  = nullptr;
    FloatVector    *p_primary_KineticEnergyFracFromShowerCharge            = nullptr;
    FloatVector    *p_primary_VertexX                                      = nullptr;
    FloatVector    *p_primary_VertexY                                      = nullptr;
    FloatVector    *p_primary_VertexZ                                      = nullptr;
    FloatVector    *p_primary_DirectionCosineX                             = nullptr;
    FloatVector    *p_primary_DirectionCosineY                             = nullptr;
    FloatVector    *p_primary_DirectionCosineZ                             = nullptr;
    StringVector   *p_primary_TypeTree                                     = nullptr;
    BoolVector     *p_primary_IsShower                                     = nullptr;
    BoolVector     *p_primary_IsTrack                                      = nullptr;
    BoolVector     *p_primary_IsProton                                     = nullptr;
    BoolVector     *p_primary_IsPionOrMuon                                 = nullptr;
    UnsignedVector *p_primary_NumberOf3dHits                               = nullptr;
    UnsignedVector *p_primary_NumberOfCollectionPlaneHits                  = nullptr;
    UnsignedVector *p_primary_NumberOfDownstreamParticles                  = nullptr;
    UInt64Vector   *p_primary_mc_McParticleUid                             = nullptr;
    BoolVector     *p_primary_mc_IsParticleSplitByReco                     = nullptr;
    FloatVector    *p_primary_mc_Energy                                    = nullptr;
    FloatVector    *p_primary_mc_KineticEnergy                             = nullptr;
    FloatVector    *p_primary_mc_Mass                                      = nullptr;
    FloatVector    *p_primary_mc_VertexX                                   = nullptr;
    FloatVector    *p_primary_mc_VertexY                                   = nullptr;
    FloatVector    *p_primary_mc_VertexZ                                   = nullptr;
    FloatVector    *p_primary_mc_DirectionCosineX                          = nullptr;
    FloatVector    *p_primary_mc_DirectionCosineY                          = nullptr;
    FloatVector    *p_primary_mc_DirectionCosineZ                          = nullptr;
    FloatVector    *p_primary_mc_Momentum                                  = nullptr;
    FloatVector    *p_primary_mc_MomentumX                                 = nullptr;
    FloatVector    *p_primary_mc_MomentumY                                 = nullptr;
    FloatVector    *p_primary_mc_MomentumZ                                 = nullptr;
    BoolVector     *p_primary_mc_IsVertexFiducial                          = nullptr;
    BoolVector     *p_primary_mc_IsContained                               = nullptr;
    FloatVector    *p_primary_mc_ContainmentFraction                       = nullptr;
    StringVector   *p_primary_mc_TypeTree                                  = nullptr;
    BoolVector     *p_primary_mc_IsShower                                  = nullptr;
    BoolVector     *p_primary_mc_IsTrack                                   = nullptr;
    BoolVector     *p_primary_mc_IsProton                                  = nullptr;
    BoolVector     *p_primary_mc_IsPionOrMuon                              = nullptr;
    BoolVector     *p_primary_mc_IsCosmicRay                               = nullptr;
    IntVector      *p_primary_mc_PdgCode                                   = nullptr;
    FloatVector    *p_primary_mc_HitPurity                                 = nullptr;
    FloatVector    *p_primary_mc_HitCompleteness                           = nullptr;
    FloatVector    *p_primary_mc_CollectionPlaneHitPurity                  = nullptr;
    FloatVector    *p_primary_mc_CollectionPlaneHitCompleteness            = nullptr;
    UInt_t          cr_Number                                              = 0U;
    BoolVector     *p_cr_WasReconstructed                                  = nullptr;
    BoolVector     *p_cr_IsVertexFiducial                                  = nullptr;
    BoolVector     *p_cr_IsContained                                       = nullptr;
    FloatVector    *p_cr_FiducialHitFraction                               = nullptr;
    BoolVector     *p_cr_HasMcInfo                                         = nullptr;
    FloatVector    *p_cr_KineticEnergy                                     = nullptr;
    FloatVector    *p_cr_KineticEnergyFracFromRange                        = nullptr;
    FloatVector    *p_cr_KineticEnergyFracFromCorrectedTrackCharge         = nullptr;
    FloatVector    *p_cr_KineticEnergyFracFromUncorrectedTrackCharge       = nullptr;
    FloatVector    *p_cr_KineticEnergyFracFromShowerCharge                 = nullptr;
    FloatVector    *p_cr_VertexX                                           = nullptr;
    FloatVector    *p_cr_VertexY                                           = nullptr;
    FloatVector    *p_cr_VertexZ                                           = nullptr;
    FloatVector    *p_cr_DirectionCosineX                                  = nullptr;
    FloatVector    *p_cr_DirectionCosineY                                  = nullptr;
    FloatVector    *p_cr_DirectionCosineZ                                  = nullptr;
    StringVector   *p_cr_TypeTree                                          = nullptr;
    UnsignedVector *p_cr_NumberOf3dHits                                    = nullptr;
    UnsignedVector *p_cr_NumberOfCollectionPlaneHits                       = nullptr;
    UnsignedVector *p_cr_NumberOfDownstreamParticles                       = nullptr;
    UInt64Vector   *p_cr_mc_McParticleUid                                  = nullptr;
    BoolVector     *p_cr_mc_IsParticleSplitByReco                          = nullptr;
    FloatVector    *p_cr_mc_Energy                                         = nullptr;
    FloatVector    *p_cr_mc_KineticEnergy                                  = nullptr;
    FloatVector    *p_cr_mc_Mass                                           = nullptr;
    FloatVector    *p_cr_mc_VertexX                                        = nullptr;
    FloatVector    *p_cr_mc_VertexY                                        = nullptr;
    FloatVector    *p_cr_mc_VertexZ                                        = nullptr;
    FloatVector    *p_cr_mc_DirectionCosineX                               = nullptr;
    FloatVector    *p_cr_mc_DirectionCosineY                               = nullptr;
    FloatVector    *p_cr_mc_DirectionCosineZ                               = nullptr;
    FloatVector    *p_cr_mc_Momentum                                       = nullptr;
    FloatVector    *p_cr_mc_MomentumX                                      = nullptr;
    FloatVector    *p_cr_mc_MomentumY                                      = nullptr;
    FloatVector    *p_cr_mc_MomentumZ                                      = nullptr;
    BoolVector     *p_cr_mc_IsVertexFiducial                               = nullptr;
    BoolVector     *p_cr_mc_IsContained                                    = nullptr;
    FloatVector    *p_cr_mc_ContainmentFraction                            = nullptr;
    StringVector   *p_cr_mc_TypeTree                                       = nullptr;
    BoolVector     *p_cr_mc_IsShower                                       = nullptr;
    BoolVector     *p_cr_mc_IsTrack                                        = nullptr;
    BoolVector     *p_cr_mc_IsProton                                       = nullptr;
    BoolVector     *p_cr_mc_IsPionOrMuon                                   = nullptr;
    BoolVector     *p_cr_mc_IsCosmicRay                                    = nullptr;
    IntVector      *p_cr_mc_PdgCode                                        = nullptr;
    FloatVector    *p_cr_mc_HitPurity                                      = nullptr;
    FloatVector    *p_cr_mc_HitCompleteness                                = nullptr;
    FloatVector    *p_cr_mc_CollectionPlaneHitPurity                       = nullptr;
    FloatVector    *p_cr_mc_CollectionPlaneHitCompleteness                 = nullptr;
    
    pTree->SetBranchAddress("nu_WasReconstructed",                                 &nu_WasReconstructed);
    pTree->SetBranchAddress("nu_IsVertexFiducial",                                 &nu_IsVertexFiducial);
    pTree->SetBranchAddress("nu_IsContained",                                      &nu_IsContained);
    pTree->SetBranchAddress("nu_FiducialHitFraction",                              &nu_FiducialHitFraction);
    pTree->SetBranchAddress("nu_HasMcInfo",                                        &nu_HasMcInfo);
    pTree->SetBranchAddress("nu_VisibleEnergy",                                    &nu_VisibleEnergy);
    pTree->SetBranchAddress("nu_VisibleLongitudinalEnergy",                        &nu_VisibleLongitudinalEnergy);
    pTree->SetBranchAddress("nu_VisibleTransverseEnergy",                          &nu_VisibleTransverseEnergy);
    pTree->SetBranchAddress("nu_VisibleEnergyFracFromRange",                       &nu_VisibleEnergyFracFromRange);
    pTree->SetBranchAddress("nu_VisibleEnergyFracFromCorrectedTrackCharge",        &nu_VisibleEnergyFracFromCorrectedTrackCharge);
    pTree->SetBranchAddress("nu_VisibleEnergyFracFromUncorrectedTrackCharge",      &nu_VisibleEnergyFracFromUncorrectedTrackCharge);
    pTree->SetBranchAddress("nu_VisibleEnergyFracFromShowerCharge",                &nu_VisibleEnergyFracFromShowerCharge);
    pTree->SetBranchAddress("nu_VertexX",                                          &nu_VertexX);
    pTree->SetBranchAddress("nu_VertexY",                                          &nu_VertexY);
    pTree->SetBranchAddress("nu_VertexZ",                                          &nu_VertexZ);
    pTree->SetBranchAddress("nu_DirectionCosineX",                                 &nu_DirectionCosineX);
    pTree->SetBranchAddress("nu_DirectionCosineY",                                 &nu_DirectionCosineY);
    pTree->SetBranchAddress("nu_DirectionCosineZ",                                 &nu_DirectionCosineZ);
    pTree->SetBranchAddress("nu_TypeTree",                                         &p_nu_TypeTree);
    pTree->SetBranchAddress("nu_NumberOf3dHits",                                   &nu_NumberOf3dHits);
    pTree->SetBranchAddress("nu_NumberOfCollectionPlaneHits",                      &nu_NumberOfCollectionPlaneHits);
    pTree->SetBranchAddress("nu_NumberOfDownstreamParticles",                      &nu_NumberOfDownstreamParticles);
    pTree->SetBranchAddress("nu_mc_McParticleUid",                                 &nu_mc_McParticleUid);
    pTree->SetBranchAddress("nu_mc_Energy",                                        &nu_mc_Energy);
    pTree->SetBranchAddress("nu_mc_LongitudinalEnergy",                            &nu_mc_LongitudinalEnergy);
    pTree->SetBranchAddress("nu_mc_TransverseEnergy",                              &nu_mc_TransverseEnergy);
    pTree->SetBranchAddress("nu_mc_VisibleEnergy",                                 &nu_mc_VisibleEnergy);
    pTree->SetBranchAddress("nu_mc_VisibleLongitudinalEnergy",                     &nu_mc_VisibleLongitudinalEnergy);
    pTree->SetBranchAddress("nu_mc_VisibleTransverseEnergy",                       &nu_mc_VisibleTransverseEnergy);
    pTree->SetBranchAddress("nu_mc_VertexX",                                       &nu_mc_VertexX);
    pTree->SetBranchAddress("nu_mc_VertexY",                                       &nu_mc_VertexY);
    pTree->SetBranchAddress("nu_mc_VertexZ",                                       &nu_mc_VertexZ);
    pTree->SetBranchAddress("nu_mc_DirectionCosineX",                              &nu_mc_DirectionCosineX);
    pTree->SetBranchAddress("nu_mc_DirectionCosineY",                              &nu_mc_DirectionCosineY);
    pTree->SetBranchAddress("nu_mc_DirectionCosineZ",                              &nu_mc_DirectionCosineZ);
    pTree->SetBranchAddress("nu_mc_Momentum",                                      &nu_mc_Momentum);
    pTree->SetBranchAddress("nu_mc_MomentumX",                                     &nu_mc_MomentumX);
    pTree->SetBranchAddress("nu_mc_MomentumY",                                     &nu_mc_MomentumY);
    pTree->SetBranchAddress("nu_mc_MomentumZ",                                     &nu_mc_MomentumZ);
    pTree->SetBranchAddress("nu_mc_IsVertexFiducial",                              &nu_mc_IsVertexFiducial);
    pTree->SetBranchAddress("nu_mc_IsContained",                                   &nu_mc_IsContained);
    pTree->SetBranchAddress("nu_mc_ContainmentFraction",                           &nu_mc_ContainmentFraction);
    pTree->SetBranchAddress("nu_mc_TypeTree",                                      &p_nu_mc_TypeTree);
    pTree->SetBranchAddress("nu_mc_InteractionType",                               &p_nu_mc_InteractionType);
    pTree->SetBranchAddress("nu_mc_IsChargedCurrent",                              &nu_mc_IsChargedCurrent);
    pTree->SetBranchAddress("nu_mc_VisibleEnergyFraction",                         &nu_mc_VisibleEnergyFraction);
    pTree->SetBranchAddress("nu_mc_PdgCode",                                       &nu_mc_PdgCode);
    pTree->SetBranchAddress("nu_mc_HitPurity",                                     &nu_mc_HitPurity);
    pTree->SetBranchAddress("nu_mc_HitCompleteness",                               &nu_mc_HitCompleteness);
    pTree->SetBranchAddress("nu_mc_CollectionPlaneHitPurity",                      &nu_mc_CollectionPlaneHitPurity);
    pTree->SetBranchAddress("nu_mc_CollectionPlaneHitCompleteness",                &nu_mc_CollectionPlaneHitCompleteness);
    pTree->SetBranchAddress("primary_Number",                                      &primary_Number);
    pTree->SetBranchAddress("primary_WasReconstructed",                            &p_primary_WasReconstructed);
    pTree->SetBranchAddress("primary_IsVertexFiducial",                            &p_primary_IsVertexFiducial);
    pTree->SetBranchAddress("primary_IsContained",                                 &p_primary_IsContained);
    pTree->SetBranchAddress("primary_FiducialHitFraction",                         &p_primary_FiducialHitFraction);
    pTree->SetBranchAddress("primary_HasMcInfo",                                   &p_primary_HasMcInfo);
    pTree->SetBranchAddress("primary_KineticEnergy",                               &p_primary_KineticEnergy);
    pTree->SetBranchAddress("primary_KineticEnergyFracFromRange",                  &p_primary_KineticEnergyFracFromRange);
    pTree->SetBranchAddress("primary_KineticEnergyFracFromCorrectedTrackCharge",   &p_primary_KineticEnergyFracFromCorrectedTrackCharge);
    pTree->SetBranchAddress("primary_KineticEnergyFracFromUncorrectedTrackCharge", &p_primary_KineticEnergyFracFromUncorrectedTrackCharge);
    pTree->SetBranchAddress("primary_KineticEnergyFracFromShowerCharge",           &p_primary_KineticEnergyFracFromShowerCharge);
    pTree->SetBranchAddress("primary_VertexX",                                     &p_primary_VertexX);
    pTree->SetBranchAddress("primary_VertexY",                                     &p_primary_VertexY);
    pTree->SetBranchAddress("primary_VertexZ",                                     &p_primary_VertexZ);
    pTree->SetBranchAddress("primary_DirectionCosineX",                            &p_primary_DirectionCosineX);
    pTree->SetBranchAddress("primary_DirectionCosineY",                            &p_primary_DirectionCosineY);
    pTree->SetBranchAddress("primary_DirectionCosineZ",                            &p_primary_DirectionCosineZ);
    pTree->SetBranchAddress("primary_TypeTree",                                    &p_primary_TypeTree);
    pTree->SetBranchAddress("primary_IsShower",                                    &p_primary_IsShower);
    pTree->SetBranchAddress("primary_IsTrack",                                     &p_primary_IsTrack);
    pTree->SetBranchAddress("primary_IsProton",                                    &p_primary_IsProton);
    pTree->SetBranchAddress("primary_IsPionOrMuon",                                &p_primary_IsPionOrMuon);
    pTree->SetBranchAddress("primary_NumberOf3dHits",                              &p_primary_NumberOf3dHits);
    pTree->SetBranchAddress("primary_NumberOfCollectionPlaneHits",                 &p_primary_NumberOfCollectionPlaneHits);
    pTree->SetBranchAddress("primary_NumberOfDownstreamParticles",                 &p_primary_NumberOfDownstreamParticles);
    pTree->SetBranchAddress("primary_mc_McParticleUid",                            &p_primary_mc_McParticleUid);
    pTree->SetBranchAddress("primary_mc_IsParticleSplitByReco",                    &p_primary_mc_IsParticleSplitByReco);
    pTree->SetBranchAddress("primary_mc_Energy",                                   &p_primary_mc_Energy);
    pTree->SetBranchAddress("primary_mc_KineticEnergy",                            &p_primary_mc_KineticEnergy);
    pTree->SetBranchAddress("primary_mc_Mass",                                     &p_primary_mc_Mass);
    pTree->SetBranchAddress("primary_mc_VertexX",                                  &p_primary_mc_VertexX);
    pTree->SetBranchAddress("primary_mc_VertexY",                                  &p_primary_mc_VertexY);
    pTree->SetBranchAddress("primary_mc_VertexZ",                                  &p_primary_mc_VertexZ);
    pTree->SetBranchAddress("primary_mc_DirectionCosineX",                         &p_primary_mc_DirectionCosineX);
    pTree->SetBranchAddress("primary_mc_DirectionCosineY",                         &p_primary_mc_DirectionCosineY);
    pTree->SetBranchAddress("primary_mc_DirectionCosineZ",                         &p_primary_mc_DirectionCosineZ);
    pTree->SetBranchAddress("primary_mc_Momentum",                                 &p_primary_mc_Momentum);
    pTree->SetBranchAddress("primary_mc_MomentumX",                                &p_primary_mc_MomentumX);
    pTree->SetBranchAddress("primary_mc_MomentumY",                                &p_primary_mc_MomentumY);
    pTree->SetBranchAddress("primary_mc_MomentumZ",                                &p_primary_mc_MomentumZ);
    pTree->SetBranchAddress("primary_mc_IsVertexFiducial",                         &p_primary_mc_IsVertexFiducial);
    pTree->SetBranchAddress("primary_mc_IsContained",                              &p_primary_mc_IsContained);
    pTree->SetBranchAddress("primary_mc_ContainmentFraction",                      &p_primary_mc_ContainmentFraction);
    pTree->SetBranchAddress("primary_mc_TypeTree",                                 &p_primary_mc_TypeTree);
    pTree->SetBranchAddress("primary_mc_IsShower",                                 &p_primary_mc_IsShower);
    pTree->SetBranchAddress("primary_mc_IsTrack",                                  &p_primary_mc_IsTrack);
    pTree->SetBranchAddress("primary_mc_IsProton",                                 &p_primary_mc_IsProton);
    pTree->SetBranchAddress("primary_mc_IsPionOrMuon",                             &p_primary_mc_IsPionOrMuon);
    pTree->SetBranchAddress("primary_mc_IsCosmicRay",                              &p_primary_mc_IsCosmicRay);
    pTree->SetBranchAddress("primary_mc_PdgCode",                                  &p_primary_mc_PdgCode);
    pTree->SetBranchAddress("primary_mc_HitPurity",                                &p_primary_mc_HitPurity);
    pTree->SetBranchAddress("primary_mc_HitCompleteness",                          &p_primary_mc_HitCompleteness);
    pTree->SetBranchAddress("primary_mc_CollectionPlaneHitPurity",                 &p_primary_mc_CollectionPlaneHitPurity);
    pTree->SetBranchAddress("primary_mc_CollectionPlaneHitCompleteness",           &p_primary_mc_CollectionPlaneHitCompleteness);
    pTree->SetBranchAddress("cr_Number",                                           &cr_Number);
    pTree->SetBranchAddress("cr_WasReconstructed",                                 &p_cr_WasReconstructed);
    pTree->SetBranchAddress("cr_IsVertexFiducial",                                 &p_cr_IsVertexFiducial);
    pTree->SetBranchAddress("cr_IsContained",                                      &p_cr_IsContained);
    pTree->SetBranchAddress("cr_FiducialHitFraction",                              &p_cr_FiducialHitFraction);
    pTree->SetBranchAddress("cr_HasMcInfo",                                        &p_cr_HasMcInfo);
    pTree->SetBranchAddress("cr_KineticEnergy",                                    &p_cr_KineticEnergy);
    pTree->SetBranchAddress("cr_KineticEnergyFracFromRange",                       &p_cr_KineticEnergyFracFromRange);
    pTree->SetBranchAddress("cr_KineticEnergyFracFromCorrectedTrackCharge",        &p_cr_KineticEnergyFracFromCorrectedTrackCharge);
    pTree->SetBranchAddress("cr_KineticEnergyFracFromUncorrectedTrackCharge",      &p_cr_KineticEnergyFracFromUncorrectedTrackCharge);
    pTree->SetBranchAddress("cr_KineticEnergyFracFromShowerCharge",                &p_cr_KineticEnergyFracFromShowerCharge);
    pTree->SetBranchAddress("cr_VertexX",                                          &p_cr_VertexX);
    pTree->SetBranchAddress("cr_VertexY",                                          &p_cr_VertexY);
    pTree->SetBranchAddress("cr_VertexZ",                                          &p_cr_VertexZ);
    pTree->SetBranchAddress("cr_DirectionCosineX",                                 &p_cr_DirectionCosineX);
    pTree->SetBranchAddress("cr_DirectionCosineY",                                 &p_cr_DirectionCosineY);
    pTree->SetBranchAddress("cr_DirectionCosineZ",                                 &p_cr_DirectionCosineZ);
    pTree->SetBranchAddress("cr_TypeTree",                                         &p_cr_TypeTree);
    pTree->SetBranchAddress("cr_NumberOf3dHits",                                   &p_cr_NumberOf3dHits);
    pTree->SetBranchAddress("cr_NumberOfCollectionPlaneHits",                      &p_cr_NumberOfCollectionPlaneHits);
    pTree->SetBranchAddress("cr_NumberOfDownstreamParticles",                      &p_cr_NumberOfDownstreamParticles);
    pTree->SetBranchAddress("cr_mc_McParticleUid",                                 &p_cr_mc_McParticleUid);
    pTree->SetBranchAddress("cr_mc_IsParticleSplitByReco",                         &p_cr_mc_IsParticleSplitByReco);
    pTree->SetBranchAddress("cr_mc_Energy",                                        &p_cr_mc_Energy);
    pTree->SetBranchAddress("cr_mc_KineticEnergy",                                 &p_cr_mc_KineticEnergy);
    pTree->SetBranchAddress("cr_mc_Mass",                                          &p_cr_mc_Mass);
    pTree->SetBranchAddress("cr_mc_VertexX",                                       &p_cr_mc_VertexX);
    pTree->SetBranchAddress("cr_mc_VertexY",                                       &p_cr_mc_VertexY);
    pTree->SetBranchAddress("cr_mc_VertexZ",                                       &p_cr_mc_VertexZ);
    pTree->SetBranchAddress("cr_mc_DirectionCosineX",                              &p_cr_mc_DirectionCosineX);
    pTree->SetBranchAddress("cr_mc_DirectionCosineY",                              &p_cr_mc_DirectionCosineY);
    pTree->SetBranchAddress("cr_mc_DirectionCosineZ",                              &p_cr_mc_DirectionCosineZ);
    pTree->SetBranchAddress("cr_mc_Momentum",                                      &p_cr_mc_Momentum);
    pTree->SetBranchAddress("cr_mc_MomentumX",                                     &p_cr_mc_MomentumX);
    pTree->SetBranchAddress("cr_mc_MomentumY",                                     &p_cr_mc_MomentumY);
    pTree->SetBranchAddress("cr_mc_MomentumZ",                                     &p_cr_mc_MomentumZ);
    pTree->SetBranchAddress("cr_mc_IsVertexFiducial",                              &p_cr_mc_IsVertexFiducial);
    pTree->SetBranchAddress("cr_mc_IsContained",                                   &p_cr_mc_IsContained);
    pTree->SetBranchAddress("cr_mc_ContainmentFraction",                           &p_cr_mc_ContainmentFraction);
    pTree->SetBranchAddress("cr_mc_TypeTree",                                      &p_cr_mc_TypeTree);
    pTree->SetBranchAddress("cr_mc_IsShower",                                      &p_cr_mc_IsShower);
    pTree->SetBranchAddress("cr_mc_IsTrack",                                       &p_cr_mc_IsTrack);
    pTree->SetBranchAddress("cr_mc_IsProton",                                      &p_cr_mc_IsProton);
    pTree->SetBranchAddress("cr_mc_IsPionOrMuon",                                  &p_cr_mc_IsPionOrMuon);
    pTree->SetBranchAddress("cr_mc_IsCosmicRay",                                   &p_cr_mc_IsCosmicRay);
    pTree->SetBranchAddress("cr_mc_PdgCode",                                       &p_cr_mc_PdgCode);
    pTree->SetBranchAddress("cr_mc_HitPurity",                                     &p_cr_mc_HitPurity);
    pTree->SetBranchAddress("cr_mc_HitCompleteness",                               &p_cr_mc_HitCompleteness);
    pTree->SetBranchAddress("cr_mc_CollectionPlaneHitPurity",                      &p_cr_mc_CollectionPlaneHitPurity);
    pTree->SetBranchAddress("cr_mc_CollectionPlaneHitCompleteness",                &p_cr_mc_CollectionPlaneHitCompleteness);
    
    
    // Define some new ntuples to plot later.
    TFile *const pFile = new TFile("tmp.root", "RECREATE");
    
    
    TNtuple *const pNtuple_Neutrino_0 = new TNtuple("Neutrino_0", "Neutrino_0",
        "nu_VisibleEnergy:"
        "nu_VisibleLongitudinalEnergy:"
        "nu_VisibleTransverseEnergy:"
        "nu_VisibleEnergyFracFromRange:"
        "nu_VisibleEnergyFracFromCorrectedTrackCharge:"
        "nu_VisibleEnergyFracFromUncorrectedTrackCharge:"
        "nu_VisibleEnergyFracFromShowerCharge:"
        "nu_mc_Energy:"
        "nu_mc_LongitudinalEnergy:"
        "nu_mc_TransverseEnergy:"
        "nu_mc_VisibleEnergy:"
        "nu_mc_VisibleLongitudinalEnergy:"
        "nu_mc_VisibleTransverseEnergy");
        
    TNtuple *const pNtuple_Neutrino_1 = new TNtuple("Neutrino_1", "Neutrino_1",
        "nu_DirectionCosineX:"
        "nu_DirectionCosineY:"
        "nu_DirectionCosineZ:"
        "nu_mc_DirectionCosineX:"
        "nu_mc_DirectionCosineY:"
        "nu_mc_DirectionCosineZ:"
        "nu_mc_Momentum:"
        "nu_mc_MomentumX:"
        "nu_mc_MomentumY:"
        "nu_mc_MomentumZ");
        
    TNtuple *const pNtuple_Neutrino_2 = new TNtuple("Neutrino_2", "Neutrino_2",
        "nu_WasReconstructed:"
        "nu_IsVertexFiducial:"
        "nu_IsContained:"
        "nu_FiducialHitFraction:"
        "nu_mc_IsVertexFiducial:"
        "nu_mc_IsContained:"
        "nu_mc_ContainmentFraction:"
        "nu_mc_McParticleUid:"
        "nu_mc_IsChargedCurrent:"
        "nu_mc_VisibleEnergyFraction:"
        "nu_mc_PdgCode:"
        "nu_mc_HitPurity:"
        "nu_mc_HitCompleteness:"
        "nu_mc_CollectionPlaneHitPurity:"
        "nu_mc_CollectionPlaneHitCompleteness");
        
    TNtuple *const pNtuple_Neutrino_3 = new TNtuple("Neutrino_3", "Neutrino_3",
        "nu_HasMcInfo:"
        "primary_Number:"
        "cr_Number:"
        "nu_NumberOf3dHits:"
        "nu_NumberOfCollectionPlaneHits:"
        "nu_NumberOfDownstreamParticles:"
        "nu_VertexX:"
        "nu_VertexY:"
        "nu_VertexZ:"
        "nu_mc_VertexX:"
        "nu_mc_VertexY:"
        "nu_mc_VertexZ");

    TNtuple *const pNtuple_Primary_0 = new TNtuple("Primary_0", "Primary_0",
        "primary_KineticEnergy:"
        "primary_KineticEnergyFracFromRange:"
        "primary_KineticEnergyFracFromCorrectedTrackCharge:"
        "primary_KineticEnergyFracFromUncorrectedTrackCharge:"
        "primary_KineticEnergyFracFromShowerCharge:"
        "primary_VertexX:"
        "primary_VertexY:"
        "primary_VertexZ:"
        "primary_mc_Energy:"
        "primary_mc_KineticEnergy:"
        "primary_mc_Mass:"
        "primary_mc_VertexX:"
        "primary_mc_VertexY:"
        "primary_mc_VertexZ");
        
    TNtuple *const pNtuple_Primary_1 = new TNtuple("Primary_1", "Primary_1",
        "primary_DirectionCosineX:"
        "primary_DirectionCosineY:"
        "primary_DirectionCosineZ:"
        "primary_mc_DirectionCosineX:"
        "primary_mc_DirectionCosineY:"
        "primary_mc_DirectionCosineZ:"
        "primary_mc_Momentum:"
        "primary_mc_MomentumX:"
        "primary_mc_MomentumY:"
        "primary_mc_MomentumZ");
         
    TNtuple *const pNtuple_Primary_2 = new TNtuple("Primary_2", "Primary_2",
        "primary_WasReconstructed:"
        "primary_IsVertexFiducial:"
        "primary_IsContained:"
        "primary_FiducialHitFraction:"
        "primary_HasMcInfo:"
        "primary_mc_IsVertexFiducial:"
        "primary_mc_IsContained:"
        "primary_mc_ContainmentFraction");
        
    TNtuple *const pNtuple_Primary_3 = new TNtuple("Primary_3", "Primary_3",
        "primary_IsShower:"
        "primary_IsTrack:"
        "primary_IsProton:"
        "primary_IsPionOrMuon:"
        "primary_mc_IsShower:"
        "primary_mc_IsTrack:"
        "primary_mc_IsProton:"
        "primary_mc_IsPionOrMuon:"
        "primary_mc_IsCosmicRay:"
        "primary_mc_PdgCode");
        
    TNtuple *const pNtuple_Primary_4 = new TNtuple("Primary_4", "Primary_4",
        "primary_NumberOf3dHits:"
        "primary_NumberOfCollectionPlaneHits:"
        "primary_NumberOfDownstreamParticles:"
        "primary_mc_McParticleUid:"
        "primary_mc_IsParticleSplitByReco:"
        "primary_mc_HitPurity:"
        "primary_mc_HitCompleteness:"
        "primary_mc_CollectionPlaneHitPurity:"
        "primary_mc_CollectionPlaneHitCompleteness");
        
    TNtuple *const pNtuple_CosmicRay_0 = new TNtuple("Cosmic_Ray_0", "Cosmic_Ray_0",
        "cr_KineticEnergy:"
        "cr_KineticEnergyFracFromRange:"
        "cr_KineticEnergyFracFromCorrectedTrackCharge:"
        "cr_KineticEnergyFracFromUncorrectedTrackCharge:"
        "cr_KineticEnergyFracFromShowerCharge:"
        "cr_VertexX:"
        "cr_VertexY:"
        "cr_VertexZ:"
        "cr_mc_Energy:"
        "cr_mc_KineticEnergy:"
        "cr_mc_Mass:"
        "cr_mc_VertexX:"
        "cr_mc_VertexY:"
        "cr_mc_VertexZ");
        
    TNtuple *const pNtuple_CosmicRay_1 = new TNtuple("Cosmic_Ray_1", "Cosmic_Ray_1",
        "cr_DirectionCosineX:"
        "cr_DirectionCosineY:"
        "cr_DirectionCosineZ:"
        "cr_mc_Momentum:"
        "cr_mc_MomentumX:"
        "cr_mc_MomentumY:"
        "cr_mc_MomentumZ:"
        "cr_mc_DirectionCosineX:"
        "cr_mc_DirectionCosineY:"
        "cr_mc_DirectionCosineZ");
        
    TNtuple *const pNtuple_CosmicRay_2 = new TNtuple("Cosmic_Ray_2", "Cosmic_Ray_2",
        "cr_WasReconstructed:"
        "cr_IsVertexFiducial:"
        "cr_IsContained:"
        "cr_FiducialHitFraction:"
        "cr_HasMcInfo:"
        "cr_mc_IsVertexFiducial:"
        "cr_mc_IsContained:"
        "cr_mc_ContainmentFraction:"
        "cr_mc_IsShower:"
        "cr_mc_IsTrack:"
        "cr_mc_IsProton:"
        "cr_mc_IsPionOrMuon:"
        "cr_mc_IsCosmicRay:"
        "cr_mc_PdgCode");
        
    TNtuple *const pNtuple_CosmicRay_3 = new TNtuple("Cosmic_Ray_3", "Cosmic_Ray_3",
        "cr_NumberOf3dHits:"
        "cr_NumberOfCollectionPlaneHits:"
        "cr_NumberOfDownstreamParticles:"
        "cr_mc_McParticleUid:"
        "cr_mc_IsParticleSplitByReco:"
        "cr_mc_HitPurity:"
        "cr_mc_HitCompleteness:"
        "cr_mc_CollectionPlaneHitPurity:"
        "cr_mc_CollectionPlaneHitCompleteness");
    
    // Loop over the input primaries.
    const int nEntries = pTree->GetEntries();
    
    unsigned type_RecoProton_McProton(0U);
    unsigned type_RecoProton_McPionMuon(0U);
    unsigned type_RecoProton_McOtherTrack(0U);
    unsigned type_RecoProton_McShower(0U);
    
    unsigned type_RecoPionMuon_McProton(0U);
    unsigned type_RecoPionMuon_McPionMuon(0U);
    unsigned type_RecoPionMuon_McOtherTrack(0U);
    unsigned type_RecoPionMuon_McShower(0U);
    
    unsigned type_RecoOtherTrack_McProton(0U);
    unsigned type_RecoOtherTrack_McPionMuon(0U);
    unsigned type_RecoOtherTrack_McOtherTrack(0U);
    unsigned type_RecoOtherTrack_McShower(0U);
    
    unsigned type_RecoShower_McProton(0U);
    unsigned type_RecoShower_McPionMuon(0U);
    unsigned type_RecoShower_McOtherTrack(0U);
    unsigned type_RecoShower_McShower(0U);
    
    for (int i = 0; i < nEntries; ++i)
    {
        pTree->GetEntry(i);
        
        if (nu_WasReconstructed && nu_HasMcInfo && nu_mc_IsChargedCurrent && nu_mc_IsContained && nu_mc_IsVertexFiducial)
        {
            pNtuple_Neutrino_0->Fill(
                nu_VisibleEnergy, 
                nu_VisibleLongitudinalEnergy, 
                nu_VisibleTransverseEnergy, 
                nu_VisibleEnergyFracFromRange,
                nu_VisibleEnergyFracFromCorrectedTrackCharge,
                nu_VisibleEnergyFracFromUncorrectedTrackCharge,
                nu_VisibleEnergyFracFromShowerCharge,
                nu_mc_Energy,
                nu_mc_LongitudinalEnergy,
                nu_mc_TransverseEnergy,
                nu_mc_VisibleEnergy,
                nu_mc_VisibleLongitudinalEnergy,
                nu_mc_VisibleTransverseEnergy);
        
            pNtuple_Neutrino_1->Fill(
                nu_DirectionCosineX,
                nu_DirectionCosineY,
                nu_DirectionCosineZ,
                nu_mc_DirectionCosineX,
                nu_mc_DirectionCosineY,
                nu_mc_DirectionCosineZ,
                nu_mc_Momentum,
                nu_mc_MomentumX,
                nu_mc_MomentumY,
                nu_mc_MomentumZ);
                
            pNtuple_Neutrino_2->Fill(
                nu_WasReconstructed,
                nu_IsVertexFiducial,
                nu_IsContained,
                nu_FiducialHitFraction,
                nu_mc_IsVertexFiducial,
                nu_mc_IsContained,
                nu_mc_ContainmentFraction,
                nu_mc_McParticleUid,
                nu_mc_IsChargedCurrent,
                nu_mc_VisibleEnergyFraction,
                nu_mc_PdgCode,
                nu_mc_HitPurity,
                nu_mc_HitCompleteness,
                nu_mc_CollectionPlaneHitPurity,
                nu_mc_CollectionPlaneHitCompleteness);
                
            pNtuple_Neutrino_3->Fill(
                nu_HasMcInfo,
                primary_Number,
                cr_Number,
                nu_NumberOf3dHits,
                nu_NumberOfCollectionPlaneHits,
                nu_NumberOfDownstreamParticles,
                nu_VertexX,
                nu_VertexY,
                nu_VertexZ,
                nu_mc_VertexX,
                nu_mc_VertexY,
                nu_mc_VertexZ);

            for (int i = 0; i < primary_Number; ++i)
            {
                if ((*p_primary_WasReconstructed)[i] && (*p_primary_HasMcInfo)[i] && (*p_primary_mc_IsVertexFiducial)[i] && (*p_primary_mc_IsContained)[i])
                {
                    const bool recoIsProton     = (*p_primary_IsProton)[i];
                    const bool mcIsProton       = (*p_primary_mc_IsProton)[i];
                    const bool recoIsPionMuon   = (*p_primary_IsPionOrMuon)[i];
                    const bool mcIsPionMuon     = (*p_primary_mc_IsPionOrMuon)[i];
                    const bool recoIsOtherTrack = (*p_primary_IsTrack)[i];
                    const bool mcIsOtherTrack   = (*p_primary_mc_IsTrack)[i];
                    const bool recoIsShower     = (*p_primary_IsShower)[i];
                    const bool mcIsShower       = (*p_primary_mc_IsShower)[i];

                    if      (recoIsProton && mcIsProton)         ++type_RecoProton_McProton;
                    else if (recoIsProton && mcIsPionMuon)       ++type_RecoProton_McPionMuon;
                    else if (recoIsProton && mcIsOtherTrack)     ++type_RecoProton_McOtherTrack;
                    else if (recoIsProton && mcIsShower)         ++type_RecoProton_McShower;
                                                                 
                    else if (recoIsPionMuon && mcIsProton)       ++type_RecoPionMuon_McProton;
                    else if (recoIsPionMuon && mcIsPionMuon)     ++type_RecoPionMuon_McPionMuon;
                    else if (recoIsPionMuon && mcIsOtherTrack)   ++type_RecoPionMuon_McOtherTrack;
                    else if (recoIsPionMuon && mcIsShower)       ++type_RecoPionMuon_McShower;
                    
                    else if (recoIsOtherTrack && mcIsProton)     ++type_RecoOtherTrack_McProton;
                    else if (recoIsOtherTrack && mcIsPionMuon)   ++type_RecoOtherTrack_McPionMuon;
                    else if (recoIsOtherTrack && mcIsOtherTrack) ++type_RecoOtherTrack_McOtherTrack;
                    else if (recoIsOtherTrack && mcIsShower)     ++type_RecoOtherTrack_McShower;
                    
                    else if (recoIsShower && mcIsProton)         ++type_RecoShower_McProton;
                    else if (recoIsShower && mcIsPionMuon)       ++type_RecoShower_McPionMuon;
                    else if (recoIsShower && mcIsOtherTrack)     ++type_RecoShower_McOtherTrack;
                    else if (recoIsShower && mcIsShower)         ++type_RecoShower_McShower;
                    
                    else
                        std::cout << "Unknown particle type" << std::endl;
                    
                    pNtuple_Primary_0->Fill(
                        (*p_primary_KineticEnergy)[i],
                        (*p_primary_KineticEnergyFracFromRange)[i],
                        (*p_primary_KineticEnergyFracFromCorrectedTrackCharge)[i],
                        (*p_primary_KineticEnergyFracFromUncorrectedTrackCharge)[i],
                        (*p_primary_KineticEnergyFracFromShowerCharge)[i],
                        (*p_primary_VertexX)[i],
                        (*p_primary_VertexY)[i],
                        (*p_primary_VertexZ)[i],
                        (*p_primary_mc_Energy)[i],
                        (*p_primary_mc_KineticEnergy)[i],
                        (*p_primary_mc_Mass)[i],
                        (*p_primary_mc_VertexX)[i],
                        (*p_primary_mc_VertexY)[i],
                        (*p_primary_mc_VertexZ)[i]);
                        
                    pNtuple_Primary_1->Fill(
                        (*p_primary_DirectionCosineX)[i],
                        (*p_primary_DirectionCosineY)[i],
                        (*p_primary_DirectionCosineZ)[i],
                        (*p_primary_mc_DirectionCosineX)[i],
                        (*p_primary_mc_DirectionCosineY)[i],
                        (*p_primary_mc_DirectionCosineZ)[i],
                        (*p_primary_mc_Momentum)[i],
                        (*p_primary_mc_MomentumX)[i],
                        (*p_primary_mc_MomentumY)[i],
                        (*p_primary_mc_MomentumZ)[i]);
                         
                    pNtuple_Primary_2->Fill(
                        (*p_primary_WasReconstructed)[i],
                        (*p_primary_IsVertexFiducial)[i],
                        (*p_primary_IsContained)[i],
                        (*p_primary_FiducialHitFraction)[i],
                        (*p_primary_HasMcInfo)[i],
                        (*p_primary_mc_IsVertexFiducial)[i],
                        (*p_primary_mc_IsContained)[i],
                        (*p_primary_mc_ContainmentFraction)[i]);
                        
                    pNtuple_Primary_3->Fill(
                        (*p_primary_IsShower)[i],
                        (*p_primary_IsTrack)[i],
                        (*p_primary_IsProton)[i],
                        (*p_primary_IsPionOrMuon)[i],
                        (*p_primary_mc_IsShower)[i],
                        (*p_primary_mc_IsTrack)[i],
                        (*p_primary_mc_IsProton)[i],
                        (*p_primary_mc_IsPionOrMuon)[i],
                        (*p_primary_mc_IsCosmicRay)[i],
                        (*p_primary_mc_PdgCode)[i]);
                        
                    pNtuple_Primary_4->Fill(
                        (*p_primary_NumberOf3dHits)[i],
                        (*p_primary_NumberOfCollectionPlaneHits)[i],
                        (*p_primary_NumberOfDownstreamParticles)[i],
                        (*p_primary_mc_McParticleUid)[i],
                        (*p_primary_mc_IsParticleSplitByReco)[i],
                        (*p_primary_mc_HitPurity)[i],
                        (*p_primary_mc_HitCompleteness)[i],
                        (*p_primary_mc_CollectionPlaneHitPurity)[i],
                        (*p_primary_mc_CollectionPlaneHitCompleteness)[i]);
                }
            }
            
            for (int i = 0; i < cr_Number; ++i)
            {
                if ((*p_cr_WasReconstructed)[i] && (*p_cr_HasMcInfo)[i] && (*p_cr_mc_IsCosmicRay)[i] && 
                    (*p_cr_mc_HitCompleteness)[i] > 0.8f && (*p_cr_mc_HitPurity)[i] > 0.8f)
                {
                    pNtuple_CosmicRay_0->Fill(
                        (*p_cr_KineticEnergy)[i],
                        (*p_cr_KineticEnergyFracFromRange)[i],
                        (*p_cr_KineticEnergyFracFromCorrectedTrackCharge)[i],
                        (*p_cr_KineticEnergyFracFromUncorrectedTrackCharge)[i],
                        (*p_cr_KineticEnergyFracFromShowerCharge)[i],
                        (*p_cr_VertexX)[i],
                        (*p_cr_VertexY)[i],
                        (*p_cr_VertexZ)[i],
                        (*p_cr_mc_Energy)[i],
                        (*p_cr_mc_KineticEnergy)[i],
                        (*p_cr_mc_Mass)[i],
                        (*p_cr_mc_VertexX)[i],
                        (*p_cr_mc_VertexY)[i],
                        (*p_cr_mc_VertexZ)[i]);
                        
                    pNtuple_CosmicRay_1->Fill(
                        (*p_cr_DirectionCosineX)[i],
                        (*p_cr_DirectionCosineY)[i],
                        (*p_cr_DirectionCosineZ)[i],
                        (*p_cr_mc_Momentum)[i],
                        (*p_cr_mc_MomentumX)[i],
                        (*p_cr_mc_MomentumY)[i],
                        (*p_cr_mc_MomentumZ)[i],
                        (*p_cr_mc_DirectionCosineX)[i],
                        (*p_cr_mc_DirectionCosineY)[i],
                        (*p_cr_mc_DirectionCosineZ)[i]);
                        
                    pNtuple_CosmicRay_2->Fill(
                        (*p_cr_WasReconstructed)[i],
                        (*p_cr_IsVertexFiducial)[i],
                        (*p_cr_IsContained)[i],
                        (*p_cr_FiducialHitFraction)[i],
                        (*p_cr_HasMcInfo)[i],
                        (*p_cr_mc_IsVertexFiducial)[i],
                        (*p_cr_mc_IsContained)[i],
                        (*p_cr_mc_ContainmentFraction)[i],
                        (*p_cr_mc_IsShower)[i],
                        (*p_cr_mc_IsTrack)[i],
                        (*p_cr_mc_IsProton)[i],
                        (*p_cr_mc_IsPionOrMuon)[i],
                        (*p_cr_mc_IsCosmicRay)[i],
                        (*p_cr_mc_PdgCode)[i]);
                        
                    pNtuple_CosmicRay_3->Fill(
                        (*p_cr_NumberOf3dHits)[i],
                        (*p_cr_NumberOfCollectionPlaneHits)[i],
                        (*p_cr_NumberOfDownstreamParticles)[i],
                        (*p_cr_mc_McParticleUid)[i],
                        (*p_cr_mc_IsParticleSplitByReco)[i],
                        (*p_cr_mc_HitPurity)[i],
                        (*p_cr_mc_HitCompleteness)[i],
                        (*p_cr_mc_CollectionPlaneHitPurity)[i],
                        (*p_cr_mc_CollectionPlaneHitCompleteness)[i]);
                }
            }
        }
    }
    
    if (true)
    {
//        QUICK_PLOT_1D("nu_VisibleEnergy",                                    pNtuple_Neutrino_0);
//        QUICK_PLOT_1D("nu_VisibleLongitudinalEnergy",                        pNtuple_Neutrino_0);
//        QUICK_PLOT_1D("nu_VisibleTransverseEnergy",                          pNtuple_Neutrino_0);
//        QUICK_PLOT_1D("nu_VisibleEnergyFracFromRange",                       pNtuple_Neutrino_0);
//        QUICK_PLOT_1D("nu_VisibleEnergyFracFromCorrectedTrackCharge",        pNtuple_Neutrino_0);
//        QUICK_PLOT_1D("nu_VisibleEnergyFracFromUncorrectedTrackCharge",      pNtuple_Neutrino_0);
//        QUICK_PLOT_1D("nu_VisibleEnergyFracFromShowerCharge",                pNtuple_Neutrino_0);
//        QUICK_PLOT_1D("nu_mc_Energy",                                        pNtuple_Neutrino_0);
//        QUICK_PLOT_1D("nu_mc_LongitudinalEnergy",                            pNtuple_Neutrino_0);
//        QUICK_PLOT_1D("nu_mc_TransverseEnergy",                              pNtuple_Neutrino_0);
//        QUICK_PLOT_1D("nu_mc_VisibleEnergy",                                 pNtuple_Neutrino_0);
//        QUICK_PLOT_1D("nu_mc_VisibleLongitudinalEnergy",                     pNtuple_Neutrino_0);
//        QUICK_PLOT_1D("nu_mc_VisibleTransverseEnergy",                       pNtuple_Neutrino_0);
//                                                                               
//        QUICK_PLOT_1D("nu_DirectionCosineX",                                 pNtuple_Neutrino_1);
//        QUICK_PLOT_1D("nu_DirectionCosineY",                                 pNtuple_Neutrino_1);
//        QUICK_PLOT_1D("nu_DirectionCosineZ",                                 pNtuple_Neutrino_1);
//        QUICK_PLOT_1D("nu_mc_DirectionCosineX",                              pNtuple_Neutrino_1);
//        QUICK_PLOT_1D("nu_mc_DirectionCosineY",                              pNtuple_Neutrino_1);
//        QUICK_PLOT_1D("nu_mc_DirectionCosineZ",                              pNtuple_Neutrino_1);
//        QUICK_PLOT_1D("nu_mc_Momentum",                                      pNtuple_Neutrino_1);
//        QUICK_PLOT_1D("nu_mc_MomentumX",                                     pNtuple_Neutrino_1);
//        QUICK_PLOT_1D("nu_mc_MomentumY",                                     pNtuple_Neutrino_1);
//        QUICK_PLOT_1D("nu_mc_MomentumZ",                                     pNtuple_Neutrino_1);
//                                                                               
//        QUICK_PLOT_1D("nu_WasReconstructed",                                 pNtuple_Neutrino_2);
//        QUICK_PLOT_1D("nu_IsVertexFiducial",                                 pNtuple_Neutrino_2);
//        QUICK_PLOT_1D("nu_IsContained",                                      pNtuple_Neutrino_2);
//        QUICK_PLOT_1D("nu_FiducialHitFraction",                              pNtuple_Neutrino_2);
//        QUICK_PLOT_1D("nu_mc_IsVertexFiducial",                              pNtuple_Neutrino_2);
//        QUICK_PLOT_1D("nu_mc_IsContained",                                   pNtuple_Neutrino_2);
//        QUICK_PLOT_1D("nu_mc_ContainmentFraction",                           pNtuple_Neutrino_2);
//        QUICK_PLOT_1D("nu_mc_McParticleUid",                                 pNtuple_Neutrino_2);
//        QUICK_PLOT_1D("nu_mc_IsChargedCurrent",                              pNtuple_Neutrino_2);
//        QUICK_PLOT_1D("nu_mc_VisibleEnergyFraction",                         pNtuple_Neutrino_2);
//        QUICK_PLOT_1D("nu_mc_PdgCode",                                       pNtuple_Neutrino_2);
//        QUICK_PLOT_1D("nu_mc_HitPurity",                                     pNtuple_Neutrino_2);
//        QUICK_PLOT_1D("nu_mc_HitCompleteness",                               pNtuple_Neutrino_2);
//        QUICK_PLOT_1D("nu_mc_CollectionPlaneHitPurity",                      pNtuple_Neutrino_2);
//        QUICK_PLOT_1D("nu_mc_CollectionPlaneHitCompleteness",                pNtuple_Neutrino_2);
//                                                                               
//        QUICK_PLOT_1D("nu_HasMcInfo",                                        pNtuple_Neutrino_3);
//        QUICK_PLOT_1D("primary_Number",                                      pNtuple_Neutrino_3);
//        QUICK_PLOT_1D("cr_Number",                                           pNtuple_Neutrino_3);
//        QUICK_PLOT_1D("nu_NumberOf3dHits",                                   pNtuple_Neutrino_3);
//        QUICK_PLOT_1D("nu_NumberOfCollectionPlaneHits",                      pNtuple_Neutrino_3);
//        QUICK_PLOT_1D("nu_NumberOfDownstreamParticles",                      pNtuple_Neutrino_3);
//        QUICK_PLOT_1D("nu_VertexX",                                          pNtuple_Neutrino_3);
//        QUICK_PLOT_1D("nu_VertexY",                                          pNtuple_Neutrino_3);
//        QUICK_PLOT_1D("nu_VertexZ",                                          pNtuple_Neutrino_3);
//        QUICK_PLOT_1D("nu_mc_VertexX",                                       pNtuple_Neutrino_3);
//        QUICK_PLOT_1D("nu_mc_VertexY",                                       pNtuple_Neutrino_3);
//        QUICK_PLOT_1D("nu_mc_VertexZ",                                       pNtuple_Neutrino_3);
//        
//        QUICK_PLOT_1D("primary_KineticEnergy",                               pNtuple_Primary_0);
//        QUICK_PLOT_1D("primary_KineticEnergyFracFromRange",                  pNtuple_Primary_0);
//        QUICK_PLOT_1D("primary_KineticEnergyFracFromCorrectedTrackCharge",   pNtuple_Primary_0);
//        QUICK_PLOT_1D("primary_KineticEnergyFracFromUncorrectedTrackCharge", pNtuple_Primary_0);
//        QUICK_PLOT_1D("primary_KineticEnergyFracFromShowerCharge",           pNtuple_Primary_0);
//        QUICK_PLOT_1D("primary_VertexX",                                     pNtuple_Primary_0);
//        QUICK_PLOT_1D("primary_VertexY",                                     pNtuple_Primary_0);
//        QUICK_PLOT_1D("primary_VertexZ",                                     pNtuple_Primary_0);
//        QUICK_PLOT_1D("primary_mc_Energy",                                   pNtuple_Primary_0);
//        QUICK_PLOT_1D("primary_mc_KineticEnergy",                            pNtuple_Primary_0);
//        QUICK_PLOT_1D("primary_mc_Mass",                                     pNtuple_Primary_0);
//        QUICK_PLOT_1D("primary_mc_VertexX",                                  pNtuple_Primary_0);
//        QUICK_PLOT_1D("primary_mc_VertexY",                                  pNtuple_Primary_0);
//        QUICK_PLOT_1D("primary_mc_VertexZ",                                  pNtuple_Primary_0);
//        
//        QUICK_PLOT_1D("primary_DirectionCosineX",                            pNtuple_Primary_1);
//        QUICK_PLOT_1D("primary_DirectionCosineY",                            pNtuple_Primary_1);
//        QUICK_PLOT_1D("primary_DirectionCosineZ",                            pNtuple_Primary_1);
//        QUICK_PLOT_1D("primary_mc_DirectionCosineX",                         pNtuple_Primary_1);
//        QUICK_PLOT_1D("primary_mc_DirectionCosineY",                         pNtuple_Primary_1);
//        QUICK_PLOT_1D("primary_mc_DirectionCosineZ",                         pNtuple_Primary_1);
//        QUICK_PLOT_1D("primary_mc_Momentum",                                 pNtuple_Primary_1);
//        QUICK_PLOT_1D("primary_mc_MomentumX",                                pNtuple_Primary_1);
//        QUICK_PLOT_1D("primary_mc_MomentumY",                                pNtuple_Primary_1);
//        QUICK_PLOT_1D("primary_mc_MomentumZ",                                pNtuple_Primary_1);
//                                                                         
//        QUICK_PLOT_1D("primary_WasReconstructed",                            pNtuple_Primary_2);
//        QUICK_PLOT_1D("primary_IsVertexFiducial",                            pNtuple_Primary_2);
//        QUICK_PLOT_1D("primary_IsContained",                                 pNtuple_Primary_2);
//        QUICK_PLOT_1D("primary_FiducialHitFraction",                         pNtuple_Primary_2);
//        QUICK_PLOT_1D("primary_HasMcInfo",                                   pNtuple_Primary_2);
//        QUICK_PLOT_1D("primary_mc_IsVertexFiducial",                         pNtuple_Primary_2);
//        QUICK_PLOT_1D("primary_mc_IsContained",                              pNtuple_Primary_2);
//        QUICK_PLOT_1D("primary_mc_ContainmentFraction",                      pNtuple_Primary_2);
//                                                                             
//        QUICK_PLOT_1D("primary_IsShower",                                    pNtuple_Primary_3);
//        QUICK_PLOT_1D("primary_IsTrack",                                     pNtuple_Primary_3);
//        QUICK_PLOT_1D("primary_IsProton",                                    pNtuple_Primary_3);
//        QUICK_PLOT_1D("primary_IsPionOrMuon",                                pNtuple_Primary_3);
//        QUICK_PLOT_1D("primary_mc_IsShower",                                 pNtuple_Primary_3);
//        QUICK_PLOT_1D("primary_mc_IsTrack",                                  pNtuple_Primary_3);
//        QUICK_PLOT_1D("primary_mc_IsProton",                                 pNtuple_Primary_3);
//        QUICK_PLOT_1D("primary_mc_IsPionOrMuon",                             pNtuple_Primary_3);
//        QUICK_PLOT_1D("primary_mc_IsCosmicRay",                              pNtuple_Primary_3);
//        QUICK_PLOT_1D("primary_mc_PdgCode",                                  pNtuple_Primary_3);
//                                                                             
//        QUICK_PLOT_1D("primary_NumberOf3dHits",                              pNtuple_Primary_4);
//        QUICK_PLOT_1D("primary_NumberOfCollectionPlaneHits",                 pNtuple_Primary_4);
//        QUICK_PLOT_1D("primary_NumberOfDownstreamParticles",                 pNtuple_Primary_4);
//        QUICK_PLOT_1D("primary_mc_McParticleUid",                            pNtuple_Primary_4);
//        QUICK_PLOT_1D("primary_mc_IsParticleSplitByReco",                    pNtuple_Primary_4);
//        QUICK_PLOT_1D("primary_mc_HitPurity",                                pNtuple_Primary_4);
//        QUICK_PLOT_1D("primary_mc_HitCompleteness",                          pNtuple_Primary_4);
//        QUICK_PLOT_1D("primary_mc_CollectionPlaneHitPurity",                 pNtuple_Primary_4);
//        QUICK_PLOT_1D("primary_mc_CollectionPlaneHitCompleteness",           pNtuple_Primary_4);
//                                                                         
//        QUICK_PLOT_1D("cr_KineticEnergy",                                    pNtuple_CosmicRay_0);
//        QUICK_PLOT_1D("cr_KineticEnergyFracFromRange",                       pNtuple_CosmicRay_0);
//        QUICK_PLOT_1D("cr_KineticEnergyFracFromCorrectedTrackCharge",        pNtuple_CosmicRay_0);
//        QUICK_PLOT_1D("cr_KineticEnergyFracFromUncorrectedTrackCharge",      pNtuple_CosmicRay_0);
//        QUICK_PLOT_1D("cr_KineticEnergyFracFromShowerCharge",                pNtuple_CosmicRay_0);
//        QUICK_PLOT_1D("cr_VertexX",                                          pNtuple_CosmicRay_0);
//        QUICK_PLOT_1D("cr_VertexY",                                          pNtuple_CosmicRay_0);
//        QUICK_PLOT_1D("cr_VertexZ",                                          pNtuple_CosmicRay_0);
//        QUICK_PLOT_1D("cr_mc_Energy",                                        pNtuple_CosmicRay_0);
//        QUICK_PLOT_1D("cr_mc_KineticEnergy",                                 pNtuple_CosmicRay_0);
//        QUICK_PLOT_1D("cr_mc_Mass",                                          pNtuple_CosmicRay_0);
//        QUICK_PLOT_1D("cr_mc_VertexX",                                       pNtuple_CosmicRay_0);
//        QUICK_PLOT_1D("cr_mc_VertexY",                                       pNtuple_CosmicRay_0);
//        QUICK_PLOT_1D("cr_mc_VertexZ",                                       pNtuple_CosmicRay_0);
                                                                             
        QUICK_PLOT_1D("cr_DirectionCosineX",                                 pNtuple_CosmicRay_1);
        QUICK_PLOT_1D("cr_DirectionCosineY",                                 pNtuple_CosmicRay_1);
        QUICK_PLOT_1D("cr_DirectionCosineZ",                                 pNtuple_CosmicRay_1);
        QUICK_PLOT_1D("cr_mc_Momentum",                                      pNtuple_CosmicRay_1);
        QUICK_PLOT_1D("cr_mc_MomentumX",                                     pNtuple_CosmicRay_1);
        QUICK_PLOT_1D("cr_mc_MomentumY",                                     pNtuple_CosmicRay_1);
        QUICK_PLOT_1D("cr_mc_MomentumZ",                                     pNtuple_CosmicRay_1);
        QUICK_PLOT_1D("cr_mc_DirectionCosineX",                              pNtuple_CosmicRay_1);
        QUICK_PLOT_1D("cr_mc_DirectionCosineY",                              pNtuple_CosmicRay_1);
        QUICK_PLOT_1D("cr_mc_DirectionCosineZ",                              pNtuple_CosmicRay_1);
//                                                                             
//        QUICK_PLOT_1D("cr_WasReconstructed",                                 pNtuple_CosmicRay_2);
//        QUICK_PLOT_1D("cr_IsVertexFiducial",                                 pNtuple_CosmicRay_2);
//        QUICK_PLOT_1D("cr_IsContained",                                      pNtuple_CosmicRay_2);
//        QUICK_PLOT_1D("cr_FiducialHitFraction",                              pNtuple_CosmicRay_2);
//        QUICK_PLOT_1D("cr_HasMcInfo",                                        pNtuple_CosmicRay_2);
//        QUICK_PLOT_1D("cr_mc_IsVertexFiducial",                              pNtuple_CosmicRay_2);
//        QUICK_PLOT_1D("cr_mc_IsContained",                                   pNtuple_CosmicRay_2);
//        QUICK_PLOT_1D("cr_mc_ContainmentFraction",                           pNtuple_CosmicRay_2);
//        QUICK_PLOT_1D("cr_mc_IsShower",                                      pNtuple_CosmicRay_2);
//        QUICK_PLOT_1D("cr_mc_IsTrack",                                       pNtuple_CosmicRay_2);
//        QUICK_PLOT_1D("cr_mc_IsProton",                                      pNtuple_CosmicRay_2);
//        QUICK_PLOT_1D("cr_mc_IsPionOrMuon",                                  pNtuple_CosmicRay_2);
//        QUICK_PLOT_1D("cr_mc_IsCosmicRay",                                   pNtuple_CosmicRay_2);
        QUICK_PLOT_1D("cr_mc_PdgCode",                                       pNtuple_CosmicRay_2);
//                                                                             
//        QUICK_PLOT_1D("cr_NumberOf3dHits",                                   pNtuple_CosmicRay_3);
//        QUICK_PLOT_1D("cr_NumberOfCollectionPlaneHits",                      pNtuple_CosmicRay_3);
//        QUICK_PLOT_1D("cr_NumberOfDownstreamParticles",                      pNtuple_CosmicRay_3);
//        QUICK_PLOT_1D("cr_mc_McParticleUid",                                 pNtuple_CosmicRay_3);
//        QUICK_PLOT_1D("cr_mc_IsParticleSplitByReco",                         pNtuple_CosmicRay_3);
//        QUICK_PLOT_1D("cr_mc_HitPurity",                                     pNtuple_CosmicRay_3);
//        QUICK_PLOT_1D("cr_mc_HitCompleteness",                               pNtuple_CosmicRay_3);
//        QUICK_PLOT_1D("cr_mc_CollectionPlaneHitPurity",                      pNtuple_CosmicRay_3);
//        QUICK_PLOT_1D("cr_mc_CollectionPlaneHitCompleteness",                pNtuple_CosmicRay_3);
    }
    
    if (false)
    {
        PLOT_2D("nu_mc_VisibleEnergy",               "nu_VisibleEnergy",             pNtuple_Neutrino_0, 0.0, 2.0, 0.0, 2.0);
        PLOT_2D("nu_mc_VisibleLongitudinalEnergy",   "nu_VisibleTransverseEnergy", pNtuple_Neutrino_0, 0.0, 2.0, -0.5, 2.0);
        PLOT_2D("nu_mc_VisibleTransverseEnergy",     "nu_VisibleLongitudinalEnergy",   pNtuple_Neutrino_0, 0.0, 0.6, 0.0, 1.0);
        PLOT_2D("nu_mc_Energy",                      "nu_VisibleEnergy",             pNtuple_Neutrino_0, 0.0, 2.0, 0.0, 2.0);
        PLOT_2D("nu_mc_LongitudinalEnergy",          "nu_VisibleTransverseEnergy", pNtuple_Neutrino_0, 0.0, 2.0, -0.5, 2.0);
        PLOT_2D("nu_mc_TransverseEnergy",            "nu_VisibleLongitudinalEnergy",   pNtuple_Neutrino_0, 0.0, 0.01, 0.0, 2.0);
        
        QUICK_PLOT_2D("nu_mc_VertexX",               "nu_VertexX",                   pNtuple_Neutrino_3);
        QUICK_PLOT_2D("nu_mc_VertexY",               "nu_VertexY",                   pNtuple_Neutrino_3);
        QUICK_PLOT_2D("nu_mc_VertexZ",               "nu_VertexZ",                   pNtuple_Neutrino_3);

        PLOT_2D("primary_mc_KineticEnergy",          "primary_KineticEnergy",        pNtuple_Primary_0, 0.0, 1.0, 0.0, 1.0);
        QUICK_PLOT_2D("primary_mc_VertexX",          "primary_VertexX",              pNtuple_Primary_0);
        QUICK_PLOT_2D("primary_mc_VertexY",          "primary_VertexY",              pNtuple_Primary_0);
        QUICK_PLOT_2D("primary_mc_VertexZ",          "primary_VertexZ",              pNtuple_Primary_0);
        QUICK_PLOT_2D("primary_mc_DirectionCosineX", "primary_DirectionCosineX",     pNtuple_Primary_1);
        QUICK_PLOT_2D("primary_mc_DirectionCosineY", "primary_DirectionCosineY",     pNtuple_Primary_1);
        QUICK_PLOT_2D("primary_mc_DirectionCosineZ", "primary_DirectionCosineZ",     pNtuple_Primary_1);
        
        PLOT_2D("cr_mc_KineticEnergy",               "cr_KineticEnergy",             pNtuple_CosmicRay_0, 0.0, 50.0, 0.0, 50.0);
        QUICK_PLOT_2D("cr_mc_VertexX",               "cr_VertexX",                   pNtuple_CosmicRay_0);
        QUICK_PLOT_2D("cr_mc_VertexY",               "cr_VertexY",                   pNtuple_CosmicRay_0);
        QUICK_PLOT_2D("cr_mc_VertexZ",               "cr_VertexZ",                   pNtuple_CosmicRay_0);
        QUICK_PLOT_2D("cr_mc_DirectionCosineX",      "cr_DirectionCosineX",          pNtuple_CosmicRay_1);
        QUICK_PLOT_2D("cr_mc_DirectionCosineY",      "cr_DirectionCosineY",          pNtuple_CosmicRay_1);
        QUICK_PLOT_2D("cr_mc_DirectionCosineZ",      "cr_DirectionCosineZ",          pNtuple_CosmicRay_1);
    }
    
    const float type_McProton     = static_cast<float>(type_RecoProton_McProton + type_RecoPionMuon_McProton + type_RecoOtherTrack_McProton + type_RecoShower_McProton);
    const float type_McPionMuon   = static_cast<float>(type_RecoProton_McPionMuon + type_RecoPionMuon_McPionMuon + type_RecoOtherTrack_McPionMuon + type_RecoShower_McPionMuon);
    const float type_McOtherTrack = static_cast<float>(type_RecoProton_McOtherTrack + type_RecoPionMuon_McOtherTrack + type_RecoOtherTrack_McOtherTrack + type_RecoShower_McOtherTrack);
    const float type_McShower     = static_cast<float>(type_RecoProton_McShower + type_RecoPionMuon_McShower + type_RecoOtherTrack_McShower + type_RecoShower_McShower);
    

    std::cout << " ---- TYPE ----------------------- " << std::endl;
    std::cout << "Reco proton,       MC proton      = " << 100.f * static_cast<float>(type_RecoProton_McProton)         / type_McProton     << "% " << type_RecoProton_McProton << std::endl;
    std::cout << "Reco proton,       MC pion/muon   = " << 100.f * static_cast<float>(type_RecoProton_McPionMuon)       / type_McPionMuon   << "% " << type_RecoProton_McPionMuon << std::endl;
    std::cout << "Reco proton,       MC other track = " << 100.f * static_cast<float>(type_RecoProton_McOtherTrack)     / type_McOtherTrack << "% " << type_RecoProton_McOtherTrack << std::endl;
    std::cout << "Reco proton,       MC shower      = " << 100.f * static_cast<float>(type_RecoProton_McShower)         / type_McShower     << "% " << type_RecoProton_McShower << std::endl;

    std::cout << "Reco pion/muon,    MC proton      = " << 100.f * static_cast<float>(type_RecoPionMuon_McProton)       / type_McProton     << "% " << type_RecoPionMuon_McProton << std::endl;
    std::cout << "Reco pion/muon,    MC pion/muon   = " << 100.f * static_cast<float>(type_RecoPionMuon_McPionMuon)     / type_McPionMuon   << "% " << type_RecoPionMuon_McPionMuon << std::endl;
    std::cout << "Reco pion/muon,    MC other track = " << 100.f * static_cast<float>(type_RecoPionMuon_McOtherTrack)   / type_McOtherTrack << "% " << type_RecoPionMuon_McOtherTrack << std::endl;
    std::cout << "Reco pion/muon,    MC shower      = " << 100.f * static_cast<float>(type_RecoPionMuon_McShower)       / type_McShower     << "% " << type_RecoPionMuon_McShower << std::endl;

    std::cout << "Reco other track,  MC proton      = " << 100.f * static_cast<float>(type_RecoOtherTrack_McProton)     / type_McProton     << "% " << type_RecoOtherTrack_McProton << std::endl;
    std::cout << "Reco other track,  MC pion/muon   = " << 100.f * static_cast<float>(type_RecoOtherTrack_McPionMuon)   / type_McPionMuon   << "% " << type_RecoOtherTrack_McPionMuon << std::endl;
    std::cout << "Reco other track,  MC other track = " << 100.f * static_cast<float>(type_RecoOtherTrack_McOtherTrack) / type_McOtherTrack << "% " << type_RecoOtherTrack_McOtherTrack << std::endl;
    std::cout << "Reco other track,  MC shower      = " << 100.f * static_cast<float>(type_RecoOtherTrack_McShower)     / type_McShower     << "% " << type_RecoOtherTrack_McShower << std::endl;

    std::cout << "Reco shower,       MC proton      = " << 100.f * static_cast<float>(type_RecoShower_McProton)         / type_McProton     << "% " << type_RecoShower_McProton << std::endl;
    std::cout << "Reco shower,       MC pion/muon   = " << 100.f * static_cast<float>(type_RecoShower_McPionMuon)       / type_McPionMuon   << "% " << type_RecoShower_McPionMuon << std::endl;
    std::cout << "Reco shower,       MC other track = " << 100.f * static_cast<float>(type_RecoShower_McOtherTrack)     / type_McOtherTrack << "% " << type_RecoShower_McOtherTrack << std::endl;
    std::cout << "Reco shower,       MC shower      = " << 100.f * static_cast<float>(type_RecoShower_McShower)         / type_McShower     << "% " << type_RecoShower_McShower << std::endl;

}
