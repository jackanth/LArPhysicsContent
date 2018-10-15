/**
 *  @file   larphysicscontent/LArNtuple/LArNtuple.cc
 *
 *  @brief  Implementation of the lar ntuple class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArNtuple/LArNtuple.h"
#include "larphysicscontent/LArHelpers/LArNtupleHelper.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{
void LArNtuple::AddScalarRecord(const LArNtupleRecord &record)
{
    if (!record.WriteToNtuple())
        return;

    if (m_addressesSet)
        this->ValidateAndAddRecord(m_scalarBranchMap, record);

    else
    {
        // Add the record and check that one does not already exist by the same name (this includes type-checking)
        if (!m_scalarBranchMap.emplace(record.BranchName(), LArBranchPlaceholder(record)).second)
        {
            std::cerr << "LArNtuple: cannot add multiple scalar records with the same branch name '" << record.BranchName() << "'" << std::endl;
            throw STATUS_CODE_NOT_ALLOWED;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArNtuple::AddVectorRecordElement(const LArNtupleRecord &record, const LArNtupleHelper::VECTOR_BRANCH_TYPE type)
{
    if (!record.WriteToNtuple())
        return;

    BranchMap &branchMap = this->GetVectorBranchMap(type);

    // If the vector elements are locked in, reuse the scalar record validation mechanics
    if (m_areVectorElementsLocked || m_addressesSet)
        this->ValidateAndAddRecord(branchMap, record);

    else
    {
        // Add the record and check that one does not already exist by the same name (this includes type-checking)
        if (!branchMap.emplace(record.BranchName(), LArBranchPlaceholder(record)).second)
        {
            std::cerr << "LArNtuple: cannot add multiple scalar records with the same branch name '" << record.BranchName() << "'" << std::endl;
            throw STATUS_CODE_NOT_ALLOWED;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArNtuple::FillVectors(const LArNtupleHelper::VECTOR_BRANCH_TYPE type)
{
    BranchMap &branchMap = this->GetVectorBranchMap(type);

    for (auto &entry : branchMap)
    {
        LArBranchPlaceholder &branchPlaceholder = entry.second;

        if (!branchPlaceholder.GetNtupleScalarRecord()) // the branch has not been filled
        {
            std::cerr << "LArNtuple: Could not fill vectors as the branch '" << entry.first << "' has not been populated" << std::endl;
            throw STATUS_CODE_NOT_ALLOWED;
        }

        branchPlaceholder.PushNtupleScalarRecord();
    }

    if (!branchMap.empty()) // if ther has been at least one particle, then we can lock the vector
        m_areVectorElementsLocked = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArNtuple::PushVectors(const LArNtupleHelper::VECTOR_BRANCH_TYPE type)
{
    if (!m_addressesSet)
    {
        m_vectorBranchMaps.emplace(type, std::move(m_vectorElementBranchMap));
        m_vectorElementBranchMap.clear();
        m_areVectorElementsLocked = false;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArNtuple::Fill()
{
    // Fill the ntuple
    std::size_t numBranches = this->SetScalarBranchAddresses();
    numBranches += this->SetVectorBranchAddresses();

    if (numBranches == 0UL)
    {
        std::cerr << "LArNtuple: Could not fill because no branches were set" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    if (m_pOutputTree->Fill() < 0)
    {
        std::cerr << "LArNtuple: Error filling TTree" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    // Prepare the ntuple for the next event
    m_addressesSet = true;
    m_ntupleEmpty  = false;
    this->Reset();
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArNtuple::LArNtuple(const std::string &filePath, const std::string &treeName, const std::string &treeTitle, const bool appendMode) :
    m_pOutputTFile(nullptr),
    m_pOutputTree(nullptr),
    m_scalarBranchMap(),
    m_vectorElementBranchMap(),
    m_vectorBranchMaps(),
    m_addressesSet(false),
    m_ntupleEmpty(true),
    m_areVectorElementsLocked(false),
    m_trackSlidingFitWindow(25U),
    m_cache(),
    m_cacheMCParticles(),
    m_cacheMCCosmics(),
    m_cacheMCPrimaries(),
    m_cacheMCNeutrinos(),
    m_cacheDownstreamThreeDHits(),
    m_cacheDownstreamUHits(),
    m_cacheDownstreamVHits(),
    m_cacheDownstreamWHits(),
    m_cacheDownstreamPfos(),
    m_cacheTrackFits(),
    m_spRegistry(new LArRootRegistry(filePath, appendMode ? LArRootRegistry::FILE_MODE::APPEND : LArRootRegistry::FILE_MODE::NEW))
{
    this->InstantiateTTree(appendMode, treeName, treeTitle, filePath);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArNtuple::Reset()
{
    // Reset the ntuple state to allow recovery from internal errors
    m_vectorElementBranchMap.clear();
    m_areVectorElementsLocked = false;

    if (m_addressesSet)
    {
        // We want to keep the structure but lose the shared pointers for both scalar and vector branches
        for (auto &entry : m_scalarBranchMap)
            entry.second.ClearNtupleRecords();

        for (auto &mapPair : m_vectorBranchMaps)
        {
            for (auto &entry : mapPair.second)
                entry.second.ClearNtupleRecords();
        }
    }

    else
    {
        m_scalarBranchMap.clear();
        m_vectorBranchMaps.clear();
        m_pOutputTree->ResetBranchAddresses();
        m_cache.clear();
    }

    // Clear the transient caches.
    m_cacheMCParticles.clear();
    m_cacheMCCosmics.clear();
    m_cacheMCPrimaries.clear();
    m_cacheMCNeutrinos.clear();
    m_cacheDownstreamThreeDHits.clear();
    m_cacheDownstreamUHits.clear();
    m_cacheDownstreamVHits.clear();
    m_cacheDownstreamWHits.clear();
    m_cacheDownstreamPfos.clear();
    m_cacheTrackFits.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::size_t LArNtuple::SetScalarBranchAddresses()
{
    std::size_t numBranches(0UL);

    for (auto &entry : m_scalarBranchMap)
    {
        const auto spRecord = entry.second.GetNtupleScalarRecord();

        if (!spRecord) // the branch has not been filled
        {
            std::cerr << "LArNtuple: Could not fill ntuple as the branch '" << entry.first << "' has not been populated" << std::endl;
            throw STATUS_CODE_NOT_ALLOWED;
        }

        switch (spRecord->ValueType())
        {
            case LArNtupleRecord::VALUE_TYPE::R_FLOAT:
                this->PushScalarToBranch<LArNtupleRecord::RFloat>(entry.first, entry.second, false);
                break;

            case LArNtupleRecord::VALUE_TYPE::R_INT:
                this->PushScalarToBranch<LArNtupleRecord::RInt>(entry.first, entry.second, false);
                break;

            case LArNtupleRecord::VALUE_TYPE::R_BOOL:
                this->PushScalarToBranch<LArNtupleRecord::RBool>(entry.first, entry.second, false);
                break;

            case LArNtupleRecord::VALUE_TYPE::R_UINT:
                this->PushScalarToBranch<LArNtupleRecord::RUInt>(entry.first, entry.second, false);
                break;

            case LArNtupleRecord::VALUE_TYPE::R_ULONG64:
                this->PushScalarToBranch<LArNtupleRecord::RULong64>(entry.first, entry.second, false);
                break;

            case LArNtupleRecord::VALUE_TYPE::R_TSTRING:
                this->PushScalarToBranch<LArNtupleRecord::RTString>(entry.first, entry.second, false);
                break;

            case LArNtupleRecord::VALUE_TYPE::R_FLOAT_VECTOR:
                this->PushScalarToBranch<LArNtupleRecord::RFloatVector>(entry.first, entry.second, false);
                break;

            case LArNtupleRecord::VALUE_TYPE::R_INT_VECTOR:
                this->PushScalarToBranch<LArNtupleRecord::RIntVector>(entry.first, entry.second, false);
                break;

            default:
                std::cerr << "LArNtuple: Unknown value type" << std::endl;
                throw STATUS_CODE_FAILURE;
        }

        ++numBranches;
    }

    return numBranches;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::size_t LArNtuple::SetVectorBranchAddresses()
{
    std::size_t numBranches(0UL);

    for (auto &mapPair : m_vectorBranchMaps)
    {
        for (auto &entry : mapPair.second)
        {
            LArBranchPlaceholder &branchPlaceholder = entry.second;

            switch (branchPlaceholder.ValueType())
            {
                case LArNtupleRecord::VALUE_TYPE::R_FLOAT:
                    this->PushVectorToBranch<LArNtupleRecord::RFloat>(entry.first, branchPlaceholder, false);
                    break;

                case LArNtupleRecord::VALUE_TYPE::R_INT:
                    this->PushVectorToBranch<LArNtupleRecord::RInt>(entry.first, branchPlaceholder, false);
                    break;

                case LArNtupleRecord::VALUE_TYPE::R_BOOL:
                    this->PushVectorToBranch<LArNtupleRecord::RBool>(entry.first, branchPlaceholder, false);
                    break;

                case LArNtupleRecord::VALUE_TYPE::R_UINT:
                    this->PushVectorToBranch<LArNtupleRecord::RUInt>(entry.first, branchPlaceholder, false);
                    break;

                case LArNtupleRecord::VALUE_TYPE::R_ULONG64:
                    this->PushVectorToBranch<LArNtupleRecord::RULong64>(entry.first, branchPlaceholder, false);
                    break;

                case LArNtupleRecord::VALUE_TYPE::R_TSTRING:
                    this->PushVectorToBranch<LArNtupleRecord::RTString>(entry.first, branchPlaceholder, false);
                    break;

                case LArNtupleRecord::VALUE_TYPE::R_FLOAT_VECTOR:
                    this->PushVectorToBranch<LArNtupleRecord::RFloatVector>(entry.first, branchPlaceholder, false);
                    break;

                case LArNtupleRecord::VALUE_TYPE::R_INT_VECTOR:
                    this->PushVectorToBranch<LArNtupleRecord::RIntVector>(entry.first, branchPlaceholder, false);
                    break;

                default:
                    std::cerr << "LArNtuple: Unknown value type" << std::endl;
                    throw STATUS_CODE_FAILURE;
            }

            ++numBranches;
        }
    }

    return numBranches;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArNtuple::ValidateAndAddRecord(BranchMap &branchMap, const LArNtupleRecord &record)
{
    const auto findIter = branchMap.find(record.BranchName());

    if (findIter == branchMap.end())
    {
        std::cerr << "LArNtuple: cannot add vector record element with branch name '" << record.BranchName()
                  << "' as it was not present in previous ntuple fills" << std::endl;

        throw STATUS_CODE_NOT_ALLOWED;
    }

    if (const LArBranchPlaceholder::NtupleRecordSPtr spCurrentRecord = findIter->second.GetNtupleScalarRecord()) // the branch is already filled
    {
        // Allow overwriting of records that aren't to be written to the ntuple
        if (spCurrentRecord->WriteToNtuple())
        {
            std::cerr << "LArNtuple: cannot add vector record element with branch name '" << record.BranchName()
                      << "' as it has already been populated since the last fill" << std::endl;

            throw STATUS_CODE_NOT_ALLOWED;
        }
    }

    // Add the record (includes type-checking)
    findIter->second.SetNtupleScalarRecord(record);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArNtuple::InstantiateTFile(const bool appendMode, const std::string &filePath)
{
    try
    {
        m_pOutputTFile = new TFile(filePath.c_str(), appendMode ? "UPDATE" : "NEW");
    }

    catch (...)
    {
        std::cerr << "LArNtuple: Failed to open file at " << filePath << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    if (!m_pOutputTFile->IsOpen())
    {
        std::cerr << "LArNtuple: Failed to open file at " << filePath << std::endl;
        throw STATUS_CODE_FAILURE;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArNtuple::InstantiateTTree(const bool appendMode, const std::string &treeName, const std::string &treeTitle, const std::string &filePath)
{
    try
    {
        m_spRegistry->DoAsRegistry([&]() {
            // For our branch mechanics to work, we need to know if the TTree exists and whether it's empty or not
            if (appendMode && m_pOutputTFile->GetListOfKeys()->Contains(treeName.c_str()))
            {
                m_pOutputTree     = dynamic_cast<TTree *>(m_pOutputTFile->Get(treeName.c_str()));
                TObjArray *pArray = m_pOutputTree->GetListOfBranches();

                if (pArray && pArray->GetEntries() > 0)
                {
                    m_ntupleEmpty = false;
                    std::cout << "LArNtuple: Appending new data to non-empty TTree '" << treeName << "' at " << filePath << std::endl;
                }

                else
                    std::cout << "LArNtuple: Appending new data to empty TTree '" << treeName << "' at " << filePath << std::endl;
            }

            else
            {
                m_pOutputTree = new TTree(treeName.c_str(), treeTitle.c_str());
                std::cout << "LArNtuple: Writing data to new TTree '" << treeName << "' at " << filePath << std::endl;
            }
        });
    }

    catch (...)
    {
        std::cerr << "LArNtuple: Failed to find/instantiate TTree with name '" << treeName << "' and title '" << treeTitle << "' at "
                  << filePath << std::endl;
        throw STATUS_CODE_FAILURE;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitList LArNtuple::GetAllDownstreamTwoDHits(const ParticleFlowObject *const pPfo) const
{
    // Get all the (possibly-cached) U-, V-, and W-hits and add them all up.
    CaloHitList hitsU(this->GetAllDownstreamUHits(pPfo));
    CaloHitList hitsV(this->GetAllDownstreamVHits(pPfo));
    CaloHitList hitsW(this->GetAllDownstreamWHits(pPfo));

    CaloHitList hitsTwoD;
    hitsTwoD.insert(hitsTwoD.end(), std::make_move_iterator(hitsU.begin()), std::make_move_iterator(hitsU.end()));
    hitsTwoD.insert(hitsTwoD.end(), std::make_move_iterator(hitsV.begin()), std::make_move_iterator(hitsV.end()));
    hitsTwoD.insert(hitsTwoD.end(), std::make_move_iterator(hitsW.begin()), std::make_move_iterator(hitsW.end()));

    return hitsTwoD;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const PfoList &LArNtuple::GetAllDownstreamPfos(const ParticleFlowObject *const pPfo) const
{
    return this->CacheWrapper<PfoList>(pPfo, m_cacheDownstreamPfos, [&]() {
        PfoList downstreamPfos;
        LArPfoHelper::GetAllDownstreamPfos(pPfo, downstreamPfos);
        return downstreamPfos;
    });
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHitList &LArNtuple::GetAllDownstreamHitsImpl(const ParticleFlowObject *const pPfo, const HitType hitType, PfoCache<CaloHitList> &cache) const
{
    return this->CacheWrapper<CaloHitList>(pPfo, cache, [&]() {
        const PfoList &downstreamPfos = this->GetAllDownstreamPfos(pPfo);

        CaloHitList caloHitList;
        LArPfoHelper::GetCaloHits(downstreamPfos, hitType, caloHitList);

        return caloHitList;
    });
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArBranchPlaceholder::NtupleRecordSPtr LArNtuple::GetScalarRecord(const std::string &branchName) const
{
    const LArBranchPlaceholder &                 branchPlaceholder = this->GetBranchPlaceholder(m_scalarBranchMap, branchName);
    const LArBranchPlaceholder::NtupleRecordSPtr spRecord          = branchPlaceholder.GetNtupleScalarRecord();

    if (!spRecord)
    {
        std::cerr << "LArNtuple: Found record by name '" << branchName << "' but it had not yet been filled" << std::endl;
        throw pandora::STATUS_CODE_FAILURE;
    }

    return spRecord;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitList LArNtuple::GetAllTwoDHits(const ParticleFlowObject *const pPfo) const
{
    CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, caloHitList);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, caloHitList);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, caloHitList);

    return caloHitList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const LArBranchPlaceholder &LArNtuple::GetBranchPlaceholder(const LArNtupleHelper::VECTOR_BRANCH_TYPE type, const std::string &branchName) const
{
    const auto branchMapFindIter = m_vectorBranchMaps.find(type);

    if (branchMapFindIter == m_vectorBranchMaps.end())
    {
        std::cerr << "LArNtuple: Could not find vector record set of the requested type for branch '" << branchName << "'" << std::endl;
        throw pandora::STATUS_CODE_FAILURE;
    }

    return this->GetBranchPlaceholder(branchMapFindIter->second, branchName);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const LArBranchPlaceholder &LArNtuple::GetBranchPlaceholder(const BranchMap &branchMap, const std::string &branchName) const
{
    const auto findIter = branchMap.find(branchName);

    if (findIter == m_scalarBranchMap.end())
    {
        std::cerr << "LArNtuple: Could not find record by name '" << branchName << "'" << std::endl;
        throw pandora::STATUS_CODE_FAILURE;
    }

    return findIter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArNtupleHelper::TrackFitSharedPtr LArNtuple::CalculateTrackFit(
    const Pandora &pandoraInstance, const ParticleFlowObject *const pPfo, const unsigned int slidingFitWindow) const
{
    // Get the 3D clusters and make sure there's at least one
    ClusterList threeDClusterList;
    LArPfoHelper::GetClusters(pPfo, TPC_3D, threeDClusterList);

    if (threeDClusterList.empty())
        return nullptr;

    // Get the complete 3D coordinate vector and make sure it's not too small
    CartesianPointVector coordinateVector;

    for (const Cluster *const pCluster : threeDClusterList)
    {
        CartesianPointVector clusterCoordinateVector;
        LArClusterHelper::GetCoordinateVector(pCluster, clusterCoordinateVector);

        coordinateVector.insert(coordinateVector.end(), std::make_move_iterator(clusterCoordinateVector.begin()),
            std::make_move_iterator(clusterCoordinateVector.end()));
    }

    if (coordinateVector.size() < 3UL)
        return nullptr;

    const float layerPitch(LArGeometryHelper::GetWireZPitch(pandoraInstance));

    // If the fit fails, just return a nullptr
    try
    {
        return LArNtupleHelper::TrackFitSharedPtr(new ThreeDSlidingFitResult(&coordinateVector, slidingFitWindow, layerPitch));
    }

    catch (...)
    {
    }

    return nullptr;
}
} // namespace lar_physics_content