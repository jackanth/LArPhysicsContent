/**
 *  @file   larphysicscontent/LArNtuple/LArNtuple.cc
 *
 *  @brief  Implementation of the lar ntuple class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArNtuple/LArNtuple.h"
#include "larphysicscontent/LArHelpers/LArNtupleHelper.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{

LArNtuple::~LArNtuple()
{
    // Write any remaining data in the TTree
    if (m_pOutputTree)
        m_pOutputTree->Write();

    // If we have a TFile, close and delete it - the TTree is then deleted by its owning TFile
    if (m_pOutputTFile)
    {
        if (m_pOutputTFile->IsOpen())
            m_pOutputTFile->Close();

        delete m_pOutputTFile;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArNtuple::AddScalarRecord(const LArNtupleRecord &record)
{
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

void LArNtuple::AddVectorRecordElement(const LArNtupleRecord &record)
{
    BranchMap &branchMap = this->GetCurrentVectorBranchMap();

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

void LArNtuple::FillVectors()
{
    BranchMap &branchMap = this->GetCurrentVectorBranchMap();

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

void LArNtuple::PushVectors()
{
    // Tick over the vector stack
    if (m_addressesSet)
        ++m_vectorBranchMapIndex;

    else
    {
        m_vectorBranchMaps.emplace_back(std::move(m_vectorElementBranchMap));
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
    m_vectorBranchMapIndex(0UL),
    m_addressesSet(false),
    m_ntupleEmpty(true),
    m_areVectorElementsLocked(false),
    m_cache(),
    m_cacheMCParticles(),
    m_cacheMCCosmics(),
    m_cacheMCPrimaries(),
    m_cacheMCNeutrinos(),
    m_cacheDownstreamThreeDHits(),
    m_cacheDownstreamUHits(),
    m_cacheDownstreamVHits(),
    m_cacheDownstreamWHits()
{
    this->InstantiateTFile(appendMode, filePath);
    this->InstantiateTTree(appendMode, treeName, treeTitle, filePath); // the TTree will now be 'owned' by the TFile
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArNtuple::GetMCParticle(const ParticleFlowObject *const pPfo, const MCParticleList *const pMCParticleList)
{
    if (LArNtupleHelper::GetParticleClass(pPfo) == LArNtupleHelper::PARTICLE_CLASS::NEUTRINO)
        return this->GetMCNeutrino(pPfo, pMCParticleList);

    const CaloHitList caloHitList = this->GetAllTwoDHits(pPfo);
    return this->GetMCParticleImpl(caloHitList, [](const MCParticle *const pMCParticle) { return pMCParticle; });
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArNtuple::GetMCNeutrino(const ParticleFlowObject *const pPfo, const MCParticleList *const pMCParticleList)
{
    // If we only have 0 or 1 MC neutrinos, then the answer is obvious
    MCParticleVector mcNeutrinoVector;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, mcNeutrinoVector);

    if (mcNeutrinoVector.empty())
        return nullptr;

    if (mcNeutrinoVector.size() == 1UL)
        return mcNeutrinoVector.front();

    // We have multiple MC neutrinos, so find the most fitting one
    const CaloHitList caloHitList = this->GetAllDownstreamTwoDHits(pPfo);

    return this->GetMCParticleImpl(caloHitList, [](const MCParticle *const pMCParticle) {
        const MCParticle *const pParent = LArMCParticleHelper::GetParentMCParticle(pMCParticle);
        return LArMCParticleHelper::IsNeutrino(pParent) ? pParent : nullptr; // this means 'is beam neutrino?'
    });
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArNtuple::Reset()
{
    // Reset the ntuple state to allow recovery from internal errors
    m_vectorElementBranchMap.clear();
    m_areVectorElementsLocked = false;
    m_vectorBranchMapIndex    = 0UL;

    if (m_addressesSet)
    {
        // We want to keep the structure but lose the shared pointers for both scalar and vector branches
        for (auto &entry : m_scalarBranchMap)
            entry.second.ClearNtupleRecords();

        for (auto &branchMap : m_vectorBranchMaps)
        {
            for (auto &entry : branchMap)
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
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArNtuple::BranchMap &LArNtuple::GetCurrentVectorBranchMap()
{
    if (!m_addressesSet)
        return m_vectorElementBranchMap;

    if (m_vectorBranchMapIndex >= m_vectorBranchMaps.size())
    {
        std::cerr << "LArNtuple: The vector branch map index exceeded the map size" << std::endl;
        throw STATUS_CODE_NOT_ALLOWED;
    }

    return m_vectorBranchMaps.at(m_vectorBranchMapIndex);
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
                this->PushScalarToBranch<LArNtupleRecord::RFloat>(entry.first, entry.second, true);
                break;

            case LArNtupleRecord::VALUE_TYPE::R_INT:
                this->PushScalarToBranch<LArNtupleRecord::RInt>(entry.first, entry.second, true);
                break;

            case LArNtupleRecord::VALUE_TYPE::R_BOOL:
                this->PushScalarToBranch<LArNtupleRecord::RBool>(entry.first, entry.second, true);
                break;

            case LArNtupleRecord::VALUE_TYPE::R_UINT:
                this->PushScalarToBranch<LArNtupleRecord::RUInt>(entry.first, entry.second, true);
                break;

            case LArNtupleRecord::VALUE_TYPE::R_ULONG64:
                this->PushScalarToBranch<LArNtupleRecord::RULong64>(entry.first, entry.second, true);
                break;

            case LArNtupleRecord::VALUE_TYPE::R_TSTRING:
                this->PushScalarToBranch<LArNtupleRecord::RTString>(entry.first, entry.second, true);
                break;

            case LArNtupleRecord::VALUE_TYPE::R_FLOAT_VECTOR:
                this->PushScalarToBranch<LArNtupleRecord::RFloatVector>(entry.first, entry.second, true);
                break;

            case LArNtupleRecord::VALUE_TYPE::R_INT_VECTOR:
                this->PushScalarToBranch<LArNtupleRecord::RIntVector>(entry.first, entry.second, true);
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

    for (auto &branchMap : m_vectorBranchMaps)
    {
        for (auto &entry : branchMap)
        {
            LArBranchPlaceholder &branchPlaceholder = entry.second;

            switch (branchPlaceholder.ValueType())
            {
                case LArNtupleRecord::VALUE_TYPE::R_FLOAT:
                    this->PushVectorToBranch<LArNtupleRecord::RFloat>(entry.first, branchPlaceholder, true);
                    break;

                case LArNtupleRecord::VALUE_TYPE::R_INT:
                    this->PushVectorToBranch<LArNtupleRecord::RInt>(entry.first, branchPlaceholder, true);
                    break;

                case LArNtupleRecord::VALUE_TYPE::R_BOOL:
                    this->PushVectorToBranch<LArNtupleRecord::RBool>(entry.first, branchPlaceholder, true);
                    break;

                case LArNtupleRecord::VALUE_TYPE::R_UINT:
                    this->PushVectorToBranch<LArNtupleRecord::RUInt>(entry.first, branchPlaceholder, true);
                    break;

                case LArNtupleRecord::VALUE_TYPE::R_ULONG64:
                    this->PushVectorToBranch<LArNtupleRecord::RULong64>(entry.first, branchPlaceholder, true);
                    break;

                case LArNtupleRecord::VALUE_TYPE::R_TSTRING:
                    this->PushVectorToBranch<LArNtupleRecord::RTString>(entry.first, branchPlaceholder, false);
                    break;

                case LArNtupleRecord::VALUE_TYPE::R_FLOAT_VECTOR:
                    this->PushVectorToBranch<LArNtupleRecord::RFloatVector>(entry.first, branchPlaceholder, true);
                    break;

                case LArNtupleRecord::VALUE_TYPE::R_INT_VECTOR:
                    this->PushVectorToBranch<LArNtupleRecord::RIntVector>(entry.first, branchPlaceholder, true);
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

    if (findIter->second.GetNtupleScalarRecord()) // the branch is already filled
    {
        std::cerr << "LArNtuple: cannot add vector record element with branch name '" << record.BranchName()
                  << "' as it has already been populated since the last fill" << std::endl;

        throw STATUS_CODE_NOT_ALLOWED;
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
    }

    catch (...)
    {
        std::cerr << "LArNtuple: Failed to find/instantiate TTree with name '" << treeName << "' and title '" << treeTitle << "' at "
                  << filePath << std::endl;
        throw STATUS_CODE_FAILURE;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

MCParticleWeightMap LArNtuple::GetMCParticleWeightMap(const CaloHitList &caloHitList, const MCParticleMapFn &mapFn) const
{
    MCParticleWeightMap mcParticleWeightMap;

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        // Synthesize the weights from every hit into the main map
        const MCParticleWeightMap &hitMCParticleWeightMap(pCaloHit->GetMCParticleWeightMap());
        MCParticleVector           mcParticleVector;

        for (const MCParticleWeightMap::value_type &mapEntry : hitMCParticleWeightMap)
            mcParticleVector.push_back(mapEntry.first);

        std::sort(mcParticleVector.begin(), mcParticleVector.end(), PointerLessThan<MCParticle>());

        for (const MCParticle *const pMCParticle : mcParticleVector)
        {
            try
            {
                const MCParticle *const pMCMappedParticle = mapFn(pMCParticle);

                if (!pMCMappedParticle)
                    continue;

                const auto findIter = mcParticleWeightMap.find(pMCMappedParticle);

                if (findIter == mcParticleWeightMap.end())
                    mcParticleWeightMap.emplace(pMCMappedParticle, hitMCParticleWeightMap.at(pMCParticle));

                else
                    findIter->second += hitMCParticleWeightMap.at(pMCParticle);
            }

            catch (...)
            {
                continue;
            }
        }
    }

    return mcParticleWeightMap;
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

const CaloHitList &LArNtuple::GetAllDownstreamHitsImpl(const ParticleFlowObject *const pPfo, const HitType hitType, PfoCache<CaloHitList> &cache) const
{
    // Check the cache first
    const auto findIter = cache.find(pPfo);

    if (findIter != cache.end())
        return findIter->second;

    // Not in cache, so work it out and cache it
    PfoList downstreamPfos;
    LArPfoHelper::GetAllDownstreamPfos(pPfo, downstreamPfos);

    CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(downstreamPfos, hitType, caloHitList);

    return cache.emplace(pPfo, std::move(caloHitList)).first->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArNtuple::GetMCParticleImpl(const CaloHitList &caloHitList, const MCParticleMapFn &mapFn) const
{
    const MCParticleWeightMap mcParticleWeightMap = this->GetMCParticleWeightMap(caloHitList, mapFn);

    float             bestWeight(0.f);
    const MCParticle *pBestMCParticle(nullptr);

    MCParticleVector mcParticleVector;

    for (const MCParticleWeightMap::value_type &mapEntry : mcParticleWeightMap)
        mcParticleVector.push_back(mapEntry.first);

    std::sort(mcParticleVector.begin(), mcParticleVector.end(), PointerLessThan<MCParticle>());

    for (const MCParticle *const pCurrentMCParticle : mcParticleVector)
    {
        const float currentWeight(mcParticleWeightMap.at(pCurrentMCParticle));

        if (currentWeight > bestWeight)
        {
            pBestMCParticle = pCurrentMCParticle;
            bestWeight      = currentWeight;
        }
    }

    return pBestMCParticle; // can be nullptr
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

} // namespace lar_physics_content