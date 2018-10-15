/**
 *  @file   larphysicscontent/LArNtuple/LArRootRegistry.cc
 *
 *  @brief  Implementation of the lar ROOT registry class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArObjects/LArRootRegistry.h"
#include "TROOT.h"

using namespace pandora;

namespace lar_physics_content
{
std::size_t LArRootRegistry::m_objectNameCount = 0UL;

LArRootRegistry::LArRootRegistry(const std::string &filePath, const FILE_MODE fileMode) : m_pFile(nullptr)
{
    // Get the original directory
    TDirectory *pOriginalDir = TDirectory::CurrentDirectory();

    try
    {
        switch (fileMode)
        {
            case FILE_MODE::NEW:
                m_pFile = new TFile(filePath.c_str(), "NEW");
                break;
            case FILE_MODE::APPEND:
                m_pFile = new TFile(filePath.c_str(), "UPDATE");
                break;
            case FILE_MODE::OVERWRITE:
                m_pFile = new TFile(filePath.c_str(), "RECREATE");
                break;
            default:
                throw StatusCodeException(STATUS_CODE_FAILURE);
        }
    }

    catch (...)
    {
        std::cerr << "LArRootRegistry: Failed to open file at " << filePath << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    if (!m_pFile || !m_pFile->IsOpen() || m_pFile->IsZombie() || m_pFile->TestBit(TFile::kRecovered))
    {
        std::cerr << "LArRootRegistry: Failed to open file at " << filePath << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    m_pFile->SetBit(TFile::kDevNull);
    this->ChangeDirectory(pOriginalDir);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArRootRegistry::~LArRootRegistry()
{
    if (m_pFile)
    {
        if (m_pFile->IsOpen())
        {
            m_pFile->Write();
            this->ClearMemory();
            m_pFile->Close();
        }

        gROOT->GetListOfFiles()->Remove(m_pFile);
        delete m_pFile;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArRootRegistry::DoAsRegistry(const std::function<void(void)> &fn) const
{
    if (!m_pFile || !m_pFile->IsOpen())
    {
        std::cerr << "LArRootRegistry: Failed to access internal TFile" << std::endl;
        throw pandora::STATUS_CODE_FAILURE;
    }

    // Change to this registry's directory if required
    TDirectory *pOriginalDir = TDirectory::CurrentDirectory();
    const bool  changeDir    = (pOriginalDir != m_pFile);

    if (changeDir)
        this->ChangeDirectory(m_pFile);

    fn();

    // Change back to the original directory
    if (changeDir)
        this->ChangeDirectory(pOriginalDir);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArRootRegistry::ChangeDirectory(TDirectory *pDir) const
{
    if (pDir && !pDir->cd())
    {
        std::cerr << "LArRootRegistry: Failed to change ROOT directory" << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }
}
} // namespace lar_physics_content
