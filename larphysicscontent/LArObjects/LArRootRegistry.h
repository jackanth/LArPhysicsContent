/**
 *  @file   larphysicscontent/LArObjects/LArRootRegistry.h
 *
 *  @brief  Header file for the lar ROOT registry class.
 *
 *  $Log: $
 */
#ifndef LAR_ROOT_REGISTRY_H
#define LAR_ROOT_REGISTRY_H 1

#include "Pandora/StatusCodes.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TSystem.h"

#include <functional>
#include <iostream>
#include <unordered_set>

namespace lar_physics_content
{

/**
 *  @brief  LArRootRegistry class
 */
class LArRootRegistry
{
public:
    /**
     *  @brief  ROOT file opening modes
     */
    enum class FILE_MODE
    {
        NEW,      ///< Write a new file
        APPEND,   ///< Append to a file
        OVERWRITE ///< Overwrite the file
    };

    /**
     *  @brief  Constructor
     *
     *  @param  filePath the ROOT file path
     *  @param  fileMode the file mode
     */
    LArRootRegistry(const std::string &filePath, const FILE_MODE fileMode);

    /**
     * @brief  Deleted copy constructor
     */
    LArRootRegistry(const LArRootRegistry &) = delete;

    /**
     * @brief  Default move constructor
     */
    LArRootRegistry(LArRootRegistry &&) = default;

    /**
     * @brief  Deleted copy assignment operator
     */
    LArRootRegistry &operator=(const LArRootRegistry &) = delete;

    /**
     * @brief  Default move assignment operator
     */
    LArRootRegistry &operator=(LArRootRegistry &&) = default;

    /**
     * @brief  Destructor
     */
    ~LArRootRegistry();

    /**
     *  @brief  Do something with directory changed to the registry directory
     *
     *  @param  fn the function to run
     */
    void DoAsRegistry(const std::function<void(void)> &fn) const;

    /**
     *  @brief  Create a ROOT object in this registry with a unique name
     *
     *  @param  name the name (a prefix will be added)
     *  @param  args the constructor arguments to forward
     *
     *  @return pointer to the object
     */
    template <typename T, typename... TARGS>
    std::decay_t<T> *CreateWithUniqueName(const std::string &name, TARGS &&... args);

    /**
     *  @brief  Create a ROOT object in this registry
     *
     *  @param  args the constructor arguments to forward
     *
     *  @return pointer to the object
     */
    template <typename T, typename... TARGS>
    std::decay_t<T> *Create(TARGS &&... args);

    /**
     *  @brief  Delete all objects created by the registry - in memory and on file.
     */
    void Clear();

    /**
     *  @brief  Delete all in-memory objects created by the registry.
     */
    void ClearMemory();

    /**
     *  @brief  Write all in-memory objects to the file.
     */
    void Write() const;

    /**
     *  @brief  Get the TFile
     * 
     *  @return address of the TFile
     */
    TFile * GetTFile();

private:
    TFile *                       m_pFile;           ///< Address of this registry's TFile
    static std::size_t            m_objectNameCount; ///< The object name count for creating unique names
    std::unordered_set<TObject *> m_objectList;      ///< The list of objects

    /**
     *  @brief  Change directory if the pointer is not null
     *
     *  @param  pDir address of the new directory
     */
    void ChangeDirectory(TDirectory *pDir) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T, typename... TARGS>
inline std::decay_t<T> *LArRootRegistry::CreateWithUniqueName(const std::string &name, TARGS &&... args)
{
    return this->Create<std::decay_t<T>>((std::to_string(m_objectNameCount++) + "_" + name).c_str(), std::forward<TARGS>(args)...);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T, typename... TARGS>
inline std::decay_t<T> *LArRootRegistry::Create(TARGS &&... args)
{
    std::decay_t<T> *pObj(nullptr);
    this->DoAsRegistry([&]() { pObj = new std::decay_t<T>(std::forward<TARGS>(args)...); });

    m_objectList.insert(static_cast<TObject *>(pObj));
    return pObj;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArRootRegistry::Clear()
{
    if (m_pFile && m_pFile->IsOpen())
    {
        for (TObject *pObject : m_objectList)
        {
            if (pObject)
            {
                delete pObject;
                pObject = nullptr;
            }
        }

        m_objectList.clear();
        m_pFile->Delete("*;*");
        m_pFile->DeleteAll();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArRootRegistry::ClearMemory()
{
    if (m_pFile && m_pFile->IsOpen())
    {
        for (TObject *pObject : m_objectList)
        {
            if (pObject)
            {
                delete pObject;
                pObject = nullptr;
            }
        }

        m_objectList.clear();
        m_pFile->Delete("*");
        m_pFile->DeleteAll();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArRootRegistry::Write() const
{
    if (m_pFile && m_pFile->IsOpen())
        m_pFile->Write();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TFile * LArRootRegistry::GetTFile()
{
    return m_pFile;
}

} // namespace lar_physics_content

#endif // #ifndef LAR_ROOT_REGISTRY_H