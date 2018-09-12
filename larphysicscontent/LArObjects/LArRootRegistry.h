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

#include <iostream>
#include <functional>

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
    std::decay_t<T> *CreateWithUniqueName(const std::string &name, TARGS &&... args) const;

    /**
     *  @brief  Create a ROOT object in this registry
     *
     *  @param  args the constructor arguments to forward
     *
     *  @return pointer to the object
     */
    template <typename T, typename... TARGS>
    std::decay_t<T> *Create(TARGS &&... args) const;

    /**
     *  @brief  Delete all objects created by the registry.
     */
    void Clear() const;

private:
    TFile *            m_pFile;           ///< Address of this registry's TFile
    static std::size_t m_objectNameCount; ///< The object name count for creating unique names

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
inline std::decay_t<T> *LArRootRegistry::CreateWithUniqueName(const std::string &name, TARGS &&... args) const
{
    return this->Create<std::decay_t<T>>((std::to_string(m_objectNameCount++) + "_" + name).c_str(), std::forward<TARGS>(args)...);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T, typename... TARGS>
inline std::decay_t<T> *LArRootRegistry::Create(TARGS &&... args) const
{
    std::decay_t<T> *pObj(nullptr);
    this->DoAsRegistry([&]() { pObj = new std::decay_t<T>(std::forward<TARGS>(args)...); });
    return pObj;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArRootRegistry::Clear() const
{
    if (m_pFile && m_pFile->IsOpen())
        m_pFile->Delete("*;*");
}
} // namespace lar_physics_content

#endif // #ifndef LAR_ROOT_REGISTRY_H