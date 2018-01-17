/**
 *  @file   LArPhysicsContent/include/LArAnalysisParticle.h
 *
 *  @brief  Header file for the lar analysis particle class.
 *
 *  $Log: $
 */
 
#ifndef LAR_ANALYSIS_PARTICLE_H
#define LAR_ANALYSIS_PARTICLE_H 1

#include "Objects/CartesianVector.h"
#include "Objects/ParticleFlowObject.h"

#include "Pandora/ObjectFactory.h"

#include "DebugDefinitions.h"

using namespace pandora;

namespace lar_physics_content
{

/**
 *  @brief  Forward declaration of the LArAnalysisParticleParameters class
 */
class LArAnalysisParticleParameters;
    
/**
 *  @brief  LArAnalysisParticle class
 */
class LArAnalysisParticle : public ParticleFlowObject
{
public:
    using List = std::list<LArAnalysisParticle>; ///< Alias for a list of LArAnalysisParticles

    /**
     *  @brief  Enum for the type of particle
     */
    enum class TYPE : int
    {
        UNKNOWN   = 0, ///< Unknown particle type
        TRACK     = 1, ///< A track not clearly a pion, muon or proton
        SHOWER    = 2, ///< A shower
        PROTON    = 3, ///< A proton
        PION_MUON = 4  ///< A pion or a muon
    };
    
    using PfoTypeMap = std::unordered_map<const ParticleFlowObject *, TYPE>; ///< Alias for a map from PFOs to types
    
    /**
     *  @brief  TypeTree class
     */
    class TypeTree
    {
    public:
        using List = std::list<TypeTree>; ///< Alias for a list of TypeTrees
        
        /**
         *  @brief  Constructor
         * 
         *  @param  type the particle type
         *  @param  daughterTypes the daughter TypeTrees
         * 
         */
        TypeTree(const TYPE type, TypeTree::List daughterTypes) noexcept;
        
        /**
         *  @brief  Constructor
         * 
         */
        TypeTree() noexcept;
        
        /**
         *  @brief  Constructor
         * 
         *  @param  type the particle type
         * 
         */
        TypeTree(const TYPE type) noexcept;
        
        /**
         *  @brief  Get the particle type
         * 
         *  @return the particle type
         *
         */
        TYPE Type() const noexcept;
        
        /**
         *  @brief  Get the daughter types
         * 
         *  @return the daughter types
         * 
         */
        const TypeTree::List & Daughters() const noexcept;
        
    private:
        TYPE           m_type;          ///< The particle type
        TypeTree::List m_daughterTypes; ///< The list of daughter types
    };
    
    /**
     *  @brief  Constructor
     * 
     *  @param  parameters the parameters with which to construct the object
     * 
     */
    LArAnalysisParticle(const LArAnalysisParticleParameters &parameters);

    /**
     *  @brief  Get the particle type
     * 
     *  @return the particle type
     * 
     */
    TYPE Type() const noexcept;

    /**
     *  @brief  Get the particle type tree
     * 
     *  @return the particle type tree
     * 
     */
    const TypeTree & GetTypeTree() const noexcept;

    /**
     *  @brief  Get the particle energy
     * 
     *  @return the particle energy
     * 
     */
    float AnalysisEnergy() const noexcept;
    
    /**
     *  @brief  Get the particle energy from charge
     * 
     *  @return the particle energy from charge
     * 
     */
    float EnergyFromCharge() const noexcept;

    /**
     *  @brief  Get whether the fiducial cut is satisfied
     * 
     *  @return whether the fiducial cut is satisfied
     * 
     */
    bool IsVertexFiducial() const noexcept;
    
    /**
     *  @brief  Get whether the fiducial cut is satisfied
     * 
     *  @return whether the fiducial cut is satisfied
     * 
     */
    bool AreAllHitsFiducial() const noexcept;

    /**
     *  @brief  Get the vertex position
     * 
     *  @return the vertex position
     * 
     */
    const CartesianVector & VertexPosition() const noexcept;

    /**
     *  @brief  Get the initial direction
     * 
     *  @return the initial direction
     * 
     */
    const CartesianVector & DirectionCosines() const noexcept;
    
    /**
     *  @brief  Get the initial direction
     * 
     *  @return the initial direction
     * 
     */
    const CartesianVector & AnalysisMomentum() const noexcept;

    /**
     *  @brief  Get the number of 3D hits
     * 
     *  @return the number of 3D hits
     * 
     */
    unsigned NumberOf3dHits() const noexcept;

    /**
     *  @brief  Get the number of W hits
     * 
     *  @return the number of W hits
     * 
     */
    unsigned NumberOfCollectionPlaneHits() const noexcept;

    /**
     *  @brief  Get whether the primary is a shower
     * 
     *  @return whether the primary is a shower
     * 
     */
    bool IsShower() const noexcept;
    
    /**
     *  @brief  Get whether the primary is a shower
     * 
     *  @return whether the primary is a shower
     * 
     */
    unsigned NumberOfDownstreamParticles() const noexcept;

    /**
     *  @brief  Get whether the particle has MC info
     * 
     *  @return whether the particle has MC info
     * 
     */
    bool HasMcInfo() const noexcept;

    /**
     *  @brief  Get the MC particle type
     * 
     *  @return the MC particle type
     * 
     */
    TYPE McType() const;

    /**
     *  @brief  Get the MC particle type tree
     * 
     *  @return the MC particle type tree
     * 
     */
    const TypeTree & GetMcTypeTree() const;

    /**
     *  @brief  Get the MC particle energy
     * 
     *  @return the MC particle energy
     * 
     */
    float McEnergy() const;
    
    /**
     *  @brief  Get the MC vertex position
     * 
     *  @return the MC vertex position
     * 
     */
    const CartesianVector & McMomentum() const;

    /**
     *  @brief  Get the MC vertex position
     * 
     *  @return the MC vertex position
     * 
     */
    const CartesianVector & McVertexPosition() const;
    
    /**
     *  @brief  Get the MC direction cosines
     * 
     *  @return the MC direction cosines
     * 
     */
    const CartesianVector & McDirectionCosines() const;
    
    /**
     *  @brief  Get the MC vertex position
     * 
     *  @return the MC vertex position
     * 
     */
    bool McIsVertexFiducial() const;
    
    /**
     *  @brief  Get the MC vertex position
     * 
     *  @return the MC vertex position
     * 
     */
    bool McIsContained() const;

    /**
     *  @brief  Get whether the particle is a shower (MC)
     * 
     *  @return whether the particle is a shower (MC)
     * 
     */
    bool McIsShower() const;

    /**
     *  @brief  Get whether the particle has been reconstructed correctly
     * 
     *  @return whether the particle has been reconstructed correctly
     * 
     */
    bool McIsCorrectlyReconstructed() const;
    
    /**
     *  @brief  Get whether the particle has been reconstructed correctly
     * 
     *  @return whether the particle has been reconstructed correctly
     * 
     */
    int McPdgCode() const;
    
    /**
     *  @brief  Print some information about the particle
     * 
     */
    void Print() const;
        
    /**
     *  @brief  Get the particle type as a string
     * 
     *  @param  type the type
     * 
     *  @return the string type
     * 
     */
    static std::string TypeAsString(const TYPE type);
    
    /**
     *  @brief  Get the particle type tree as a string
     * 
     *  @param  typeTree the typeTree
     * 
     *  @return the string type
     * 
     */
    static std::string TypeTreeAsString(const TypeTree &typeTree);

private:
    TYPE               m_type;                           ///< The particle type
    TypeTree           m_typeTree;                       ///< The type tree
    float              m_analysisEnergy;                 ///< The particle energy in GeV
    float              m_energyFromCharge;               ///< The particle energy in GeV, calculated only using charge
    bool               m_isVertexFiducial;               ///< Whether the vertex is fiducial
    bool               m_areAllHitsFiducial;             ///< Whether all the hits are fiducial
    CartesianVector    m_vertexPosition;                 ///< The vertex position
    CartesianVector    m_directionCosines;               ///< The direction cosines at the vertex
    CartesianVector    m_analysisMomentum;               ///< The momentum at the vertex
    unsigned           m_numberOf3dHits;                 ///< The number of 3D hits
    unsigned           m_numberOfCollectionPlaneHits;    ///< The number of collection-plane hits
    bool               m_isShower;                       ///< Whether the particle is a shower
    unsigned           m_numberOfDownstreamParticles;    ///< The number of downstream particles
    bool               m_hasMcInfo;                      ///< Whether the particle has MC info attached
    TYPE               m_mcType;                         ///< The MC type
    TypeTree           m_mcTypeTree;                     ///< The MC type tree
    float              m_mcEnergy;                       ///< The MC energy
    CartesianVector    m_mcMomentum;                     ///< The MC momentum at the vertex
    CartesianVector    m_mcVertexPosition;               ///< The MC vertex position
    CartesianVector    m_mcDirectionCosines;             ///< The MC direction cosines
    bool               m_mcIsVertexFiducial;             ///< Whether the vertex is fiducial (MC quantity)
    bool               m_mcIsContained;                  ///< Whether the particle is contained (MC quantity)
    bool               m_mcIsShower;                     ///< Whether the particle is a shower (MC quantity)
    bool               m_mcIsCorrectlyReconstructed;     ///< Whether the particle has been reconstructed correctly (MC quantity)
    int                m_mcPdgCode;                      ///< The PDG code of the particle (MC quantity)
    
    /**
     *  @brief  Get the particle type tree as a string (implementation)
     * 
     *  @param  typeTree the typeTree
     *  @param  printTrailingDelimiter whether to print the trailing delimiter
     * 
     *  @return the string type
     * 
     */
    static std::string TypeTreeAsStringImpl(const TypeTree &typeTree, const bool printTrailingDelimiter);
    
    /**
     *  @brief  Print an error message and throw if there is no MC information present
     */
    void ThrowIfNoMcInfo() const;
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LArAnalysisParticle parameters
 */
class LArAnalysisParticleParameters : public object_creation::ParticleFlowObject::Parameters
{
public:
    LArAnalysisParticleParameters() noexcept;

    LArAnalysisParticle::TYPE        m_type;                       
    LArAnalysisParticle::TypeTree    m_typeTree;                   
    float                            m_analysisEnergy;             
    float                            m_energyFromCharge;           
    bool                             m_isVertexFiducial;           
    bool                             m_areAllHitsFiducial;         
    CartesianVector                  m_vertexPosition;             
    CartesianVector                  m_directionCosines;           
    CartesianVector                  m_analysisMomentum;                   
    unsigned                         m_numberOf3dHits;             
    unsigned                         m_numberOfCollectionPlaneHits;
    bool                             m_isShower;                   
    unsigned                         m_numberOfDownstreamParticles;
    bool                             m_hasMcInfo;                  
    LArAnalysisParticle::TYPE        m_mcType;                     
    LArAnalysisParticle::TypeTree    m_mcTypeTree;                 
    float                            m_mcEnergy;                   
    CartesianVector                  m_mcMomentum;                 
    CartesianVector                  m_mcVertexPosition; 
    CartesianVector                  m_mcDirectionCosines;
    bool                             m_mcIsVertexFiducial;         
    bool                             m_mcIsContained;              
    bool                             m_mcIsShower;                 
    bool                             m_mcIsCorrectlyReconstructed; 
    int                              m_mcPdgCode;                  
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LArAnalysisParticle factory
 */
class LArAnalysisParticleFactory : public pandora::ObjectFactory<object_creation::ParticleFlowObject::Parameters, pandora::ParticleFlowObject>
{
public:
    /**
     *  @brief  Create new parameters instance on the heap (memory-management to be controlled by user)
     *
     *  @return the address of the new parameters instance
     */
    Parameters *NewParameters() const;

    /**
     *  @brief  Read any additional (derived class only) object parameters from file using the specified file reader
     *
     *  @param  parameters the parameters to pass in constructor
     *  @param  fileReader the file reader, used to extract any additional parameters from file
     */
    pandora::StatusCode Read(Parameters &parameters, pandora::FileReader &fileReader) const;

    /**
     *  @brief  Persist any additional (derived class only) object parameters using the specified file writer
     *
     *  @param  pObject the address of the object to persist
     *  @param  fileWriter the file writer
     */
    pandora::StatusCode Write(const pandora::ParticleFlowObject *const pObject, pandora::FileWriter &fileWriter) const;

    /**
     *  @brief  Create an object with the given parameters
     *
     *  @param  parameters the parameters to pass in constructor
     *  @param  pObject to receive the address of the object created
     */
    pandora::StatusCode Create(const object_creation::ParticleFlowObject::Parameters &parameters, const pandora::ParticleFlowObject *&pObject) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArAnalysisParticle::TYPE LArAnalysisParticle::Type() const noexcept
{
    return this->m_type;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArAnalysisParticle::TypeTree & LArAnalysisParticle::GetTypeTree() const noexcept
{
    return this->m_typeTree;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArAnalysisParticle::AnalysisEnergy() const noexcept
{
    return this->m_analysisEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArAnalysisParticle::EnergyFromCharge() const noexcept
{
    return this->m_energyFromCharge;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArAnalysisParticle::IsVertexFiducial() const noexcept
{
    return this->m_isVertexFiducial;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArAnalysisParticle::AreAllHitsFiducial() const noexcept
{
    return this->m_areAllHitsFiducial;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const CartesianVector & LArAnalysisParticle::VertexPosition() const noexcept
{
    return this->m_vertexPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const CartesianVector & LArAnalysisParticle::DirectionCosines() const noexcept
{
    return this->m_directionCosines;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const CartesianVector & LArAnalysisParticle::AnalysisMomentum() const noexcept
{
    return this->m_analysisMomentum;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned LArAnalysisParticle::NumberOf3dHits() const noexcept
{
    return this->m_numberOf3dHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned LArAnalysisParticle::NumberOfCollectionPlaneHits() const noexcept
{
    return this->m_numberOfCollectionPlaneHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArAnalysisParticle::IsShower() const noexcept
{
    return this->m_isShower;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned LArAnalysisParticle::NumberOfDownstreamParticles() const noexcept
{
    return this->m_numberOfDownstreamParticles;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArAnalysisParticle::HasMcInfo() const noexcept
{
    return this->m_hasMcInfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArAnalysisParticle::TYPE LArAnalysisParticle::McType() const
{
    ThrowIfNoMcInfo(); 
    return this->m_mcType;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArAnalysisParticle::TypeTree & LArAnalysisParticle::GetMcTypeTree() const
{
    ThrowIfNoMcInfo(); 
    return this->m_mcTypeTree;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArAnalysisParticle::McEnergy() const
{
    ThrowIfNoMcInfo(); 
    return this->m_mcEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const CartesianVector & LArAnalysisParticle::McMomentum() const
{
    ThrowIfNoMcInfo(); 
    return this->m_mcMomentum;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const CartesianVector & LArAnalysisParticle::McVertexPosition() const
{
    ThrowIfNoMcInfo(); 
    return this->m_mcVertexPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const CartesianVector & LArAnalysisParticle::McDirectionCosines() const
{
    ThrowIfNoMcInfo(); 
    return this->m_mcDirectionCosines;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArAnalysisParticle::McIsVertexFiducial() const
{
    ThrowIfNoMcInfo(); 
    return this->m_mcIsVertexFiducial;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArAnalysisParticle::McIsContained() const
{
    ThrowIfNoMcInfo(); 
    return this->m_mcIsContained;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArAnalysisParticle::McIsShower() const
{
    ThrowIfNoMcInfo(); 
    return this->m_mcIsShower;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArAnalysisParticle::McIsCorrectlyReconstructed() const
{
    ThrowIfNoMcInfo(); 
    return this->m_mcIsCorrectlyReconstructed;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int LArAnalysisParticle::McPdgCode() const
{
    ThrowIfNoMcInfo(); 
    return this->m_mcPdgCode;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::string LArAnalysisParticle::TypeTreeAsString(const TypeTree &typeTree)
{
    return TypeTreeAsStringImpl(typeTree, false);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArAnalysisParticle::TypeTree::TypeTree(const TYPE type, TypeTree::List daughterTypes) noexcept :
    m_type{type},
    m_daughterTypes{std::move_if_noexcept(daughterTypes)}
{
}
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArAnalysisParticle::TypeTree::TypeTree() noexcept :
    m_type{TYPE::UNKNOWN}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArAnalysisParticle::TypeTree::TypeTree(const TYPE type) noexcept :
    m_type{type},
    m_daughterTypes{}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
        
inline LArAnalysisParticle::TYPE LArAnalysisParticle::TypeTree::Type() const noexcept
{
    return this->m_type;
}
        
//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArAnalysisParticle::TypeTree::List & LArAnalysisParticle::TypeTree::Daughters() const noexcept
{
    return this->m_daughterTypes;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArAnalysisParticle::ThrowIfNoMcInfo() const
{
    if (!this->m_hasMcInfo)
    {
        std::cerr << "Could not get MC property of LArAnalysisParticle as no MC information was attached" << std::endl;
        throw STATUS_CODE_NOT_FOUND;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArAnalysisParticleFactory::Parameters *LArAnalysisParticleFactory::NewParameters() const
{
    return (new LArAnalysisParticleParameters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArAnalysisParticleFactory::Create(const Parameters &parameters, const pandora::ParticleFlowObject *&pObject) const
{
    const LArAnalysisParticleParameters &larPfoParameters(dynamic_cast<const LArAnalysisParticleParameters&>(parameters));
    pObject = new LArAnalysisParticle(larPfoParameters);

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArAnalysisParticleFactory::Read(Parameters&, pandora::FileReader&) const
{
    // TODO: Provide this functionality when necessary

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArAnalysisParticleFactory::Write(const pandora::ParticleFlowObject*, pandora::FileWriter&) const
{
    // TODO: Provide this functionality when necessary

    return pandora::STATUS_CODE_SUCCESS;
}
} // namespace lar_physics_content

#endif // #ifndef LAR_ANALYSIS_PARTICLE_H
