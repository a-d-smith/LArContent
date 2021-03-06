/**
 *  @file   larpandoracontent/LArStitching/MultiPandoraApiImpl.cc
 * 
 *  @brief  Implementation of the MultiPandoraApiImpl class.
 * 
 *  $Log: $
 */

#include "Pandora/Pandora.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArStitching/MultiPandoraApiImpl.h"

const PandoraInstanceMap &MultiPandoraApiImpl::GetPandoraInstanceMap() const
{
    return m_pandoraInstanceMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const PandoraInstanceList &MultiPandoraApiImpl::GetDaughterPandoraInstanceList(const pandora::Pandora *const pPrimaryPandora) const
{
    PandoraInstanceMap::const_iterator iter = m_pandoraInstanceMap.find(pPrimaryPandora);

    if (m_pandoraInstanceMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::Pandora *MultiPandoraApiImpl::GetDaughterPandoraInstance(const pandora::Pandora *const pPrimaryPandora, const int idNumber) const
{
    PrimaryToVolumeIdMap::const_iterator iter = m_primaryToVolumeIdMap.find(pPrimaryPandora);

    if (m_primaryToVolumeIdMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    VolumeIdMap::const_iterator idIter = iter->second.find(idNumber);

    if (iter->second.end() == idIter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return idIter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::Pandora *MultiPandoraApiImpl::GetPrimaryPandoraInstance(const pandora::Pandora *const pDaughterPandora) const
{
    PandoraRelationMap::const_iterator iter = m_pandoraRelationMap.find(pDaughterPandora);

    if (m_pandoraRelationMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const VolumeInfo &MultiPandoraApiImpl::GetVolumeInfo(const pandora::Pandora *const pPandora) const
{
    VolumeInfoMap::const_iterator iter = m_volumeInfoMap.find(pPandora);

    if (m_volumeInfoMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return *(iter->second);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MultiPandoraApiImpl::SetVolumeInfo(const pandora::Pandora *const pPandora, VolumeInfo *const pVolumeInfo)
{
    if (m_volumeInfoMap.count(pPandora))
        throw pandora::StatusCodeException(pandora::STATUS_CODE_ALREADY_PRESENT);

    if (!m_volumeInfoMap.insert(VolumeInfoMap::value_type(pPandora, pVolumeInfo)).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    const pandora::Pandora *const pPrimaryPandora((m_pandoraInstanceMap.count(pPandora)) ? pPandora : this->GetPrimaryPandoraInstance(pPandora));
    PrimaryToVolumeIdMap::iterator iter = m_primaryToVolumeIdMap.find(pPrimaryPandora);

    if (m_primaryToVolumeIdMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    if (!iter->second.insert(VolumeIdMap::value_type(pVolumeInfo->GetIdNumber(), pPandora)).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MultiPandoraApiImpl::SetParticleX0(const pandora::Pandora *const pPandora, const pandora::ParticleFlowObject *const pPfo, const float x0)
{
    VolumeInfoMap::iterator iter = m_volumeInfoMap.find(pPandora);

    if (m_volumeInfoMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    iter->second->SetParticleX0(pPfo, x0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MultiPandoraApiImpl::ClearParticleX0Map(const pandora::Pandora *const pPandora)
{
    VolumeInfoMap::iterator iter = m_volumeInfoMap.find(pPandora);

    if (m_volumeInfoMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    iter->second->ClearParticleX0Map();
}

//------------------------------------------------------------------------------------------------------------------------------------------

MultiPandoraApiImpl::MultiPandoraApiImpl()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

MultiPandoraApiImpl::~MultiPandoraApiImpl()
{
    PandoraInstanceMap pandoraInstanceMap(m_pandoraInstanceMap);

    for (const auto &mapElement : pandoraInstanceMap)
        this->DeletePandoraInstances(mapElement.first);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MultiPandoraApiImpl::AddPrimaryPandoraInstance(const pandora::Pandora *const pPrimaryPandora)
{
    if (m_pandoraInstanceMap.count(pPrimaryPandora))
        throw pandora::StatusCodeException(pandora::STATUS_CODE_ALREADY_PRESENT);

    if (!m_pandoraInstanceMap.insert(PandoraInstanceMap::value_type(pPrimaryPandora, PandoraInstanceList())).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
        
    if (!m_primaryToVolumeIdMap.insert(PrimaryToVolumeIdMap::value_type(pPrimaryPandora, VolumeIdMap())).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MultiPandoraApiImpl::AddDaughterPandoraInstance(const pandora::Pandora *const pPrimaryPandora, const pandora::Pandora *const pDaughterPandora)
{
    PandoraInstanceMap::iterator iter = m_pandoraInstanceMap.find(pPrimaryPandora);

    if (m_pandoraInstanceMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    iter->second.push_back(pDaughterPandora);

    if (!m_pandoraRelationMap.insert(PandoraRelationMap::value_type(pDaughterPandora, pPrimaryPandora)).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MultiPandoraApiImpl::DeletePandoraInstances(const pandora::Pandora *const pPrimaryPandora)
{
    PandoraInstanceList pandoraInstanceList(this->GetDaughterPandoraInstanceList(pPrimaryPandora));
    pandoraInstanceList.push_back(pPrimaryPandora);
    m_pandoraInstanceMap.erase(pPrimaryPandora);
    m_primaryToVolumeIdMap.erase(pPrimaryPandora);
    
    for (const pandora::Pandora *const pPandora : pandoraInstanceList)
    {
        VolumeInfoMap::iterator volumeInfoIter = m_volumeInfoMap.find(pPandora);

        if (m_volumeInfoMap.end() != volumeInfoIter)
        {
            delete volumeInfoIter->second;
            m_volumeInfoMap.erase(volumeInfoIter);
        }

        m_pandoraRelationMap.erase(pPandora);
        delete pPandora;
    }
}
