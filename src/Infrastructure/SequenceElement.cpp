#include "TrackingPipeline/Infrastructure/SequenceElement.hpp"

void SequenceElement::registerWriteHandle(const DataHandleBase& handle) {
    m_writeHandles.push_back(&handle);
}

void SequenceElement::registerReadHandle(const DataHandleBase& handle) {
    m_readHandles.push_back(&handle);
}

const std::vector<const DataHandleBase*>& SequenceElement::writeHandles()
    const {
        return m_writeHandles;
}

const std::vector<const DataHandleBase*>& SequenceElement::readHandles() const {
    return m_readHandles;
}
