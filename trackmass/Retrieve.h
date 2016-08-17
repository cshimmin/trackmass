#include "xAODRootAccess/TEvent.h"
#include "AthContainers/DataVector.h"

#include <algorithm>
#include <vector>

template <class T>
const DataVector<T>* Retrieve(xAOD::TEvent* evt, const std::string& container_name) {
  const DataVector<T>* container(nullptr);
  evt->retrieve(container, container_name);
  return container;
}

template <class T>
std::vector<const T* > CopyRetrieve(xAOD::TEvent* evt, const std::string& container_name) {
  const DataVector<T>* container(nullptr);
  evt->retrieve(container, container_name);
  std::vector<const T* > copied(container->size());
  std::copy(container->begin(), container->end(), copied.begin());
  return copied;
}
