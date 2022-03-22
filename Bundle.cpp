#include "Bundle.h"

Bundle::Bundle(
               const Object* const hostObject,
               unordered_map<int64, Object*>& mobileObjects)
{
    for (int i = 0; i < POINTS; ++i) {
        if (hostObject->getNumber(i) != 0) {
            lines[i] = new OneLine(hostObject, i, mobileObjects);
        }
        else {
            lines[i] = nullptr;
        }
    }
}

Bundle::~Bundle()
{
    for (int i = 0; i < POINTS; ++i) {
        if (lines[i] != nullptr) {
            delete lines[i];
        }
    }
}
