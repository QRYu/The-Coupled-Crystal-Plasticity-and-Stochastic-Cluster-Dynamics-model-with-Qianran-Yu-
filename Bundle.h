// Bundle.h -- struct Bundle
#pragma once
#include"OneLine.h"
struct Bundle {
    OneLine* lines[POINTS];  /* pointers that point to line */
    Bundle(const Object* const, unordered_map<int64, Object*>&);/* constructor */
    ~Bundle();               /* destructor */
};
