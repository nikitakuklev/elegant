#ifndef DL_PATCH_MANUAL_H
#define DL_PATCH_MANUAL_H


#ifndef TURBO_MODE
#define TURBO_MODE 0
#endif


// note that even at level 0, some code has been moved around to headers
// 0 - least possible changes
// 1 - csbend field caching
// 2 - layout changes/static fields
// 3 - static declarations


// ANY SETTINGS >= 10 affect output (as does USE_OLD_IPOW in mdb.h)

// 10 - csbendfield
// 12 - expandPower/apply_canonical_multipole_kicks/etc.
// 13 - various small equation tweaks
// 14 - matrix memory layout overhaul (why output changes is unclear)

#endif //DL_PATCH_MANUAL_H
