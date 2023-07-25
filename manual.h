#ifndef DL_PATCH_MANUAL_H
#define DL_PATCH_MANUAL_H


#ifndef TURBO_MODE
#define TURBO_MODE 10
#endif


// note that even at level 0, some code has been moved around to headers
// 0 - least possible changes
// 1 - csbend field caching
// 2 - layout changes/static fields
// 3 - static declarations
// 4 - more statics
// 5 - csbendfield
// 6 - matrix memory layout overhaul (why output changes is unclear)
// 7 - static momenta conversions
// 8 - extract non-distribution radiation kick
// 9 - try to fix alignment of rootname copy, still need to reduce amount of copying or buffer size

// ANY SETTINGS >= 10 affect output (as does USE_OLD_IPOW in mdb.h)

// 11 - csbendfield reorder calculations
// 13 - various small equation tweaks



// >= 20 settings are still work in progress, do not use
// 22 - proper caching of multipole expansion orders/apply_canonical_multipole_kicks/etc.

#endif //DL_PATCH_MANUAL_H
