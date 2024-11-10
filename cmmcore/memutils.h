//
// Created by Andrew Astakhov on 07.11.24.
//

#ifndef CMMCORE_MEMUTILS_H
#define CMMCORE_MEMUTILS_H


#ifdef _WIN32
#include <malloc.h>
#define cmmalloc(size) _alloca(size)
#else
#include <alloca.h>
#define cmmalloc(size) alloca(size)
#endif

#endif // CMMCORE_MEMUTILS_H
