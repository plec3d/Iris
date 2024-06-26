///|/ Copyright (c) Prusa Research 2022 Tomáš Mészáros @tamasmeszaros
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
// Original implementation of STEP format import created by Bambulab.
// https://github.com/bambulab/BambuStudio
// Forked off commit 1555904, modified PLEC³ᵈ.

#ifndef slic3r_Format_STEP_hpp_
#define slic3r_Format_STEP_hpp_

namespace Slic3r {

class Model;

//typedef std::function<void(int load_stage, int current, int total, bool& cancel)> ImportStepProgressFn;

// Load a step file into a provided model.
extern bool load_step(const char *path_str, Model *model /*LMBBS:, ImportStepProgressFn proFn = nullptr*/);

}; // namespace Slic3r

#endif /* slic3r_Format_STEP_hpp_ */
