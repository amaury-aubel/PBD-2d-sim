#ifndef Magnum_Examples_PBDSimulation2D_Solver2D_h
#define Magnum_Examples_PBDSimulation2D_Solver2D_h
/*
    This file is part of Magnum.

    Original authors — credit is appreciated but not required:

        2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020 —
            Vladimír Vondruš <mosra@centrum.cz>
        2019 — Nghia Truong <nghiatruong.vn@gmail.com>

    This is free and unencumbered software released into the public domain.

    Anyone is free to copy, modify, publish, use, compile, sell, or distribute
    this software, either in source code form or as a compiled binary, for any
    purpose, commercial or non-commercial, and by any means.

    In jurisdictions that recognize copyright laws, the author or authors of
    this software dedicate any and all copyright interest in the software to
    the public domain. We make this dedication for the benefit of the public
    at large and to the detriment of our heirs and successors. We intend this
    dedication to be an overt act of relinquishment in perpetuity of all
    present and future rights to this software under copyright law.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
    IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include <Corrade/Containers/Pointer.h>
#include <tbb/parallel_for.h>
#include "PBDSolver/SolverData.h"


namespace Magnum { namespace Examples {

namespace TaskScheduler {

    template<class IndexType, class Function> void forEach(IndexType endIdx, Function&& func) {
        tbb::parallel_for(tbb::blocked_range<IndexType>(IndexType(0), endIdx),
            [&](const tbb::blocked_range<IndexType>& r) {
                for (IndexType i = r.begin(), iEnd = r.end(); i < iEnd; ++i) {
                    func(i);
                }
            });
    }
}

/* 2D PBD solver */
class PBDSolver2D {
public:
    explicit PBDSolver2D(const Vector2& origin, Float cellSize, Int Ni, Int Nj, SceneObjects* sceneObjs);

    void emitParticles(const Vector2& pos, Float size, SDFObject::ObjectType type);

    /* Manipulation */
    void reset() {
        _particles.reset();
        _particles.addParticles(_particles.positionsT0, 0);
    }
    
    void setIterations(unsigned int num) { _numConstraintIteration = num; }
    void advanceFrame(Float frameDuration, Int nstep=3);

    /* Properties */
    UnsignedInt numParticles() const { return _particles.size(); }
    auto numConstraintIterations() const { return _numConstraintIteration; }
    Float particleRadius() const { return _particles.particleRadius; }

    const std::vector<Vector2>& particlePositions() const {
        return _particles.positions;
    }

    SDFObject::ObjectType getBoundaryType() { return _objects->boundary.type; }
    void updateBoundary(SceneObjects* sceneObjs) { delete _objects;  _objects = sceneObjs; }
    void orientBoundary(Float orient) { _objects->boundary.setOrient(orient); _orient = orient;  }

    // friction
    void multFriction(Float mult) {
        _friction = Math::clamp(_friction*mult, 1e-4f, 0.5f); 
        _boundaryFriction = Math::clamp(_boundaryFriction * mult, 1e-4f, 0.5f);
    }
    Float getFriction() const { return _friction; }
    Float getBoundaryFriction() const { return _boundaryFriction; }

private:
    /* Initialization */
    void generateParticles(const SDFObject& sdfObj, Float initialVelocity_y);

    /* Simulation */
    void findNeighbors(const std::vector<Vector2>& pos, std::vector<std::vector<int>>& nbors) const;
    void solveBoundaryConstraints(std::vector<Vector2> &new_pos) const;
    void solveParticleConstraints(std::vector<Vector2> &new_pos, const std::vector<std::vector<int>>& nbors) const;

    SceneObjects *_objects;
    ParticleData _particles;
    Vector2 _origin;
    Float _cellSize;
    Int _ni, _nj;
    Float gravity = 9.81f;
    Float _boundaryFriction = 0.02f;
    Float _friction = 0.04f;
    Float _orient = 0;
    unsigned int _numConstraintIteration = 2;
};

}}

#endif
