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

#include "PBDSolver/PBDSolver2D.h"
#include <Corrade/Utility/DebugStl.h>
#include <algorithm>
#include <random>
#include <chrono>



namespace Magnum { namespace Examples {

PBDSolver2D::PBDSolver2D(const Vector2& origin, Float cellSize, Int nI, Int nJ, SceneObjects* sceneObjs):
    _objects{sceneObjs},
    _particles{cellSize*0.5f},
    _cellSize{ cellSize },
    _origin{origin},
    _ni{nI}, _nj{nJ}

{
    if(nI < 1 || nJ < 1) {
        Fatal{} << "Invalid grid resolution";
    }
    if(!sceneObjs) {
        Fatal{} << "Invalid scene object";
    }

    /* Initialize data */
    generateParticles(_objects->emitterT0, 0);
}

void PBDSolver2D::emitParticles(const Vector2& pos, Float size, SDFObject::ObjectType type) {
    
    SDFObject *emitter = new SDFObject{ pos, size, type };
    SDFObject *boundary = new SDFObject{ _objects->boundary.center,
                                         _objects->boundary.radii,
                                         _objects->boundary.type };
    

    SDFObject emitterMinusBoundary{ emitter, boundary, SDFObject::ObjectType::Intersection };
    generateParticles(emitterMinusBoundary, 10);    
}

void PBDSolver2D::generateParticles(const SDFObject& sdfObj, Float initialVelocity_y) {

    /* Generate new particles */
    std::vector<Vector2> newParticles;
    for (Int i = 0; i < _ni; i++) {
        for (Int j = 0; j < _ni; j++) {

            const Vector2 gridPos{ i + 0.5f, j + 0.5f };            
            const Vector2 ppos = gridPos * _cellSize + _origin;
            if (sdfObj.signedDistance(ppos) < 0) {
                newParticles.push_back(ppos);
            }
        }
    }

    /* Insert into the system */
    _particles.addParticles(newParticles, initialVelocity_y);
}


void PBDSolver2D::findNeighbors(const std::vector<Vector2>& pos, std::vector<std::vector<int>> &nbors) const {

    
    KDTree<Vector2, 2> kd;
    kd.build(pos);

    TaskScheduler::forEach(pos.size(), [&](const std::size_t i) {

        std::vector<int> candidates = kd.radiusSearch(pos[i], _cellSize * 1.75);

        for (auto n = candidates.begin(); n != candidates.end(); ++n) {

            // eliminate self
            if (*n == i) continue;
            nbors[i].push_back(*n);            
        }
    });
    kd.clear();    
}


void PBDSolver2D::advanceFrame(Float frameDuration, Int nstep) {

    Float substep = frameDuration / static_cast<float>(nstep);

    const Vector2 gravityForce{ 0.0f, -gravity*substep };

    auto start = std::chrono::steady_clock::now();    
    std::vector<std::vector<int>> nbors(_particles.positions.size());
    findNeighbors(_particles.positions, nbors);
    auto stop = std::chrono::steady_clock::now();
    using FpMilliseconds = std::chrono::duration<float, std::chrono::milliseconds::period>;
    Float f_ms = FpMilliseconds(stop - start).count();
    
    //Debug{} << "Kd-tree took " << f_ms << " ms.\n";

    // make friction independent of number of constraint iterations
    Float boundaryFriction = Math::pow(1.0f - _boundaryFriction, 1 / static_cast<Float>(_numConstraintIteration));
    Float friction = Math::pow(1.0f - _friction, 1 / static_cast<Float>(_numConstraintIteration));

    // each frame is decomposed into a number of substeps
    start = std::chrono::steady_clock::now();
    std::vector<Vector2> new_p(_particles.positions.size());
    for (auto step = 0; step < nstep; step++) {
        

        TaskScheduler::forEach(_particles.positions.size(), [&](const std::size_t i) {
            // new estimated velocities by integrating forces over substep
            _particles.velocities[i] += gravityForce;

            // new estimated positions by using estimated velocities at the end of the time step                
            new_p[i] = _particles.positions[i] + _particles.velocities[i] * substep;
            }
        );
        
        // iteratively constrain estimated positions
        for (auto i = 0; i < _numConstraintIteration; i++) {
            solveBoundaryConstraints(new_p, boundaryFriction);
            solveParticleConstraints(new_p, nbors, friction);
        }

        // true up velocities based on previous and new position
        TaskScheduler::forEach(_particles.positions.size(), [&](const std::size_t i) {
                _particles.velocities[i] = (new_p[i] - _particles.positions[i]) / substep;
            }
        );
        _particles.positions = new_p;
    }
    stop = std::chrono::steady_clock::now();    
    f_ms = FpMilliseconds(stop - start).count();

    //Debug{} << "One Frame took " << f_ms << " ms.\n";
}


void PBDSolver2D::solveParticleConstraints(std::vector<Vector2>& new_p, const std::vector<std::vector<int>>& nbors, Float friction) const {

    // record collisions with other particles (as indices)
    std::vector<std::vector<int>> collision(new_p.size());
 
    
    // loop over each particle and test collision with other particles
    //
    // While it is tempting to multithread the following loop, it leads to instability 
    // as the symmetry of corrective displacements is no longer guaranteed.
    // As a result, oscillations can occur for neighboring particles
    //
    // By looping iteratively without multithreading, we ensure that when two 
    // particles collide, they are each moved by the same amount in opposite directions.    
    for (unsigned int idx = 0; idx < new_p.size(); idx++) {

        // inter particles collisions
        // loop over neighbors
        for (auto nbor = nbors[idx].begin(); nbor != nbors[idx].end(); ++nbor) {

            // use symmetry for speed and stability
            if (*nbor < idx) continue;

            Vector2 dir = new_p[idx] - new_p[*nbor];
            Float dist = dir.length() - _cellSize;
            if (dist < 0) {
                collision[idx].push_back(*nbor);
                collision[*nbor].push_back(idx);
                Float strength = -0.25f * dist; // -0.25 = 0.5 (constrait weight) * 0.5 (half on each particle) * -1 (dist<0)
                dir = strength * dir.normalized();
                new_p[idx] += dir;
                new_p[*nbor] -= dir;
            }
        }
    }

    // loop over each particle and add friction in the tangential direction
    TaskScheduler::forEach(new_p.size(), [&](const std::size_t idx) {
        //int idx = indices[i];
        Vector2& p = new_p[idx];

        for (auto col = collision[idx].begin(); col != collision[idx].end(); ++col) {
            Vector2 N = (p - new_p[*col]).normalized();            
            Vector2 delta = p - _particles.positions[idx];
            Vector2 normalDelta = dot(delta, N) * N;
            Vector2 tangentialDelta = delta - normalDelta;
         
            p = _particles.positions[idx] + normalDelta + friction * tangentialDelta;
        }
    }
    );
}

void PBDSolver2D::solveBoundaryConstraints(std::vector<Vector2>& new_p, Float friction) const {

    
    // loop over each particle and test collision with boundary
    TaskScheduler::forEach(new_p.size(), [&](const std::size_t idx) {
        
        Vector2& p = new_p[idx];
        
            // record collision with boundary
        bool boundaryCollision = false;
        Float dist;

        // collision with the boundary
        if (_objects->boundary.type == SDFObject::ObjectType::Circle) {
            dist = _objects->boundary.signedDistance(p);
            if (dist < 0) {
                boundaryCollision = true;
                Vector2 dir = (p - _objects->boundary.center).normalized();
                p += 0.5f * dist * dir;
            }
        }
        else {
            Float s = std::sin(-_orient);
            Float c = std::cos(-_orient);
            Vector2 rot_p{ c * p[0] - s * p[1], s * p[0] + c * p[1] };
            Vector2 d{ 0,0 };

            // collision with box, solve along each axis independently
            for (unsigned int axis = 0; axis < 2; axis++) {                
                dist = Math::abs(rot_p[axis] - _objects->boundary.center[axis]) - _objects->boundary.radii[axis];
                if (dist < 0) continue;
                if (rot_p[axis] > _objects->boundary.center[axis]) dist = -dist;
                
                d[axis] += 0.5f * dist;
                boundaryCollision = true;
            }
            s = std::sin(_orient);
            c = std::cos(_orient);            
            p[0] += c * d[0] - s * d[1];
            p[1] += s * d[0] + c * d[1];
            
        }
        // add friction in the tangential direction
        if (boundaryCollision) {

            Vector2 N = (p - _objects->boundary.center).normalized();
            Vector2 delta = p - _particles.positions[idx];
            Vector2 normalDelta = dot(delta, N) * N;
            Vector2 tangentialDelta = delta - normalDelta;
            
            p = _particles.positions[idx] + normalDelta + friction * tangentialDelta;
        }
    });
}

}}
