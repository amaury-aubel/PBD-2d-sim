/*
    This file is part of Magnum.

    Original authors — credit is appreciated but not required:

        2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020 —
            Vladimír Vondruš <mosra@centrum.cz>
        2019 — Nghia Truong <nghiatruong.vn@gmail.com>
        2020 — Amaury Aubel <Amaury.Aubel@gmail.com>

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
#include <Corrade/Utility/DebugStl.h>
#include <Corrade/Utility/StlMath.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/GL/PixelFormat.h>
#include <Magnum/GL/Context.h>
#include <Magnum/GL/Version.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Math/FunctionsBatch.h>

#include <Magnum/MeshTools/Compile.h>
#include <Magnum/Platform/Sdl2Application.h>
#include <Magnum/Primitives/Circle.h>
#include <Magnum/Primitives/Plane.h>
#include <Magnum/SceneGraph/Camera.h>
#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>
#include <Magnum/SceneGraph/Scene.h>
#include <Magnum/Timeline.h>
#include <Magnum/Trade/MeshData.h>

#include "DrawableObjects/ParticleGroup2D.h"
#include "DrawableObjects/WireframeObject2D.h"
#include "PBDSolver/PBDSolver2D.h"
#include <chrono>
#include <random>

namespace Magnum { namespace Examples {

class PBDSimulation2DExample: public Platform::Application {
    public:
        explicit PBDSimulation2DExample(const Arguments& arguments);

    protected:
        void viewportEvent(ViewportEvent& event) override;
        void keyPressEvent(KeyEvent& event) override;
        void keyReleaseEvent(KeyEvent& event) override;
        void mousePressEvent(MouseEvent& event) override;
        void drawEvent() override;

        
        /* Simulation helper functions */
        void setupSolver(SDFObject::ObjectType boundaryType, bool firstTime = false);
        void resetSimulation();
        void changeResolution();
        void showInfo();
        void toggleBoundary();        

        /* Window control */        
        Vector2 windowPos2WorldPos(const Vector2i& winPos);

        /* Scene and drawable group must be constructed before camera and other
           drawble objects */
        Containers::Pointer<Scene2D> _scene;
        Containers::Pointer<SceneGraph::DrawableGroup2D> _drawableGroup;

        /* Camera helpers */
        Containers::Pointer<Object2D> _objCamera;
        Containers::Pointer<SceneGraph::Camera2D> _camera;

        /* PBD simulation system */
        Containers::Pointer<PBDSolver2D> _pbdSolver;
        Containers::Pointer<ParticleGroup2D> _drawableParticles;
        Containers::Pointer<WireframeObject2D> _drawableBoundary;
        Float _speed = 3.0f;
        Float _evolvedTime = 0.0f;
        bool _pausedSimulation = false;

        /* Mouse-Particles interaction */
        Containers::Pointer<WireframeObject2D> _drawablePointer;
        Timeline _timeline;
        Vector2 _lastMousePressedWorldPos;
        Float _mouseInteractionRadius = 5.0f;
        Float _orient = 0.0f;
        Float _orientInc = 0.0f;
};

namespace {


Vector2i NumGridCells{100};       /* number of cells */
Float GridCellLength = 1.0f;      /* length of 1 grid cell */
constexpr Vector2 GridStart{-50.0f, -50.0f}; /* lower corner of the grid */
constexpr Int RadiusCircleBoundary = 45; /* radius of the boundary circle */

/* Viewport will display this window */
constexpr Float ProjectionScale = 1.05f;
const Vector2i DomainDisplaySize = NumGridCells*GridCellLength*ProjectionScale;

Vector2 gridCenter() {
    return Vector2{NumGridCells}*GridCellLength*0.5f + GridStart;
}

}

using namespace Math::Literals;

PBDSimulation2DExample::PBDSimulation2DExample(const Arguments& arguments): Platform::Application{arguments, NoCreate} {

    _pbdSolver = nullptr;

    /* Setup window */
    {
        const Vector2 dpiScaling = this->dpiScaling({});
        Configuration conf;
        conf.setTitle("2D PBD Simulation - Multithreaded")
            .setSize(conf.size(), dpiScaling)
            .setWindowFlags(Configuration::WindowFlag::Resizable);
        GLConfiguration glConf;
        glConf.setSampleCount(dpiScaling.max() < 2.0f ? 8 : 2);
        if(!tryCreate(conf, glConf)) {
            create(conf, glConf.setSampleCount(0));
        }
    }

    /* Setup scene objects and camera */
    {
        /* Setup scene objects */
        _scene.emplace();
        _drawableGroup.emplace();

        /* Configure camera */
        _objCamera.emplace(_scene.get());
        _objCamera->setTransformation(Matrix3::translation(gridCenter()));

        _camera.emplace(*_objCamera);
        _camera->setAspectRatioPolicy(SceneGraph::AspectRatioPolicy::Extend)
            .setProjectionMatrix(Matrix3::projection(Vector2{DomainDisplaySize}))
            .setViewport(GL::defaultFramebuffer.viewport().size());
    }

    /* Setup solver */
    {
        setupSolver(SDFObject::ObjectType::Box, true);
    }

    /* Enable depth test, render particles as sprites */
    GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
    GL::Renderer::enable(GL::Renderer::Feature::ProgramPointSize);

    /* Start the timer, loop at 60 Hz max */
    setSwapInterval(1);
    setMinimalLoopPeriod(16);
    _timeline.start();

    showInfo();
}

void PBDSimulation2DExample::setupSolver(SDFObject::ObjectType boundaryType, bool firstTime) {

    SceneObjects* sceneObjs = new SceneObjects;
    sceneObjs->emitterT0 = SDFObject{ gridCenter() + Vector2(10.0f, -5.0f), 12.0f, SDFObject::ObjectType::Circle };

    if (boundaryType == SDFObject::ObjectType::Circle) 
        sceneObjs->boundary = SDFObject{ gridCenter(), Float(RadiusCircleBoundary), SDFObject::ObjectType::Circle, false };
    else sceneObjs->boundary = SDFObject{ gridCenter(), Vector2(RadiusCircleBoundary,RadiusCircleBoundary), SDFObject::ObjectType::Box, false };

    if (firstTime) {
        _pbdSolver.emplace(GridStart, GridCellLength, NumGridCells.x(), NumGridCells.y(), sceneObjs);
     
        /* Drawable particles */
        _drawableParticles.emplace(_pbdSolver->particlePositions(),
            _pbdSolver->particleRadius() * 1.25f);
        //_drawableParticles->setColor(0x55c8f5_rgbf);
    }
    else _pbdSolver->updateBoundary(sceneObjs);
        
    /* Drawable boundary*/
    if (boundaryType == SDFObject::ObjectType::Circle) 
        _drawableBoundary.emplace(_scene.get(), _drawableGroup.get(),
            MeshTools::compile(Primitives::circle2DWireframe(128)));
    else
        _drawableBoundary.emplace(_scene.get(), _drawableGroup.get(),
            MeshTools::compile(Primitives::planeWireframe()));

    Matrix3 mat = Matrix3::rotation( Rad(_orient) );
    mat[0][0] *= RadiusCircleBoundary + _pbdSolver->particleRadius();
    mat[0][1] *= RadiusCircleBoundary + _pbdSolver->particleRadius();
    mat[1][0] *= RadiusCircleBoundary + _pbdSolver->particleRadius();
    mat[1][1] *= RadiusCircleBoundary + _pbdSolver->particleRadius();
    _drawableBoundary->setTransformation(mat);
    _drawableBoundary->setColor(0xffffff_rgbf);
}

void PBDSimulation2DExample::drawEvent() {
    GL::defaultFramebuffer.clear(GL::FramebufferClear::Color | GL::FramebufferClear::Depth);    

    /* Draw objects */
    {
        /* orient boundary */
        Matrix3 mat = Matrix3::rotation(Rad(_orient));
        mat[0][0] *= RadiusCircleBoundary + _pbdSolver->particleRadius();
        mat[0][1] *= RadiusCircleBoundary + _pbdSolver->particleRadius();
        mat[1][0] *= RadiusCircleBoundary + _pbdSolver->particleRadius();
        mat[1][1] *= RadiusCircleBoundary + _pbdSolver->particleRadius();
        _drawableBoundary->setTransformation(mat);

        /* Trigger drawable object to update the particles to the GPU */
        _drawableParticles->setDirty();
        _drawableParticles->draw(_camera, GL::defaultFramebuffer.viewport().size().y(), DomainDisplaySize.y());

        /* Draw other objects (boundary mesh, pointer mesh) */
        _camera->draw(*_drawableGroup);
    }

    if(!_pausedSimulation) {
        _orient += _orientInc;
        _pbdSolver->orientBoundary(_orient);

        constexpr Float frameTime = 1.0f/60.0f;
        /* pause for a while before starting simulation */
        if(_evolvedTime > 0.65f) _pbdSolver->advanceFrame(frameTime*_speed);
        _evolvedTime += frameTime;
    }

    swapBuffers();

    /* Run next frame immediately */
    redraw();
}

void PBDSimulation2DExample::viewportEvent(ViewportEvent& event) {
    /* Resize the main framebuffer */
    GL::defaultFramebuffer.setViewport({{}, event.framebufferSize()});
  
    /* Recompute the camera's projection matrix */
    _camera->setViewport(event.framebufferSize());
}

void PBDSimulation2DExample::keyPressEvent(KeyEvent& event) {
    // obtain a seed from the system clock:
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed);
    std::uniform_real_distribution<float> dis(-30.0, 30.0);

    switch(event.key()) {
        case KeyEvent::Key::F:
            _pbdSolver->multFriction(1.5f);
            break;
        case KeyEvent::Key::D:
            _pbdSolver->multFriction(0.75f);
            break;
        case KeyEvent::Key::One:
            _pbdSolver->setIterations(1);
            break;
        case KeyEvent::Key::Two:
            _pbdSolver->setIterations(2);
            break;
        case KeyEvent::Key::Three:
            _pbdSolver->setIterations(3);
            break;
        case KeyEvent::Key::Four:
            _pbdSolver->setIterations(4);
            break;
        case KeyEvent::Key::Five:
            _pbdSolver->setIterations(5);
            break;
        case KeyEvent::Key::Six:
            _pbdSolver->setIterations(6);
            break;
        case KeyEvent::Key::Seven:
            _pbdSolver->setIterations(7);
            break;
        case KeyEvent::Key::Eight:
            _pbdSolver->setIterations(8);
            break;
        case KeyEvent::Key::Nine:
            _pbdSolver->setIterations(9);
            break;
        case KeyEvent::Key::Left:
            _orientInc += 0.002f;
            break;
        case KeyEvent::Key::Right:
            _orientInc -= 0.002f;
            break;
        case KeyEvent::Key::B:
            toggleBoundary();
            break;
        case KeyEvent::Key::W:            
            _pbdSolver->emitParticles(gridCenter() + Vector2(Float(dis(gen)), 20.0f), 15.0f, SDFObject::ObjectType::Circle);
            break;
        case KeyEvent::Key::E:
            _pbdSolver->emitParticles(gridCenter() + Vector2(15.0f, 20.0f), 15.0f, SDFObject::ObjectType::Circle);
            break;
        case KeyEvent::Key::H:
        case KeyEvent::Key::I:
            showInfo();
            break;
        case KeyEvent::Key::R:
            resetSimulation();
            break;
        case KeyEvent::Key::Y:
            changeResolution();
            break;
        case KeyEvent::Key::Space:
            _pausedSimulation ^= true;
            break;        
        default: event.setAccepted(true);
    }
}
void PBDSimulation2DExample::toggleBoundary() {

    resetSimulation();
    if (_pbdSolver->getBoundaryType() == SDFObject::ObjectType::Circle)
        setupSolver(SDFObject::ObjectType::Box);
    else 
        setupSolver(SDFObject::ObjectType::Circle);
}


void PBDSimulation2DExample::keyReleaseEvent(KeyEvent& event) {
        event.setAccepted(true);
    
}

void PBDSimulation2DExample::mousePressEvent(MouseEvent& event) {
    
    _lastMousePressedWorldPos = windowPos2WorldPos(event.position());
    _timeline.nextFrame();
    _pbdSolver->emitParticles(gridCenter() + _lastMousePressedWorldPos, _mouseInteractionRadius, SDFObject::ObjectType::Circle);
    event.setAccepted();
}


void PBDSimulation2DExample::showInfo() {
    Debug info{};
    info << "\n";
    info << _pbdSolver->numParticles() << " particles simulated with " << _pbdSolver->numConstraintIterations() << " constraint iterations\n";
    info << "Inter-particles Friction : "<< _pbdSolver->getFriction() << "\n";
    info << "Boundary Friction : " << _pbdSolver->getBoundaryFriction() << "\n";
    info << "\n\n";
    info << "R: reset simulation\n";
    info << "Y: change resolution\n";
    info << "B: change boundary\n";
    info << "E: emit into simulation\n";
    info << "W: emit randomly into simulation\n";
    info << "left/right: rotate container (box container)\n";
    info << "left mouse button: emit into simulation\n";
    info << "H/I: display info\n";
    info << "F: Increase friction\n";
    info << "D: Decrease friction\n";
    info << "1-9: Set number of constraint iterations\n";
}

void PBDSimulation2DExample::resetSimulation() {
    _pbdSolver->reset();
    _pausedSimulation = false;
    _evolvedTime = 0.0f;
    _orient = 0.0f;
    _orientInc = 0.0f;
    showInfo();
}


void PBDSimulation2DExample::changeResolution() {

    if (abs(GridCellLength - 1.0f) < 1e-3) {
        GridCellLength = 100.f / 160.f;
        NumGridCells = Vector2i{ 160 };
    }
    else {
        GridCellLength = 1.0f;
        NumGridCells = Vector2i{ 100 };
    }
    resetSimulation();
    setupSolver(_pbdSolver->getBoundaryType(), true);
}


Vector2 PBDSimulation2DExample::windowPos2WorldPos(const Vector2i& windowPosition) {
    /* First scale the position from being relative to window size to being
       relative to framebuffer size as those two can be different on HiDPI
       systems */
    const Vector2i position = windowPosition*Vector2{framebufferSize()}/Vector2{windowSize()};

    /* Compute inverted model view projection matrix */
    const Matrix3 invViewProjMat = (_camera->projectionMatrix()*_camera->cameraMatrix()).inverted();

    /* Compute the world coordinate from window coordinate */
    const Vector2i flippedPos = Vector2i(position.x(), framebufferSize().y() - position.y());
    const Vector2 ndcPos = Vector2(flippedPos) / Vector2(framebufferSize())*Vector2{2.0f} - Vector2{1.0f};
    return invViewProjMat.transformPoint(ndcPos);
}

}}

MAGNUM_APPLICATION_MAIN(Magnum::Examples::PBDSimulation2DExample)
