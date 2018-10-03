/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#include <core/integrator.h>

#include "tiny_obj_loader.h"

TR_NAMESPACE_BEGIN

Integrator::Integrator(const Scene& scene) : scene(scene) { }

bool Integrator::init() {
    rgb = std::unique_ptr<RenderBuffer>(new RenderBuffer(scene.config.width, scene.config.height));
    rgb->clear();
    return true;
}

void Integrator::cleanUp() {
    save();
}

bool Integrator::save() {
    fs::path p = scene.config.tomlFile;
    saveEXR(rgb->data, p.replace_extension("exr").string(), scene.config.width, scene.config.height);
    return true;
}

const Emitter& Integrator::getEmitterByID(const int emitterID) const {
    return scene.emitters[emitterID];
}
const BSDF* Integrator::getBSDF(const SurfaceInteraction& hit) const {
    assert(hit.shapeID < scene.worldData.shapes.size());
    const tinyobj::shape_t& shape = scene.worldData.shapes[hit.shapeID];
    return scene.bsdfs[shape.mesh.material_ids[hit.primID]].get();
}

size_t Integrator::selectEmitter(float sample, float& pdf) const {
    size_t id = size_t(sample * scene.emitters.size());
    id = min(id, scene.emitters.size() - 1); //todo @nico : how can this happen ? (sample ==1)
    pdf = 1.f / scene.emitters.size();
    return id;
}

size_t Integrator::getEmitterIDByShapeID(size_t shapeID) const {
    auto it = std::find_if(scene.emitters.begin(), scene.emitters.end(),
                           [shapeID](const Emitter& s) { return s.shapeID == shapeID; });
    assert(it != scene.emitters.end());
    return (size_t) std::distance(scene.emitters.begin(), it);
}

float Integrator::getEmitterPdf(const Emitter& emitter) const {
    return 1.f / scene.emitters.size();
}

void Integrator::sampleEmitterDirection(Sampler& sampler,
                                        const Emitter& emitter,
                                        const v3f& n,
                                        v3f& d,
                                        float& pdf) const {
}

void Integrator::sampleEmitterPosition(Sampler& sampler, const Emitter& emitter, v3f& n, v3f& pos, float& pdf) const {
}

TR_NAMESPACE_END