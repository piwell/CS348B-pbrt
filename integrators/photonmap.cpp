
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// integrators/photonmap.cpp*
#include "stdafx.h"
#include "integrators/photonmap.h"
#include "scene.h"
#include "montecarlo.h"
#include "sampler.h"
#include "progressreporter.h"
#include "intersection.h"
#include "paramset.h"
#include "camera.h"
#include "spectrum.h"
#include "photonshooter.h"


// PhotonIntegrator Local Declarations
inline float kernel(const Photon *photon, const Point &p, float maxDist2);
static Spectrum LPhoton(KdTree<Photon> *map, int nPaths, int nLookup,
    ClosePhoton *lookupBuf, BSDF *bsdf, RNG &rng, const Intersection &isect,
    const Vector &w, float maxDistSquared);
// static Spectrum EPhoton(KdTree<Photon> *map, int count, int nLookup,
    // ClosePhoton *lookupBuf, float maxDist2, const Point &p, const Normal &n);

// PhotonIntegrator Local Definitions
inline float kernel(const Photon *photon, const Point &p,
                    float maxDist2) {
    float s = (1.f - DistanceSquared(photon->p, p) / maxDist2);
    return 3.f * INV_PI * s * s;
}

Spectrum LPhoton(KdTree<Photon> *map, int nPaths, int nLookup,
      ClosePhoton *lookupBuf, BSDF *bsdf, RNG &rng,
      const Intersection &isect, const Vector &wo, float maxDist2) {
    Spectrum L(0.);
    BxDFType nonSpecular = BxDFType(BSDF_REFLECTION |
        BSDF_TRANSMISSION | BSDF_DIFFUSE | BSDF_GLOSSY);
    if (map && bsdf->NumComponents(nonSpecular) > 0) {
        PBRT_PHOTON_MAP_STARTED_LOOKUP(const_cast<DifferentialGeometry *>(&isect.dg));
        // Do photon map lookup at intersection point
        PhotonProcess proc(nLookup, lookupBuf);
        map->Lookup(isect.dg.p, proc, maxDist2);

        // Estimate reflected radiance due to incident photons
        ClosePhoton *photons = proc.photons;
        int nFound = proc.nFound;
        Normal Nf = Faceforward(bsdf->dgShading.nn, wo);
        if (bsdf->NumComponents(BxDFType(BSDF_REFLECTION |
                BSDF_TRANSMISSION | BSDF_GLOSSY)) > 0) {
            // Compute exitant radiance from photons for glossy surface
            for (int i = 0; i < nFound; ++i) {
                const Photon *p = photons[i].photon;
                float k = kernel(p, isect.dg.p, maxDist2);
                L += (k / (nPaths * maxDist2)) * bsdf->f(wo, p->wi) *
                     p->alpha;
            }
        }
        else {
            // Compute exitant radiance from photons for diffuse surface
            Spectrum Lr(0.), Lt(0.);
            for (int i = 0; i < nFound; ++i) {
                if (Dot(Nf, photons[i].photon->wi) > 0.f) {
                    float k = kernel(photons[i].photon, isect.dg.p, maxDist2);
                    Lr += (k / (nPaths * maxDist2)) * photons[i].photon->alpha;
                }
                else {
                    float k = kernel(photons[i].photon, isect.dg.p, maxDist2);
                    Lt += (k / (nPaths * maxDist2)) * photons[i].photon->alpha;
                }
            }
            L += Lr * bsdf->rho(wo, rng, BSDF_ALL_REFLECTION) * INV_PI +
                 Lt * bsdf->rho(wo, rng, BSDF_ALL_TRANSMISSION) * INV_PI;
        }
        PBRT_PHOTON_MAP_FINISHED_LOOKUP(const_cast<DifferentialGeometry *>(&isect.dg),
            proc.nFound, proc.nLookup, &L);
    }
    return L;
}

PhotonIntegrator::~PhotonIntegrator() {
    delete[] lightSampleOffsets;
    delete[] bsdfSampleOffsets;
}
// PhotonIntegrator Method Definitions
PhotonIntegrator::PhotonIntegrator(int nl, int mdepth, float mdist, bool fg,
        int gs, float ga, PhotonShooter* psh) {

    nLookup = nl;
    maxSpecularDepth = mdepth;
    maxDistSquared = mdist * mdist;
    finalGather = fg;
    cosGatherAngle = cos(Radians(ga));
    gatherSamples = gs;
    photonShooter = psh;
    lightSampleOffsets = NULL;
    bsdfSampleOffsets = NULL;
}


void PhotonIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
        const Scene *scene){
    // photonShooter->RequestSamples(sampler,sample,scene);
    // Allocate and request samples for sampling all lights
    uint32_t nLights = scene->lights.size();
    lightSampleOffsets = new LightSampleOffsets[nLights];
    bsdfSampleOffsets = new BSDFSampleOffsets[nLights];
    for (uint32_t i = 0; i < nLights; ++i) {
        const Light *light = scene->lights[i];
        int nSamples = light->nSamples;
        if (sampler) nSamples = sampler->RoundSize(nSamples);
        lightSampleOffsets[i] = LightSampleOffsets(nSamples, sample);
        bsdfSampleOffsets[i] = BSDFSampleOffsets(nSamples, sample);
    }

    // Request samples for final gathering
    if (finalGather) {
        gatherSamples = max(1, gatherSamples/2);
        if (sampler) gatherSamples = sampler->RoundSize(gatherSamples);
        bsdfGatherSampleOffsets = BSDFSampleOffsets(gatherSamples, sample);
        indirGatherSampleOffsets = BSDFSampleOffsets(gatherSamples, sample);
    }
}

Spectrum PhotonIntegrator::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const {
    Spectrum L(0.);

    KdTree<Photon>* causticMap  = photonShooter->causticMap;
    KdTree<Photon>* indirectMap = photonShooter->indirectMap;
    KdTree<RadiancePhoton> *radianceMap = photonShooter->radianceMap;

    int nCausticPaths = photonShooter->nCausticPaths;
    int nIndirectPaths = photonShooter->nIndirectPaths;

    Vector wo = -ray.d;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);

    // Evaluate BSDF at hit point
    BSDF *bsdf = isect.GetBSDF(ray, arena);
    const Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    L += UniformSampleAllLights(scene, renderer, arena, p, n,
        wo, isect.rayEpsilon, ray.time, bsdf, sample, rng,
        lightSampleOffsets, bsdfSampleOffsets);
    // Compute caustic lighting for photon map integrator
    ClosePhoton *lookupBuf = arena.Alloc<ClosePhoton>(nLookup);
    L += LPhoton(causticMap, nCausticPaths, nLookup, lookupBuf, bsdf,
                 rng, isect, wo, maxDistSquared);

    // Compute indirect lighting for photon map integrator
    if (finalGather && indirectMap != NULL) {
    #if 1
        // Do one-bounce final gather for photon map
        BxDFType nonSpecular = BxDFType(BSDF_REFLECTION |
            BSDF_TRANSMISSION | BSDF_DIFFUSE | BSDF_GLOSSY);
        if (bsdf->NumComponents(nonSpecular) > 0) {
            // Find indirect photons around point for importance sampling
            const uint32_t nIndirSamplePhotons = 50;
            PhotonProcess proc(nIndirSamplePhotons,
                               arena.Alloc<ClosePhoton>(nIndirSamplePhotons));
            float searchDist2 = maxDistSquared;
            while (proc.nFound < nIndirSamplePhotons) {
                float md2 = searchDist2;
                proc.nFound = 0;
                indirectMap->Lookup(p, proc, md2);
                searchDist2 *= 2.f;
            }

            // Copy photon directions to local array
            Vector *photonDirs = arena.Alloc<Vector>(nIndirSamplePhotons);
            for (uint32_t i = 0; i < nIndirSamplePhotons; ++i)
                photonDirs[i] = proc.photons[i].photon->wi;

            // Use BSDF to do final gathering
            Spectrum Li = 0.;
            for (int i = 0; i < gatherSamples; ++i) {
                // Sample random direction from BSDF for final gather ray
                Vector wi;
                float pdf;
                BSDFSample bsdfSample(sample, bsdfGatherSampleOffsets, i);
                Spectrum fr = bsdf->Sample_f(wo, &wi, bsdfSample,
                                             &pdf, BxDFType(BSDF_ALL & ~BSDF_SPECULAR));
                if (fr.IsBlack() || pdf == 0.f) continue;
                Assert(pdf >= 0.f);

                // Trace BSDF final gather ray and accumulate radiance
                RayDifferential bounceRay(p, wi, ray, isect.rayEpsilon);
                Intersection gatherIsect;
                if (scene->Intersect(bounceRay, &gatherIsect)) {
                    // Compute exitant radiance _Lindir_ using radiance photons
                    Spectrum Lindir = 0.f;
                    Normal nGather = gatherIsect.dg.nn;
                    nGather = Faceforward(nGather, -bounceRay.d);
                    RadiancePhotonProcess proc(nGather);
                    float md2 = INFINITY;
                    radianceMap->Lookup(gatherIsect.dg.p, proc, md2);
                    if (proc.photon != NULL)
                        Lindir = proc.photon->Lo;
                    Lindir *= renderer->Transmittance(scene, bounceRay, NULL, rng, arena);

                    // Compute MIS weight for BSDF-sampled gather ray

                    // Compute PDF for photon-sampling of direction _wi_
                    float photonPdf = 0.f;
                    float conePdf = UniformConePdf(cosGatherAngle);
                    for (uint32_t j = 0; j < nIndirSamplePhotons; ++j)
                        if (Dot(photonDirs[j], wi) > .999f * cosGatherAngle)
                            photonPdf += conePdf;
                    photonPdf /= nIndirSamplePhotons;
                    float wt = PowerHeuristic(gatherSamples, pdf, gatherSamples, photonPdf);
                    Li += fr * Lindir * (AbsDot(wi, n) * wt / pdf);
                }
            }
            L += Li / gatherSamples;

            // Use nearby photons to do final gathering
            Li = 0.;
            for (int i = 0; i < gatherSamples; ++i) {
                // Sample random direction using photons for final gather ray
                BSDFSample gatherSample(sample, indirGatherSampleOffsets, i);
                int photonNum = min((int)nIndirSamplePhotons - 1,
                    Floor2Int(gatherSample.uComponent * nIndirSamplePhotons));

                // Sample gather ray direction from _photonNum_
                Vector vx, vy;
                CoordinateSystem(photonDirs[photonNum], &vx, &vy);
                Vector wi = UniformSampleCone(gatherSample.uDir[0], gatherSample.uDir[1],
                                              cosGatherAngle, vx, vy, photonDirs[photonNum]);

                // Trace photon-sampled final gather ray and accumulate radiance
                Spectrum fr = bsdf->f(wo, wi);
                if (fr.IsBlack()) continue;
                RayDifferential bounceRay(p, wi, ray, isect.rayEpsilon);
                Intersection gatherIsect;
                PBRT_PHOTON_MAP_STARTED_GATHER_RAY(&bounceRay);
                if (scene->Intersect(bounceRay, &gatherIsect)) {
                    // Compute exitant radiance _Lindir_ using radiance photons
                    Spectrum Lindir = 0.f;
                    Normal nGather = gatherIsect.dg.nn;
                    nGather = Faceforward(nGather, -bounceRay.d);
                    RadiancePhotonProcess proc(nGather);
                    float md2 = INFINITY;
                    radianceMap->Lookup(gatherIsect.dg.p, proc, md2);
                    if (proc.photon != NULL)
                        Lindir = proc.photon->Lo;
                    Lindir *= renderer->Transmittance(scene, bounceRay, NULL, rng, arena);

                    // Compute PDF for photon-sampling of direction _wi_
                    float photonPdf = 0.f;
                    float conePdf = UniformConePdf(cosGatherAngle);
                    for (uint32_t j = 0; j < nIndirSamplePhotons; ++j)
                        if (Dot(photonDirs[j], wi) > .999f * cosGatherAngle)
                            photonPdf += conePdf;
                    photonPdf /= nIndirSamplePhotons;

                    // Compute MIS weight for photon-sampled gather ray
                    float bsdfPdf = bsdf->Pdf(wo, wi);
                    float wt = PowerHeuristic(gatherSamples, photonPdf, gatherSamples, bsdfPdf);
                    Li += fr * Lindir * AbsDot(wi, n) * wt / photonPdf;
                }
                PBRT_PHOTON_MAP_FINISHED_GATHER_RAY(&bounceRay);
            }
            L += Li / gatherSamples;
        }
    #else
        // for debugging / examples: use the photon map directly
        Normal nn = Faceforward(n, -ray.d);
        RadiancePhotonProcess proc(nn);
        float md2 = INFINITY;
        radianceMap->Lookup(p, proc, md2);
        if (proc.photon)
            L += proc.photon->Lo;
    #endif
    }
    else
        L += LPhoton(indirectMap, nIndirectPaths, nLookup, lookupBuf,
                     bsdf, rng, isect, wo, maxDistSquared);
    if (ray.depth+1 < maxSpecularDepth) {
        Vector wi;
        // Trace rays for specular reflection and refraction
        L += SpecularReflect(ray, bsdf, rng, isect, renderer, scene, sample,
                             arena);
        L += SpecularTransmit(ray, bsdf, rng, isect, renderer, scene, sample,
                              arena);
    }
    return L;
}


PhotonIntegrator *CreatePhotonMapSurfaceIntegrator(const ParamSet &params, PhotonShooter* psh) {
    int nUsed = params.FindOneInt("nused", 50);
    if (PbrtOptions.quickRender) nUsed = max(1, nUsed / 10);
    int maxSpecularDepth = params.FindOneInt("maxspeculardepth", 5);
    bool finalGather = params.FindOneBool("finalgather", true);
    int gatherSamples = params.FindOneInt("finalgathersamples", 32);
    if (PbrtOptions.quickRender) gatherSamples = max(1, gatherSamples / 4);
    float maxDist = params.FindOneFloat("maxdist", .1f);
    float gatherAngle = params.FindOneFloat("gatherangle", 10.f);
    return new PhotonIntegrator(nUsed, maxSpecularDepth, maxDist, finalGather, gatherSamples,
        gatherAngle, psh);
}


