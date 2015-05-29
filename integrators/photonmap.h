
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_PHOTONMAP_H
#define PBRT_INTEGRATORS_PHOTONMAP_H

// integrators/photonmap.h*
#include "pbrt.h"
#include "integrator.h"
#include "kdtree.h"
#include "volume.h"
#include "scene.h"
#include "montecarlo.h"

struct Photon {
    Photon(const Point &pp, const Spectrum &wt, const Vector &w)
        : p(pp), alpha(wt), wi(w) { }
    Photon() { }
    Point p;
    Spectrum alpha;
    Vector wi;
};


struct RadiancePhoton {
    RadiancePhoton(const Point &pp, const Normal &nn)
        : p(pp), n(nn), Lo(0.f) { }
    RadiancePhoton() { }
    Point p;
    Normal n;
    Spectrum Lo;
};


struct ClosePhoton {
    // ClosePhoton Public Methods
    ClosePhoton(const Photon *p = NULL, float md2 = INFINITY)
        : photon(p), distanceSquared(md2) { }
    bool operator<(const ClosePhoton &p2) const {
        return distanceSquared == p2.distanceSquared ?
            (photon < p2.photon) : (distanceSquared < p2.distanceSquared);
    }
    const Photon *photon;
    float distanceSquared;
};


struct PhotonProcess {
    // PhotonProcess Public Methods
    PhotonProcess(uint32_t mp, ClosePhoton *buf);
    void operator()(const Point &p, const Photon &photon, float dist2,
                    float &maxDistSquared);
    ClosePhoton *photons;
    uint32_t nLookup, nFound;
};


struct RadiancePhotonProcess {
    // RadiancePhotonProcess Methods
    RadiancePhotonProcess(const Normal &nn)
        :  n(nn) {
        photon = NULL;
    }
    void operator()(const Point &p, const RadiancePhoton &rp,
                    float distSquared, float &maxDistSquared) {
        if (Dot(rp.n, n) > 0) {
            photon = &rp;
            maxDistSquared = distSquared;
        }
    }
    const Normal &n;
    const RadiancePhoton *photon;
};


// PhotonIntegrator Declarations
class PhotonIntegrator : public SurfaceIntegrator {
public:
    // PhotonIntegrator Public Methods
    PhotonIntegrator(int ncaus, int nindir, int nvol, float steps, int nLookup, int maxspecdepth,
        int maxphotondepth, float maxdist, bool finalGather, int gatherSamples,
        float ga);
    ~PhotonIntegrator();
    Spectrum Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect, const Sample *sample,
        RNG &rng, MemoryArena &arena) const;
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    void Preprocess(const Scene *scene, const Camera *camera, const Renderer *renderer);

    uint32_t nVolumePhotonsWanted;
    KdTree<Photon> *volumeMap;
    float stepSize; //for photon-marching to find media interactions


private:
    // PhotonIntegrator Private Methods
    friend class PhotonShootingTask;

    // PhotonIntegrator Private Data
    uint32_t nCausticPhotonsWanted, nIndirectPhotonsWanted, nLookup;
    float maxDistSquared;
    int maxSpecularDepth, maxPhotonDepth;
    bool finalGather;
    int gatherSamples;
    float cosGatherAngle;

    // Declare sample parameters for light source sampling
    LightSampleOffsets *lightSampleOffsets;
    BSDFSampleOffsets *bsdfSampleOffsets;
    BSDFSampleOffsets bsdfGatherSampleOffsets, indirGatherSampleOffsets;
    int nCausticPaths, nIndirectPaths, nVolumePaths;
    KdTree<Photon> *causticMap;
    KdTree<Photon> *indirectMap;
    KdTree<RadiancePhoton> *radianceMap;
};


PhotonIntegrator *CreatePhotonMapSurfaceIntegrator(const ParamSet &params);

#endif // PBRT_INTEGRATORS_PHOTONMAP_H
