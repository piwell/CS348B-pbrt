#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_PHOTONVOLUME_H
#define PBRT_INTEGRATORS_PHOTONVOLUME_H

#include "volume.h"
#include "integrator.h"
#include "scene.h"
#include "integrators/photonmap.h" 

// PhotonVolumeIntegrator Devlarations
class PhotonVolumeIntegrator : public VolumeIntegrator{
public:
	// PhotonVolumeIntegrator Public Methods
	PhotonVolumeIntegrator(float ss, int nu, float md, PhotonShooter* psh = NULL){
		stepSize=ss; nUsed =nu;  maxDist = md; maxDistSquared =md*md;
		photonShooter = psh;
	}
    Spectrum Transmittance(const Scene *, const Renderer *,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena) const;
    void RequestSamples(Sampler *sampler, Sample *sample,
        const Scene *scene);
    Spectrum Li(const Scene *, const Renderer *, const RayDifferential &ray,
         const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena) const;

private:
	Spectrum LPhoton(KdTree<Photon> *map, int nLookup, ClosePhoton *lookupBuf,
		const Vector &w, const Point &pt, VolumeRegion *vr, float maxDistSquare, float t) const;

	float stepSize, maxDist, maxDistSquared;
	int tauSampleOffset, scatterSampleOffset, nUsed;
	PhotonShooter* photonShooter;
};

PhotonVolumeIntegrator *CreatePhotonVolumeIntegrator(const ParamSet &params, PhotonShooter* psh = NULL);

#endif //  PBRT_INTEGRATORS_PHOTONVOLUME_H