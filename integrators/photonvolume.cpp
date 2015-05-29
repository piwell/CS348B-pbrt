#include "stdafx.h"
#include "integrators/photonvolume.h"
#include "renderers/samplerrenderer.h"
// #include "integrators/photonmap.cpp"
#include "paramset.h"
#include "montecarlo.h"

void PhotonVolumeIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
        const Scene *scene){
	tauSampleOffset = sample->Add1D(1);
    scatterSampleOffset = sample->Add1D(1);
}

Spectrum PhotonVolumeIntegrator::Transmittance(const Scene *scene,
        const Renderer *renderer, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena) const {
    if (!scene->volumeRegion) return Spectrum(1.f);
    float step, offset;
    if (sample) {
        step = stepSize;
        offset = sample->oneD[tauSampleOffset][0];
    }
    else {
        step = 4.f * stepSize;
        offset = rng.RandomFloat();
    }
    Spectrum tau = scene->volumeRegion->tau(ray, step, offset);
    return Exp(-tau);
}

// Spectrum PhotonVolumeIntegrator::LPhoton(
// 		KdTree<Photon> *map, int nLookup, VolumeRegion *vr,
// 		const Point &pt, const Vector &wo, float maxDistSquare) const {

Spectrum PhotonVolumeIntegrator::LPhoton(KdTree<Photon> *map, int nLookup, ClosePhoton *lookupBuf,
		const Vector &w, const Point &pt, VolumeRegion *vr, float maxDistSquare, float t) const {

	// printf("LPhoton\n");

	Spectrum L(0.);
	if (!map) return L;

	// Initialize _PhotonProcess_ object, _proc_, for photon map lookups
	PhotonProcess proc(nLookup, lookupBuf);
	proc.photons =
	 	(ClosePhoton *)alloca(nLookup * sizeof(ClosePhoton));
	// // Do photon map lookup

	map->Lookup(pt, proc, maxDistSquare);
	// // Accumulate light from nearby photons

	// // Estimate reflected light from photons
	ClosePhoton *photons = proc.photons;
	int nFound = proc.nFound;

	if (nFound==0)
	 	return L;

	Spectrum totalFlux(0.);

	float maxmd = 0.0;

	for (int i = 0; i < nFound; ++i) {
	 	const Photon *p = photons[i].photon;
	 	Point pt_i = p->p;

	 	float distSq = photons[i].distanceSquared;
	 	if (distSq>maxmd) maxmd = distSq;

	 	totalFlux += p->alpha * vr->p(pt_i,p->wi,-w,t);
	 	//printf("%f !!!\n",distSq);
	}
	// printf("totalflux from photon y() radiance est.: %f\n",totalFlux.y());
	//float distSq = photons[nFound-1].distanceSquared;
	//float distSq = photons[0].distanceSquared;
	float distSq = maxmd;
	float dV = distSq * sqrt(distSq);
	Spectrum scale = vr->sigma_s(pt, w,t);
	if (dV!=0.0 && scale!=0.0)
	 	L += totalFlux/(4.0/3.0*M_PI*dV*scale);
		//if (scale!=0.0)
		//	L += totalFlux/(scale);

	//printf("L.y() for volume: %f\n",L.y());
	return L;
	//return 0;
}


 Spectrum PhotonVolumeIntegrator::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        Spectrum *T, MemoryArena &arena) const {
 	VolumeRegion *vr = scene->volumeRegion;

// 	//this requires exphotonmap to be the surface integrator. Or else.
 	// SamplerRenderer* sr = (SamplerRenderer*)renderer;
 	KdTree<Photon>* volumeMap = photonmap->volumeMap; 
 	// ((PhotonIntegrator*)sr->surfaceIntegrator)->volumeMap;
	// printf("number of nodez in kdtrizzle: %d\n",volumeMap->nNodes);

 	float t0, t1;
 	if (!vr || !vr->IntersectP(ray, &t0, &t1) || (t1-t0) == 0.f){
 		*T = 1.f;
 	 	return 0.f;
 	 }
// 	// Do single scattering & photon multiple scattering volume integration in _vr_
 	Spectrum Lv(0.);


 	// Prepare for volume integration stepping
 	int nSamples = Ceil2Int((t1-t0) / stepSize);
 	// int nSamples = 4*Ceil2Int((t1-t0) / stepSize);
 	float step = (t1 - t0) / nSamples;
 	Spectrum Tr(1.f);
 	Point p = ray(t0), pPrev;
 	// Point p = ray(t1), pPrev;
 	Vector w = -ray.d;
 	t0 += sample->oneD[scatterSampleOffset][0] * step;

 	//if (sample)
 	//	t0 += sample->oneD[scatterSampleOffset][0] * step;
 	//else
 	//	t0 += rng.RandomFloat() * step;
 	// if (sample)
 		// t1 -= sample->oneD[scatterSampleOffset][0] * step;
 	// else
 		// t1 -= rng.RandomFloat() * step;
 	// Compute sample patterns for single scattering samples
 	// float *samp = (float *)alloca(3 * N * sizeof(float));
 	// LatinHypercube(samp, N, 3, rng);
 	float *lightNum = arena.Alloc<float>(nSamples);
    LDShuffleScrambled1D(1, nSamples, lightNum, rng);
    float *lightComp = arena.Alloc<float>(nSamples);
    LDShuffleScrambled1D(1, nSamples, lightComp, rng);
    float *lightPos = arena.Alloc<float>(2*nSamples);
    LDShuffleScrambled2D(1, nSamples, lightPos, rng);
 	int sampOffset = 0;

 	ClosePhoton *lookupBuf = new ClosePhoton[nSamples];

 	for (int i = 0; i < nSamples; ++i, t0 += step) {
 	// for (int i = 0; i < N; ++i, t1 -= step) {
 		// Advance to sample at _t0_ and update _T_
 		pPrev = p;
 		p = ray(t0);

 		//p = ray(t1);
 		Ray tauRay(pPrev, p - pPrev, 0.f, 1.f, ray.time, ray.depth);
 		// Spectrum stepTau = vr->tau(Ray(pPrev, p - pPrev, 0, 1),
 			// .5f * stepSize, rng.RandomFloat());

 		Spectrum stepTau = vr->tau(tauRay,
                                   .5f * stepSize, rng.RandomFloat());
 		// Tr *= Exp(-stepTau);
 		Tr = Exp(-stepTau);


 		// Possibly terminate raymarching if transmittance is small.
 		// if (Tr.y() < 1e-3) {
 			// const float continueProb = .5f;
 			// if (rng.RandomFloat() > continueProb){
 				// Tr = 0.f;
 				// break;
 			// }
 			// Tr /= continueProb;
 		// }
		
		
 		// Compute single-scattering source term at _p_ & photon mapped MS
 		Spectrum L_i(0.);
 		Spectrum L_d(0.);
 		Spectrum L_ii(0.);
 		
 		// Lv += Tr*vr->Lve(p, w, ray.time);
 		Spectrum ss = vr->sigma_s(p, w, ray.time);
 		Spectrum sa = vr->sigma_a(p, w, ray.time);

 		if (!ss.IsBlack() && scene->lights.size() > 0) {
 			int nLights = scene->lights.size();
 			int ln =
 				min(Floor2Int(lightNum[sampOffset] * nLights),
 				    nLights-1);
 			Light *light = scene->lights[ln];
 			// Add contribution of _light_ due to scattering at _p_
 			float pdf;
 			VisibilityTester vis;
 			Vector wo;
 			// float u1 = samp[sampOffset+1], u2 = samp[sampOffset+2];
 			// Spectrum L = light->Sample_L(p, u1, u2, 0, &wo, &pdf, &vis);

 			LightSample ls(lightComp[sampOffset], lightPos[2*sampOffset],
                           lightPos[2*sampOffset+1]);
            Spectrum L = light->Sample_L(p, 0.f, ls, ray.time, &wo, &pdf, &vis);
            

 			if (!L.IsBlack() && pdf > 0.f && vis.Unoccluded(scene)) {
 				Spectrum Ld = L * vis.Transmittance(scene,renderer, NULL, rng, arena);
 				// L_d = vr->p(p, w, -wo, 0) * L_d * float(nLights) / pdf;
 				//printf("L_d.y(): %f\n",L_d.y());
 				L_d = vr->p(p, w, -wo, ray.time) * Ld * float(nLights) /
                        pdf;
 			}
 		}
		// Compute 'indirect' in-scattered radiance from photon map
		 L_ii += LPhoton(volumeMap, nUsed, lookupBuf, w, p, vr, maxDistSquared,ray.time);
		
		// Compute total in-scattered radiance
		if (sa.y()!=0.0 || ss.y()!=0.0)
			L_i = L_d + (ss/(sa+ss))*L_ii;
		else
			L_i = L_d;

		// Lv += Tr * vr->Lve(p, w);

		Spectrum nLv = (sa*vr->Lve(p,w,ray.time)*step) + (ss*L_i*step) + (Tr * Lv);

// 		if (nLv.IsNaN())
// 		{
// 			printf("NaN'ed. Ye olde infodump: L_d.y() = %f, Lve.y() = %f, Tr.y() = %f\n",L_d.y(),(sa*vr->Lve(p,w)).y(),Tr.y());
// 		}

		Lv = nLv;
 		sampOffset++;
 	}
 	*T = Tr;
	return Lv;
}

PhotonVolumeIntegrator *CreatePhotonVolumeIntegrator(const ParamSet &params, const PhotonIntegrator* pm){
	float stepSize  = params.FindOneFloat("stepsize", 1.f);
	int nUsed  = params.FindOneInt("nused", 250);
	float maxDist = params.FindOneFloat("maxdist", 0.1f);
	return new PhotonVolumeIntegrator(stepSize, nUsed, maxDist, pm);
}