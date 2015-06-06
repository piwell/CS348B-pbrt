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


int counter = 0;
float rainbowWavelength(const Vector &w, const Vector &wi){
    //1. Calc wavelength depending on theta
    
    float cosTheta = Dot(wi, -w);
    const float radToDeg = 57.2957;
    float theta = radToDeg * acosf(cosTheta);
    

    //printf("%f\n", theta);
    //printf("%f %f %f\n", p->wi.x, p->wi.y, p->wi.z);
    
    
    //Water droplets modelled as spheres have produce a
    //primary rainbow when angle between incoming and outgoing
    //rays are [40.4, 42.3] degrees.
    //If not in this range, don't contribute to total flux
    const float minTheta = 40.4, minWavelength = 400.0;
    const float maxTheta = 42.3, maxWavelength = 700.0;
    
    
    if (theta < minTheta || maxTheta < theta){
        return 0;
    }
    
    const float thetaRange = maxTheta-minTheta;
    const float wavelengthRange = maxWavelength-minWavelength;
    
    float wavelength = minWavelength + (theta-minTheta) * wavelengthRange / thetaRange;
    
    return wavelength;
}

Spectrum PhotonVolumeIntegrator::LPhoton(KdTree<Photon> *map, int nLookup, ClosePhoton *lookupBuf,
		const Vector &w, const Point &pt, VolumeRegion *vr, float maxDistSquare, float t) const {

	Spectrum L(0.);
	if (!map) return L;

	// Initialize _PhotonProcess_ object, _proc_, for photon map lookups
	PhotonProcess proc(nLookup, lookupBuf);
	proc.photons = (ClosePhoton *)alloca(nLookup * sizeof(ClosePhoton));
	// Do photon map lookup

	map->Lookup(pt, proc, maxDistSquare);
	// // Accumulate light from nearby photons

	// // Estimate reflected light from photons
	ClosePhoton *photons = proc.photons;
	int nFound = proc.nFound;

	if (nFound<10)
	 	return L;

	Spectrum totalFlux(0.);

	float maxmd = 0.0;

	for (int i = 0; i < nFound; ++i) {
	 	const Photon *p = photons[i].photon;
	 	Point pt_i = p->p;

	 	float distSq = photons[i].distanceSquared;
	 	if (distSq>maxmd) maxmd = distSq;
        
        totalFlux += p->alpha * vr->p(pt_i,p->wi,-w,t);
	}
	float distSq = maxmd;
	float dV = distSq * sqrt(distSq);
	Spectrum scale = vr->sigma_s(pt, w,t);
    if (dV!=0.0 && scale!=0.0){
        L += totalFlux/(4.0/3.0*M_PI*dV*scale);
    }
    

	return L;
}



Spectrum PhotonVolumeIntegrator::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        Spectrum *T, MemoryArena &arena) const {
 	
 	VolumeRegion *vr = scene->volumeRegion;
 	KdTree<Photon>* volumeMap = photonShooter->volumeMap; 

 	float t0, t1;
 	if (!vr || !vr->IntersectP(ray, &t0, &t1) || (t1-t0) == 0.f){
 		*T = 1.f;
 	 	return 0.f;
 	 }
 	// Do single scattering & photon multiple scattering volume integration in _vr_
 	Spectrum Lv(0.);


 	// Prepare for volume integration stepping
 	int nSamples = Ceil2Int((t1-t0) / stepSize);
 	float step = (t1 - t0) / nSamples;
 	Spectrum Tr(1.f);
 	Point p = ray(t0), pPrev;
 	Vector w = -ray.d;
 	t0 += sample->oneD[scatterSampleOffset][0] * step;

 	float *lightNum = arena.Alloc<float>(nSamples);
    LDShuffleScrambled1D(1, nSamples, lightNum, rng);
    float *lightComp = arena.Alloc<float>(nSamples);
    LDShuffleScrambled1D(1, nSamples, lightComp, rng);
    float *lightPos = arena.Alloc<float>(2*nSamples);
    LDShuffleScrambled2D(1, nSamples, lightPos, rng);
 	int sampOffset = 0;

 	ClosePhoton *lookupBuf = new ClosePhoton[nSamples];

 	for (int i = 0; i < nSamples; ++i, t0 += step) {
 		// Advance to sample at _t0_ and update _T_
 		pPrev = p;
 		p = ray(t0);

 		Ray tauRay(pPrev, p - pPrev, 0.f, 1.f, ray.time, ray.depth);

 		Spectrum stepTau = vr->tau(tauRay,.5f * stepSize, rng.RandomFloat());
 		Tr = Exp(-stepTau);

 		// Possibly terminate raymarching if transmittance is small.
 		if (Tr.y() < 1e-3) {
 			const float continueProb = .5f;
 			if (rng.RandomFloat() > continueProb){
 				Tr = 0.f;
 				break;
 			}
 			Tr /= continueProb;
 		}
		
		
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

 			LightSample ls(lightComp[sampOffset], lightPos[2*sampOffset],
                           lightPos[2*sampOffset+1]);
            Spectrum L = light->Sample_L(p, 0.f, ls, ray.time, &wo, &pdf, &vis);
            

 			if (!L.IsBlack() && pdf > 0.f && vis.Unoccluded(scene)) {
 				Spectrum Ld = L * vis.Transmittance(scene,renderer, NULL, rng, arena);
 				L_d = vr->p(p, w, -wo, ray.time) * Ld * float(nLights)/pdf;
                
                /* OUR CODE STARTS HERE */
                
                // float wavelength = rainbowWavelength(ray.d, wo);
                // L_d = L_d.filter(wavelength);

                /* OUR CODE ENDS HERE */
 			}
 		}
		// Compute 'indirect' in-scattered radiance from photon map
        
        /* OUR CODE HERE: disabled indirect photon volume integration */
        L_ii += LPhoton(volumeMap, nUsed, lookupBuf, w, p, vr, maxDistSquared,ray.time);
		
        
		// Compute total in-scattered radiance
		if (sa.y()!=0.0 || ss.y()!=0.0)
			L_i = L_d + (ss/(sa+ss))*L_ii;
		else
			L_i = L_d;

		Spectrum nLv = (sa*vr->Lve(p,w,ray.time)*step) + (ss*L_i*step) + (Tr * L)v;

		Lv = nLv;
 		sampOffset++;
 	}
 	*T = Tr;
	return Lv;
}

PhotonVolumeIntegrator *CreatePhotonVolumeIntegrator(const ParamSet &params, PhotonShooter* phs){
	float stepSize  = params.FindOneFloat("stepsize", 1.f);
	int nUsed  = params.FindOneInt("nused", 250);
	float maxDist = params.FindOneFloat("maxdist", 0.1f);
	return new PhotonVolumeIntegrator(stepSize, nUsed, maxDist, phs);
}