
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


// PhotonIntegrator Local Declarations



class PhotonShootingTask : public Task {
public:
    PhotonShootingTask(int tn, float ti, Mutex &m, PhotonIntegrator *in,
        ProgressReporter &prog, bool &at, int &ndp,
        vector<Photon> &direct, vector<Photon> &indir, vector<Photon> &caustic, vector<Photon> &volume,
        vector<RadiancePhoton> &rps, vector<Spectrum> &rpR, vector<Spectrum> &rpT,
        uint32_t &ns, Distribution1D *distrib, const Scene *sc,
        const Renderer *sr)
    : taskNum(tn), time(ti), mutex(m), integrator(in), progress(prog),
      abortTasks(at), nDirectPaths(ndp),
      directPhotons(direct), indirectPhotons(indir), causticPhotons(caustic),volumePhotons(volume),
      radiancePhotons(rps), rpReflectances(rpR), rpTransmittances(rpT),
      nshot(ns), lightDistribution(distrib), scene(sc), renderer (sr) { }
    void Run();

    void followPhoton(RayDifferential photonRay, Intersection photonIsect, Spectrum alpha, int nIntersections, bool specularPath,
        vector<Photon>& localDirectPhotons, vector<Photon>& localIndirectPhotons, vector<Photon>& localCausticPhotons,
        vector<Photon>& localVolumePhotons, vector<RadiancePhoton>& localRadiancePhotons,
        bool& causticDone, bool& indirectDone, bool& volumeDone,
        MemoryArena& arena, RNG& rng,
        vector<Spectrum>& localRpReflectances, vector<Spectrum>& localRpTransmittances);

    int taskNum;
    float time;
    Mutex &mutex;
    PhotonIntegrator *integrator;
    ProgressReporter &progress;
    bool &abortTasks;
    int &nDirectPaths;
    vector<Photon> &directPhotons, &indirectPhotons, &causticPhotons, &volumePhotons;
    vector<RadiancePhoton> &radiancePhotons;
    vector<Spectrum> &rpReflectances, &rpTransmittances;
    uint32_t &nshot;
    const Distribution1D *lightDistribution;
    const Scene *scene;
    const Renderer *renderer;
};


class ComputeRadianceTask : public Task {
public:
    ComputeRadianceTask(ProgressReporter &prog, uint32_t tn, uint32_t nt,
        vector<RadiancePhoton> &rps, const vector<Spectrum> &rhor,
        const vector<Spectrum> &rhot,
        uint32_t nlookup, float md2,
        int ndirect, KdTree<Photon> *direct,
        int nindirect, KdTree<Photon> *indirect,
        int ncaus, KdTree<Photon> *caustic,
        int nvol, KdTree<Photon> *volume)
        : progress(prog), taskNum(tn), numTasks(nt), radiancePhotons(rps),
          rpReflectances(rhor), rpTransmittances(rhot), nLookup(nlookup),
          maxDistSquared(md2),
          nDirectPaths(ndirect), nIndirectPaths(nindirect), nCausticPaths(ncaus),nVolumePaths(nvol),
          directMap(direct), indirectMap(indirect), causticMap(caustic),volumeMap(volume) { }
    void Run();

private:
    ProgressReporter &progress;
    uint32_t taskNum, numTasks;
    vector<RadiancePhoton> &radiancePhotons;
    const vector<Spectrum> &rpReflectances, &rpTransmittances;
    uint32_t nLookup;
    float maxDistSquared;
    int nDirectPaths, nIndirectPaths, nCausticPaths, nVolumePaths;
    KdTree<Photon> *directMap, *indirectMap, *causticMap, *volumeMap;
};


PhotonProcess::PhotonProcess(uint32_t mp, ClosePhoton *buf) {
    photons = buf;
    nLookup = mp;
    nFound = 0;
}


inline float kernel(const Photon *photon, const Point &p, float maxDist2);
static Spectrum LPhoton(KdTree<Photon> *map, int nPaths, int nLookup,
    ClosePhoton *lookupBuf, BSDF *bsdf, RNG &rng, const Intersection &isect,
    const Vector &w, float maxDistSquared);
static Spectrum EPhoton(KdTree<Photon> *map, int count, int nLookup,
    ClosePhoton *lookupBuf, float maxDist2, const Point &p, const Normal &n);

// PhotonIntegrator Local Definitions
inline bool unsuccessful(uint32_t needed, uint32_t found, uint32_t shot) {
    return (found < needed && (found == 0 || found < shot / 1024));
}


inline void PhotonProcess::operator()(const Point &p,
        const Photon &photon, float distSquared, float &maxDistSquared) {
    if (nFound < nLookup) {
        // Add photon to unordered array of photons
        photons[nFound++] = ClosePhoton(&photon, distSquared);
        if (nFound == nLookup) {
            std::make_heap(&photons[0], &photons[nLookup]);
            maxDistSquared = photons[0].distanceSquared;
        }
    }
    else {
        // Remove most distant photon from heap and add new photon
        std::pop_heap(&photons[0], &photons[nLookup]);
        photons[nLookup-1] = ClosePhoton(&photon, distSquared);
        std::push_heap(&photons[0], &photons[nLookup]);
        maxDistSquared = photons[0].distanceSquared;
    }
}


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


Spectrum EPhoton(KdTree<Photon> *map, int count, int nLookup,
        ClosePhoton *lookupBuf, float maxDist2, const Point &p,
        const Normal &n) {
    if (!map) return 0.f;
    // Lookup nearby photons at irradiance computation point
    PhotonProcess proc(nLookup, lookupBuf);
    float md2 = maxDist2;
    map->Lookup(p, proc, md2);
    Assert(md2 > 0.f);

    // Accumulate irradiance value from nearby photons
    if (proc.nFound == 0) return Spectrum(0.f);
    ClosePhoton *photons = proc.photons;
    Spectrum E(0.);
    for (uint32_t i = 0; i < proc.nFound; ++i)
        if (Dot(n, photons[i].photon->wi) > 0.)
            E += photons[i].photon->alpha;
    return E / (count * md2 * M_PI);
}



// PhotonIntegrator Method Definitions
PhotonIntegrator::PhotonIntegrator(int ncaus, int nind, int nvol, float steps,
        int nl, int mdepth, int mphodepth, float mdist, bool fg,
        int gs, float ga) {

    nVolumePhotonsWanted = nvol;
    stepSize = steps;

    nCausticPhotonsWanted = ncaus;
    nIndirectPhotonsWanted = nind;
    nLookup = nl;
    maxSpecularDepth = mdepth;
    maxPhotonDepth = mphodepth;
    maxDistSquared = mdist * mdist;
    finalGather = fg;
    cosGatherAngle = cos(Radians(ga));
    gatherSamples = gs;
    nCausticPaths = nIndirectPaths = nVolumePaths =  0;
    causticMap = indirectMap = volumeMap =  NULL;
    radianceMap = NULL;
    lightSampleOffsets = NULL;
    bsdfSampleOffsets = NULL;
}


PhotonIntegrator::~PhotonIntegrator() {
    delete[] lightSampleOffsets;
    delete[] bsdfSampleOffsets;
    delete causticMap;
    delete indirectMap;
    delete radianceMap;
    delete volumeMap;
}


void PhotonIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
        const Scene *scene) {
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


void PhotonIntegrator::Preprocess(const Scene *scene,
        const Camera *camera, const Renderer *renderer) {
    if (scene->lights.size() == 0) return;
    // Declare shared variables for photon shooting
    Mutex *mutex = Mutex::Create();
    int nDirectPaths = 0;
    vector<Photon> causticPhotons, directPhotons, indirectPhotons;
    vector<RadiancePhoton> radiancePhotons;
    bool abortTasks = false;
    causticPhotons.reserve(nCausticPhotonsWanted);
    indirectPhotons.reserve(nIndirectPhotonsWanted);
    uint32_t nshot = 0;
    vector<Spectrum> rpReflectances, rpTransmittances;

    vector<Photon> volumePhotons;
    volumePhotons.reserve(nVolumePhotonsWanted);

    // Compute light power CDF for photon shooting
    Distribution1D *lightDistribution = ComputeLightSamplingCDF(scene);

    // Run parallel tasks for photon shooting
    ProgressReporter progress(nCausticPhotonsWanted+nIndirectPhotonsWanted+nVolumePhotonsWanted, "Shooting photons");
    vector<Task *> photonShootingTasks;
    int nTasks = NumSystemCores();
    for (int i = 0; i < nTasks; ++i)
        photonShootingTasks.push_back(new PhotonShootingTask(
            i, camera ? camera->shutterOpen : 0.f, *mutex, this, progress, abortTasks, nDirectPaths,
            directPhotons, indirectPhotons, causticPhotons, volumePhotons, radiancePhotons,
            rpReflectances, rpTransmittances,
            nshot, lightDistribution, scene, renderer));
    EnqueueTasks(photonShootingTasks);
    WaitForAllTasks();
    for (uint32_t i = 0; i < photonShootingTasks.size(); ++i)
        delete photonShootingTasks[i];
    Mutex::Destroy(mutex);
    progress.Done();

    // printf("Volume photons: %d, Custic photons: %d \n",volumePhotons.size(), causticPhotons.size());
    // Build kd-trees for indirect and caustic photons
    KdTree<Photon> *directMap = NULL;
    if (directPhotons.size() > 0)
        directMap = new KdTree<Photon>(directPhotons);
    if (causticPhotons.size() > 0)
        causticMap = new KdTree<Photon>(causticPhotons);
    if (indirectPhotons.size() > 0)
        indirectMap = new KdTree<Photon>(indirectPhotons);
    if (volumePhotons.size() > 0)
        volumeMap = new KdTree<Photon>(volumePhotons);

    // Precompute radiance at a subset of the photons
    if (finalGather && radiancePhotons.size()) {
        // Launch tasks to compute photon radiances
        vector<Task *> radianceTasks;
        uint32_t numTasks = 64;
        ProgressReporter progRadiance(numTasks, "Computing photon radiances");
        for (uint32_t i = 0; i < numTasks; ++i)
            radianceTasks.push_back(new ComputeRadianceTask(progRadiance,
                i, numTasks, radiancePhotons, rpReflectances, rpTransmittances,
                nLookup, maxDistSquared, nDirectPaths, directMap,
                nIndirectPaths, indirectMap,
                nCausticPaths, causticMap,
                nVolumePaths, volumeMap));
        EnqueueTasks(radianceTasks);
        WaitForAllTasks();
        for (uint32_t i = 0; i < radianceTasks.size(); ++i)
            delete radianceTasks[i];
        progRadiance.Done();
        radianceMap = new KdTree<RadiancePhoton>(radiancePhotons);
    }
    delete directMap;
}


void PhotonShootingTask::followPhoton(RayDifferential photonRay, Intersection photonIsect, Spectrum alpha, int nIntersections, bool specularPath,
        vector<Photon>& localDirectPhotons, vector<Photon>& localIndirectPhotons, vector<Photon>& localCausticPhotons,
        vector<Photon>& localVolumePhotons, vector<RadiancePhoton>& localRadiancePhotons,
        bool& causticDone, bool& indirectDone, bool& volumeDone,
        MemoryArena& arena, RNG& rng,
        vector<Spectrum>& localRpReflectances, vector<Spectrum>& localRpTransmittances){
        
        if (scene->Intersect(photonRay, &photonIsect)) {
            ++nIntersections;
            //---------------------------------------------------------------------
            //figure out where the volume is in the current photon's path
            float t0, t1;
            float length = photonRay.d.Length();
            if (length == 0.f) return;  //shouldn't happen
            Ray rn(photonRay.o, photonRay.d / length, photonRay.mint * length, photonRay.maxt * length);
            if (!scene->volumeRegion->IntersectP(rn, &t0, &t1)) {/*printf("failure to volumize!\n");*/t0 = 1.0; t1 = 0.0;} //no volumes were intersected by ray
            //else printf("Didn't fail\n");
                Spectrum tau(0.);
                t0 += rng.RandomFloat() * integrator->stepSize; //I am always stopping myself from typing "RootBeerFloat()"
                float t_i = t0;


                // //march through the volume and find the pdf of where the event will occur
                float xi = rng.RandomFloat();
                bool interaction = false;
                
                while (t0 < t1) {
                    //tau += scene->volumeRegion->sigma_t(rn(t0), rn.d*integrator->stepSize,rn.time);
                    //printf("1-exp(-blau) = %f\n",exp(-tau.y()));
                    //if (xi > ( exp(-tau.y()) )) //use the Y in XYZ (luminance) as a termination quantity
                    RayDifferential shortRay(photonRay.o, rn.d, t_i, t0);
                    Spectrum tr = renderer->Transmittance(scene, shortRay, NULL, rng, arena);
                    //scene->Transmittance(shortRay);
                    if (xi > tr.y()){
                        interaction = true;
                        break;
                    }
                    t0 += integrator->stepSize;
                }
                
                if (interaction){
                //     //printf("Interaction! Not human interaction, obviously, but interaction nonetheless!\n");
                     Point interactPt = rn(t0);
                    //figure out if it's absorbed or scattered
                    Spectrum sig_s = scene->volumeRegion->sigma_s(interactPt, rn.d, rn.time);
                    Spectrum sig_a = scene->volumeRegion->sigma_a(interactPt, rn.d, rn.time);
                    bool scatter = (rng.RandomFloat() > (sig_s.y())/(sig_a.y()+sig_s.y()));
                    
                    //if it's absorbed, terminate with extreme prejudice
                    if (!scatter){
                         //store the photon? no
                         return;
                    }
                    
                    //if it's scattered, importance sample the phase function and multiply it by the phase function and change its direction
                    if (scatter && !volumeDone){
                        //if nIntersections>1, store. Else: discard, it will be accounted for in the volume integrator in direct lighting single scattering.
                        if (nIntersections>1){
                             Photon photon(interactPt, alpha, rn.d);
                             localVolumePhotons.push_back(photon);
                        }else{
                            integrator->nVolumePaths++;
                        }


                        // get a uniform sphere direction sample
                        float u1 = rng.RandomFloat();
                        float u2 = rng.RandomFloat();
                        Vector direction = UniformSampleSphere(u1,u2);
                        float pdf = UniformSpherePdf();
                        
                        //check out what the phase function thinks of your new direction
                        Spectrum ref = scene->volumeRegion->p(interactPt,rn.d,direction,rn.time);

                        //make a new ray for the photon; continue outer loop.
                        if (ref.IsBlack() || pdf == 0.f)
                            return;
                        alpha *= ref;
                        alpha /= pdf;
                        
                        
                        photonRay = RayDifferential(interactPt, direction,0);
                        followPhoton(photonRay, photonIsect, alpha, nIntersections, specularPath,
                        localDirectPhotons, localIndirectPhotons, localCausticPhotons, localVolumePhotons, localRadiancePhotons,
                        causticDone, indirectDone, volumeDone, arena, rng, rpReflectances, rpTransmittances);
                    }
                }
                //--------------------------------------------------------------------


                vector<Spectrum> spectrums;
                // Handle photon/surface intersection
                alpha *= renderer->Transmittance(scene, photonRay, NULL, rng, arena);
                BSDF *photonBSDF = photonIsect.GetBSDF(photonRay, arena);
                BxDFType specularType = BxDFType(BSDF_REFLECTION |
                                        BSDF_TRANSMISSION | BSDF_SPECULAR);
                bool hasNonSpecular = (photonBSDF->NumComponents() >
                                       photonBSDF->NumComponents(specularType));
                
                bool hasTransmission = (photonBSDF->NumComponents(BxDFType(BSDF_ALL_TRANSMISSION))>0);
                if(hasTransmission && !alpha.monochromatic){
                    alpha.splitSpectrum(spectrums);
                    // alpha = spectrums[2];
                }else{
                    spectrums.push_back(alpha);
                }

                Vector wo = -photonRay.d;
                if (hasNonSpecular) {
                    // Deposit photon at surface
                    Photon photon(photonIsect.dg.p, alpha, wo);
                    bool depositedPhoton = false;
                    if (specularPath && nIntersections > 1) {
                        if (!causticDone) {
                            PBRT_PHOTON_MAP_DEPOSITED_CAUSTIC_PHOTON(&photonIsect.dg, &alpha, &wo);
                            depositedPhoton = true;
                            localCausticPhotons.push_back(photon);
                        }
                    }
                    else {
                        // Deposit either direct or indirect photon
                        // stop depositing direct photons once indirectDone is true; don't
                        // want to waste memory storing too many if we're going a long time
                        // trying to get enough caustic photons desposited.
                        if (nIntersections == 1 && !indirectDone && integrator->finalGather) {
                            PBRT_PHOTON_MAP_DEPOSITED_DIRECT_PHOTON(&photonIsect.dg, &alpha, &wo);
                            depositedPhoton = true;
                            localDirectPhotons.push_back(photon);
                        }
                        else if (nIntersections > 1 && !indirectDone) {
                            PBRT_PHOTON_MAP_DEPOSITED_INDIRECT_PHOTON(&photonIsect.dg, &alpha, &wo);
                            depositedPhoton = true;
                            localIndirectPhotons.push_back(photon);
                        }
                    }

                    // Possibly create radiance photon at photon intersection point
                    if (depositedPhoton && integrator->finalGather &&
                            rng.RandomFloat() < .125f) {
                        Normal n = photonIsect.dg.nn;
                        n = Faceforward(n, -photonRay.d);
                        localRadiancePhotons.push_back(RadiancePhoton(photonIsect.dg.p, n));
                        Spectrum rho_r = photonBSDF->rho(rng, BSDF_ALL_REFLECTION);
                        localRpReflectances.push_back(rho_r);
                        Spectrum rho_t = photonBSDF->rho(rng, BSDF_ALL_TRANSMISSION);
                        localRpTransmittances.push_back(rho_t);
                    }
                }
                if (nIntersections >= integrator->maxPhotonDepth) return;

                // Sample new photon ray direction
                    Vector wi;
                    float pdf;
                    BxDFType flags;

                    // for(uint32_t i=0; i<spectrums.size(); ++i){
                    // alpha = spectrums[i]
                    alpha = spectrums[0];

                    Spectrum fr = photonBSDF->Sample_f(wo, &wi, BSDFSample(rng),
                                                       &pdf, BSDF_ALL, &flags, &alpha);



                    if (fr.IsBlack() || pdf == 0.f) return;
                    Spectrum anew = alpha * fr *
                        AbsDot(wi, photonBSDF->dgShading.nn) / pdf;

                    // Possibly terminate photon path with Russian roulette
                    float continueProb = min(1.f, anew.y() / alpha.y());
                    if (rng.RandomFloat() > continueProb)
                        return;

                    alpha = anew / continueProb;
                    specularPath &= ((flags & BSDF_SPECULAR) != 0);
                    
                    if (indirectDone && !specularPath) return;
                    photonRay = RayDifferential(photonIsect.dg.p, wi, photonRay,
                                                photonIsect.rayEpsilon);
                    followPhoton(photonRay, photonIsect, alpha, nIntersections, specularPath,
                        localDirectPhotons, localIndirectPhotons, localCausticPhotons, localVolumePhotons, localRadiancePhotons,
                        causticDone, indirectDone, volumeDone, arena, rng, rpReflectances, rpTransmittances);
                // }

                // followPhoton(RayDifferential photonRay, Intersection photonIsect, Spectrum alpha, int nIntersections, bool specularPath,
        // vector<Photon>& localDirectPhotons, vector<Photon>& localIndirectPhotons, vector<Photon>& localCausticPhotons,
        // vector<Photon>& localVolumePhotons, vector<RadiancePhoton>& localRadiancePhotons,
        // bool& causticDone, bool& indirectDone, bool& volumeDone, int& totalPaths,
        // MemoryArena* arena, RNG* rng)
        }


    }


void PhotonShootingTask::Run() {
    // Declare local variables for _PhotonShootingTask_
    // printf("PhotonShootingTask\n");
    MemoryArena arena;
    RNG rng(31 * taskNum);
    vector<Photon> localDirectPhotons, localIndirectPhotons, localCausticPhotons, localVolumePhotons;
    vector<RadiancePhoton> localRadiancePhotons;
    uint32_t totalPaths = 0;
    bool causticDone    = (integrator->nCausticPhotonsWanted == 0);
    bool indirectDone   = (integrator->nIndirectPhotonsWanted == 0);
    bool volumeDone     = (integrator->nVolumePhotonsWanted == 0);

    PermutedHalton halton(6, rng);
    vector<Spectrum> localRpReflectances, localRpTransmittances;
    while (true) {
        // Follow photon paths for a block of samples
        const uint32_t blockSize = 4096;
        for (uint32_t i = 0; i < blockSize; ++i) {
            float u[6];
            halton.Sample(++totalPaths, u);
            // Choose light to shoot photon from
            float lightPdf;
            int lightNum = lightDistribution->SampleDiscrete(u[0], &lightPdf);
            const Light *light = scene->lights[lightNum];

            // Generate _photonRay_ from light source and initialize _alpha_
            RayDifferential photonRay;
            float pdf;
            LightSample ls(u[1], u[2], u[3]);
            Normal Nl;
            Spectrum Le = light->Sample_L(scene, ls, u[4], u[5],
                                          time, &photonRay, &Nl, &pdf);
            if (pdf == 0.f || Le.IsBlack()) continue;
            Spectrum alpha = (AbsDot(Nl, photonRay.d) * Le) / (pdf * lightPdf);
            if (!alpha.IsBlack()) {
                // Follow photon path through scene and record intersections
                PBRT_PHOTON_MAP_STARTED_RAY_PATH(&photonRay, &alpha);
                bool specularPath = true;
                Intersection photonIsect;
                int nIntersections = 0;
                followPhoton(photonRay, photonIsect, alpha, nIntersections, specularPath,
                    localDirectPhotons, localIndirectPhotons, localCausticPhotons, localVolumePhotons, localRadiancePhotons,
                    causticDone, indirectDone, volumeDone, arena, rng, rpReflectances, rpTransmittances);
                
                // while (scene->Intersect(photonRay, &photonIsect)) {
                //     ++nIntersections;


                //     //---------------------------------------------------------------------
                //     //figure out where the volume is in the current photon's path
                //     float t0, t1;
                //     float length = photonRay.d.Length();
                //     if (length == 0.f) break;  //shouldn't happen
                //     Ray rn(photonRay.o, photonRay.d / length, photonRay.mint * length, photonRay.maxt * length);
                //     if (!scene->volumeRegion->IntersectP(rn, &t0, &t1)) {/*printf("failure to volumize!\n");*/t0 = 1.0; t1 = 0.0;} //no volumes were intersected by ray
                //     //else printf("Didn't fail\n");
                //         Spectrum tau(0.);
                //         t0 += rng.RandomFloat() * integrator->stepSize; //I am always stopping myself from typing "RootBeerFloat()"
                //         float t_i = t0;


                //         // //march through the volume and find the pdf of where the event will occur
                //         float xi = rng.RandomFloat();
                //         bool interaction = false;
                        
                //         while (t0 < t1) {
                //             //tau += scene->volumeRegion->sigma_t(rn(t0), rn.d*integrator->stepSize,rn.time);
                //             //printf("1-exp(-blau) = %f\n",exp(-tau.y()));
                //             //if (xi > ( exp(-tau.y()) )) //use the Y in XYZ (luminance) as a termination quantity
                //             RayDifferential shortRay(photonRay.o, rn.d, t_i, t0);
                //             Spectrum tr = renderer->Transmittance(scene, shortRay, NULL, rng, arena);
                //             //scene->Transmittance(shortRay);
                //             if (xi > tr.y()){
                //                 interaction = true;
                //                 break;
                //             }
                //             t0 += integrator->stepSize;
                //         }
                        
                //         if (interaction){
                //         //     //printf("Interaction! Not human interaction, obviously, but interaction nonetheless!\n");
                //              Point interactPt = rn(t0);
                //             //figure out if it's absorbed or scattered
                //             Spectrum sig_s = scene->volumeRegion->sigma_s(interactPt, rn.d, rn.time);
                //             Spectrum sig_a = scene->volumeRegion->sigma_a(interactPt, rn.d, rn.time);
                //             bool scatter = (rng.RandomFloat() > (sig_s.y())/(sig_a.y()+sig_s.y()));
                            
                //             //if it's absorbed, terminate with extreme prejudice
                //             if (!scatter){
                //                  //store the photon? no
                //                  break;
                //             }
                            
                //             //if it's scattered, importance sample the phase function and multiply it by the phase function and change its direction
                //             if (scatter && !volumeDone){
                //                 //if nIntersections>1, store. Else: discard, it will be accounted for in the volume integrator in direct lighting single scattering.
                //                 if (nIntersections>1){
                //                      Photon photon(interactPt, alpha, rn.d);
                //                      localVolumePhotons.push_back(photon);
                //                 }else{
                //                     integrator->nVolumePaths++;
                //                 }


                //                 // get a uniform sphere direction sample
                //                 float u1 = rng.RandomFloat();
                //                 float u2 = rng.RandomFloat();
                //                 Vector direction = UniformSampleSphere(u1,u2);
                //                 float pdf = UniformSpherePdf();
                                
                //                 //check out what the phase function thinks of your new direction
                //                 Spectrum ref = scene->volumeRegion->p(interactPt,rn.d,direction,rn.time);

                //                 //make a new ray for the photon; continue outer loop.
                //                 if (ref.IsBlack() || pdf == 0.f)
                //                     break;
                //                 alpha *= ref;
                //                 alpha /= pdf;
                                
                                
                //                 photonRay = RayDifferential(interactPt, direction,0);
                //                 continue;
                //             }
                //         }
                //         //--------------------------------------------------------------------


                //         // Handle photon/surface intersection
                //         alpha *= renderer->Transmittance(scene, photonRay, NULL, rng, arena);
                //         BSDF *photonBSDF = photonIsect.GetBSDF(photonRay, arena);
                //         BxDFType specularType = BxDFType(BSDF_REFLECTION |
                //                                 BSDF_TRANSMISSION | BSDF_SPECULAR);
                //         bool hasNonSpecular = (photonBSDF->NumComponents() >
                //                                photonBSDF->NumComponents(specularType));
                        
                //         bool hasTransmission = (photonBSDF->NumComponents(BxDFType(BSDF_ALL_TRANSMISSION))>0);
                //         if(hasTransmission){
                //             vector<Spectrum> spectrums;
                //                 if(!alpha.monochromatic){
                //                     alpha.splitSpectrum(spectrums);
                //                     alpha = spectrums[2];
                //             }else{
                //                 spectrums.push_back(alpha);
                //             }
                //         }

                //         Vector wo = -photonRay.d;
                //         if (hasNonSpecular) {
                //             // Deposit photon at surface
                //             Photon photon(photonIsect.dg.p, alpha, wo);
                //             bool depositedPhoton = false;
                //             if (specularPath && nIntersections > 1) {
                //                 if (!causticDone) {
                //                     PBRT_PHOTON_MAP_DEPOSITED_CAUSTIC_PHOTON(&photonIsect.dg, &alpha, &wo);
                //                     depositedPhoton = true;
                //                     localCausticPhotons.push_back(photon);
                //                 }
                //             }
                //             else {
                //                 // Deposit either direct or indirect photon
                //                 // stop depositing direct photons once indirectDone is true; don't
                //                 // want to waste memory storing too many if we're going a long time
                //                 // trying to get enough caustic photons desposited.
                //                 if (nIntersections == 1 && !indirectDone && integrator->finalGather) {
                //                     PBRT_PHOTON_MAP_DEPOSITED_DIRECT_PHOTON(&photonIsect.dg, &alpha, &wo);
                //                     depositedPhoton = true;
                //                     localDirectPhotons.push_back(photon);
                //                 }
                //                 else if (nIntersections > 1 && !indirectDone) {
                //                     PBRT_PHOTON_MAP_DEPOSITED_INDIRECT_PHOTON(&photonIsect.dg, &alpha, &wo);
                //                     depositedPhoton = true;
                //                     localIndirectPhotons.push_back(photon);
                //                 }
                //             }

                //             // Possibly create radiance photon at photon intersection point
                //             if (depositedPhoton && integrator->finalGather &&
                //                     rng.RandomFloat() < .125f) {
                //                 Normal n = photonIsect.dg.nn;
                //                 n = Faceforward(n, -photonRay.d);
                //                 localRadiancePhotons.push_back(RadiancePhoton(photonIsect.dg.p, n));
                //                 Spectrum rho_r = photonBSDF->rho(rng, BSDF_ALL_REFLECTION);
                //                 localRpReflectances.push_back(rho_r);
                //                 Spectrum rho_t = photonBSDF->rho(rng, BSDF_ALL_TRANSMISSION);
                //                 localRpTransmittances.push_back(rho_t);
                //             }
                //         }
                //         if (nIntersections >= integrator->maxPhotonDepth) break;

                //         // Sample new photon ray direction
                //         Vector wi;
                //         float pdf;
                //         BxDFType flags;
        
                //         Spectrum fr = photonBSDF->Sample_f(wo, &wi, BSDFSample(rng),
                //                                            &pdf, BSDF_ALL, &flags, &alpha);
                //         if (fr.IsBlack() || pdf == 0.f) break;
                //         Spectrum anew = alpha * fr *
                //             AbsDot(wi, photonBSDF->dgShading.nn) / pdf;

                //         // Possibly terminate photon path with Russian roulette
                //         float continueProb = min(1.f, anew.y() / alpha.y());
                //         if (rng.RandomFloat() > continueProb)
                //             break;
                //         alpha = anew / continueProb;
                //         specularPath &= ((flags & BSDF_SPECULAR) != 0);
                        
                //         if (indirectDone && !specularPath) break;
                //         photonRay = RayDifferential(photonIsect.dg.p, wi, photonRay,
                //                                     photonIsect.rayEpsilon);
                // }
                PBRT_PHOTON_MAP_FINISHED_RAY_PATH(&photonRay, &alpha);
            }
            arena.FreeAll();
        }

        // Merge local photon data with data in _PhotonIntegrator_
        { MutexLock lock(mutex);

        // Give up if we're not storing enough photons
        if (abortTasks)
            return;
        if (nshot > 500000 &&
            (unsuccessful(integrator->nCausticPhotonsWanted,
                                      causticPhotons.size(), blockSize) ||
             unsuccessful(integrator->nIndirectPhotonsWanted,
                                      indirectPhotons.size(), blockSize)||
             unsuccessful(integrator->nVolumePhotonsWanted,
                                      volumePhotons.size(), blockSize))) {
            Error("Unable to store enough photons.  Giving up.\n");
            causticPhotons.erase(causticPhotons.begin(), causticPhotons.end());
            indirectPhotons.erase(indirectPhotons.begin(), indirectPhotons.end());
            volumePhotons.erase(volumePhotons.begin(),volumePhotons.end());
            radiancePhotons.erase(radiancePhotons.begin(), radiancePhotons.end());
            abortTasks = true;
            return;
        }
        progress.Update(localIndirectPhotons.size() + localCausticPhotons.size() + localVolumePhotons.size());
        nshot += blockSize;

        // Merge indirect photons into shared array
        if (!indirectDone) {
            integrator->nIndirectPaths += blockSize;
            for (uint32_t i = 0; i < localIndirectPhotons.size(); ++i)
                indirectPhotons.push_back(localIndirectPhotons[i]);
            localIndirectPhotons.erase(localIndirectPhotons.begin(),
                                       localIndirectPhotons.end());
            if (indirectPhotons.size() >= integrator->nIndirectPhotonsWanted)
                indirectDone = true;
            nDirectPaths += blockSize;
            for (uint32_t i = 0; i < localDirectPhotons.size(); ++i)
                directPhotons.push_back(localDirectPhotons[i]);
            localDirectPhotons.erase(localDirectPhotons.begin(),
                                     localDirectPhotons.end());
        }

        // Merge direct, caustic, and radiance photons into shared array
        if (!causticDone) {
            integrator->nCausticPaths += blockSize;
            for (uint32_t i = 0; i < localCausticPhotons.size(); ++i)
                causticPhotons.push_back(localCausticPhotons[i]);
            localCausticPhotons.erase(localCausticPhotons.begin(), localCausticPhotons.end());
            if (causticPhotons.size() >= integrator->nCausticPhotonsWanted)
                causticDone = true;
        }

        if(!volumeDone){
            integrator->nVolumePaths += blockSize;
            for(uint32_t i = 0; i < localVolumePhotons.size(); ++i){
                localVolumePhotons[i].alpha /= float(nshot);
                volumePhotons.push_back(localVolumePhotons[i]);
            }
            // printf("Volume photons: %d \n", volumePhotons.size());
            localVolumePhotons.erase(localVolumePhotons.begin(), localVolumePhotons.end());
            if(volumePhotons.size() >= integrator->nVolumePhotonsWanted)
                volumeDone = true;
        }
        
        for (uint32_t i = 0; i < localRadiancePhotons.size(); ++i)
            radiancePhotons.push_back(localRadiancePhotons[i]);
        localRadiancePhotons.erase(localRadiancePhotons.begin(), localRadiancePhotons.end());
        for (uint32_t i = 0; i < localRpReflectances.size(); ++i)
            rpReflectances.push_back(localRpReflectances[i]);
        localRpReflectances.erase(localRpReflectances.begin(), localRpReflectances.end());
        for (uint32_t i = 0; i < localRpTransmittances.size(); ++i)
            rpTransmittances.push_back(localRpTransmittances[i]);
        localRpTransmittances.erase(localRpTransmittances.begin(), localRpTransmittances.end());
        }

        // Exit task if enough photons have been found
        if (indirectDone && causticDone && volumeDone)
            break;
    }
}

void ComputeRadianceTask::Run() {
    // Compute range of radiance photons to process in task
    uint32_t taskSize = radiancePhotons.size() / numTasks;
    uint32_t excess = radiancePhotons.size() % numTasks;
    uint32_t rpStart = min(taskNum, excess) * (taskSize+1) +
                       max(0, (int)taskNum-(int)excess) * taskSize;
    uint32_t rpEnd = rpStart + taskSize + (taskNum < excess ? 1 : 0);
    if (taskNum == numTasks-1) Assert(rpEnd == radiancePhotons.size());
    ClosePhoton *lookupBuf = new ClosePhoton[nLookup];
    for (uint32_t i = rpStart; i < rpEnd; ++i) {
        // Compute radiance for radiance photon _i_
        RadiancePhoton &rp = radiancePhotons[i];
        const Spectrum &rho_r = rpReflectances[i], &rho_t = rpTransmittances[i];
        if (!rho_r.IsBlack()) {
            // Accumulate outgoing radiance due to reflected irradiance
            Spectrum E = EPhoton(directMap, nDirectPaths, nLookup, lookupBuf,
                                 maxDistSquared, rp.p, rp.n) +
                         EPhoton(indirectMap, nIndirectPaths, nLookup, lookupBuf,
                                 maxDistSquared, rp.p, rp.n) +
                         EPhoton(causticMap, nCausticPaths, nLookup, lookupBuf,
                                 maxDistSquared, rp.p, rp.n);
            rp.Lo += INV_PI * rho_r * E;
        }
        if (!rho_t.IsBlack()) {
            // Accumulate outgoing radiance due to transmitted irradiance
            Spectrum E = EPhoton(directMap, nDirectPaths, nLookup, lookupBuf,
                                 maxDistSquared, rp.p, -rp.n) +
                         EPhoton(indirectMap, nIndirectPaths, nLookup, lookupBuf,
                                 maxDistSquared, rp.p, -rp.n) +
                         EPhoton(causticMap, nCausticPaths, nLookup, lookupBuf,
                                 maxDistSquared, rp.p, -rp.n);
            rp.Lo += INV_PI * rho_t * E;
        }
    }
    delete[] lookupBuf;
    progress.Update();
}


Spectrum PhotonIntegrator::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const {
    Spectrum L(0.);
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


PhotonIntegrator *CreatePhotonMapSurfaceIntegrator(const ParamSet &params) {
    int nCaustic = params.FindOneInt("causticphotons", 20000);
    int nIndirect = params.FindOneInt("indirectphotons", 100000);
    int nVolume = params.FindOneInt("volumephotons", 0);
    float stepSize = params.FindOneFloat("stepsize", 0.1f);
    int nUsed = params.FindOneInt("nused", 50);
    if (PbrtOptions.quickRender) nCaustic = nCaustic / 10;
    if (PbrtOptions.quickRender) nIndirect = nIndirect / 10;
    if (PbrtOptions.quickRender) nUsed = max(1, nUsed / 10);
    int maxSpecularDepth = params.FindOneInt("maxspeculardepth", 5);
    int maxPhotonDepth = params.FindOneInt("maxphotondepth", 5);
    bool finalGather = params.FindOneBool("finalgather", true);
    int gatherSamples = params.FindOneInt("finalgathersamples", 32);
    if (PbrtOptions.quickRender) gatherSamples = max(1, gatherSamples / 4);
    float maxDist = params.FindOneFloat("maxdist", .1f);
    float gatherAngle = params.FindOneFloat("gatherangle", 10.f);
    return new PhotonIntegrator(nCaustic, nIndirect, nVolume, stepSize,
        nUsed, maxSpecularDepth, maxPhotonDepth, maxDist, finalGather, gatherSamples,
        gatherAngle);
}


