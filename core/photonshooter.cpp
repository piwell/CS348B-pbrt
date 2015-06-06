// core/photnshooter.cpp*
#include "stdafx.h"
#include "scene.h"
#include "montecarlo.h"
#include "sampler.h"
#include "progressreporter.h"
#include "intersection.h"
#include "paramset.h"
#include "camera.h"
#include "spectrum.h"
#include "photonshooter.h"


// static Spectrum EPhoton(KdTree<Photon> *map, int count, int nLookup,
//     ClosePhoton *lookupBuf, float maxDist2, const Point &p, const Normal &n);

static Spectrum EPhoton(KdTree<Photon> *map, int count, int nLookup,
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

inline bool unsuccessful(uint32_t needed, uint32_t found, uint32_t shot) {
    return (found < needed && (found == 0 || found < shot / 1024));
}

PhotonProcess::PhotonProcess(uint32_t mp, ClosePhoton *buf) {
    photons = buf;
    nLookup = mp;
    nFound = 0;
}

void PhotonShootingTask::followPhoton(RayDifferential photonRay, Intersection photonIsect, Spectrum alpha, int nIntersections, bool specularPath,
vector<Photon>& localDirectPhotons, vector<Photon>& localIndirectPhotons, vector<Photon>& localCausticPhotons,
vector<Photon>& localVolumePhotons, vector<RadiancePhoton>& localRadiancePhotons,
bool& causticDone, bool& indirectDone, bool& volumeDone,
MemoryArena* arena, RNG* rng,
vector<Spectrum>& localRpReflectances, vector<Spectrum>& localRpTransmittances){

if (scene->Intersect(photonRay, &photonIsect)) {
    ++nIntersections;
    //figure out where the volume is in the current photon's path
    float t0, t1;
    float length = photonRay.d.Length();
    if (length == 0.f) return;  //shouldn't happen
    Ray rn(photonRay.o, photonRay.d / length, photonRay.mint * length, photonRay.maxt * length);
    if (!scene->volumeRegion->IntersectP(rn, &t0, &t1)) {/*printf("failure to volumize!\n");*/t0 = 1.0; t1 = 0.0;} 
        Spectrum tau(0.);
        t0 += rng->RandomFloat() * shooter->stepSize;
        float t_i = t0;


        // //march through the volume and find the pdf of where the event will occur
        float xi = rng->RandomFloat();
        bool interaction = false;
        
        while (t0 < t1) {
            RayDifferential shortRay(photonRay.o, rn.d, t_i, t0);
            Spectrum tr = renderer->Transmittance(scene, shortRay, NULL, *rng, *arena);
            //scene->Transmittance(shortRay);
            if (xi > tr.y()){
                interaction = true;
                break;
            }
            t0 += shooter->stepSize;
        }
        
        if (interaction){
            Point interactPt = rn(t0);

            //figure out if it's absorbed or scattered
            Spectrum sig_s = scene->volumeRegion->sigma_s(interactPt, rn.d, rn.time);
            Spectrum sig_a = scene->volumeRegion->sigma_a(interactPt, rn.d, rn.time);
            bool scatter = (rng->RandomFloat() > (sig_s.y())/(sig_a.y()+sig_s.y()));
            
            //if it's absorbed, terminate with extreme prejudice
            if (!scatter){
                 return;
            }
            
            //if it's scattered, importance sample the phase function and multiply it by the phase function and change its direction
            if (scatter && !volumeDone){
                //if nIntersections>1, store. Else: discard, it will be accounted for in the volume integrator in direct lighting single scattering.
                if (nIntersections>1){
                    Photon photon(interactPt, alpha, rn.d);
                    { MutexLock lock(mutex);
                     localVolumePhotons.push_back(photon);
                    }
                }else{
                    shooter->nVolumePaths++;
                }

                // get a uniform sphere direction sample
                float u1 = rng->RandomFloat();
                float u2 = rng->RandomFloat();
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


        vector<Spectrum> spectrums;
        // Handle photon/surface intersection
        alpha *= renderer->Transmittance(scene, photonRay, NULL, *rng, *arena);
        BSDF *photonBSDF = photonIsect.GetBSDF(photonRay, *arena);
        BxDFType specularType = BxDFType(BSDF_REFLECTION |
                                BSDF_TRANSMISSION | BSDF_SPECULAR);
        bool hasNonSpecular = (photonBSDF->NumComponents() >
                               photonBSDF->NumComponents(specularType));
        
        bool hasTransmission = (photonBSDF->NumComponents(BxDFType(BSDF_ALL_TRANSMISSION))>0);
        if(hasTransmission && alpha.lambda<0 && photonIsect.primitive->dispersive()){
            alpha.splitSpectrum(spectrums);
        }else{
            spectrums.push_back(alpha);
        }

        Vector wo = -photonRay.d;
        { MutexLock lock(mutex);
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
                    if (nIntersections == 1 && !indirectDone && shooter->finalGather) {
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
                if (depositedPhoton && shooter->finalGather &&
                        rng->RandomFloat() < .125f) {
                    Normal n = photonIsect.dg.nn;
                    n = Faceforward(n, -photonRay.d);
                    localRadiancePhotons.push_back(RadiancePhoton(photonIsect.dg.p, n));
                    Spectrum rho_r = photonBSDF->rho(*rng, BSDF_ALL_REFLECTION);
                    localRpReflectances.push_back(rho_r);
                    Spectrum rho_t = photonBSDF->rho(*rng, BSDF_ALL_TRANSMISSION);
                    localRpTransmittances.push_back(rho_t);
                }
            }
        }//Mutex lock

        if (nIntersections >= shooter->maxPhotonDepth) 
            return;

        // Sample new photon ray direction

        Vector wi;
        float pdf;
        BxDFType flags;
        for(uint32_t i=0; i<spectrums.size(); ++i){
            alpha = spectrums[i];
            
            Spectrum fr = photonBSDF->Sample_f(wo, &wi, BSDFSample(*rng),
                                               &pdf, BSDF_ALL, &flags, &alpha);

            if (fr.IsBlack() || pdf == 0.f) 
                continue;

            Spectrum anew = alpha * fr *
                AbsDot(wi, photonBSDF->dgShading.nn) / pdf;

            // Possibly terminate photon path with Russian roulette
            float continueProb = min(1.f, anew.y() / alpha.y());
            if (rng->RandomFloat() > continueProb)
                continue;

            alpha = anew / continueProb;
            specularPath &= ((flags & BSDF_SPECULAR) != 0);
            
            if (indirectDone && !specularPath) 
                continue;

            photonRay = RayDifferential(photonIsect.dg.p, wi, photonRay,
                                        photonIsect.rayEpsilon);
            followPhoton(photonRay, photonIsect, alpha, nIntersections, specularPath,
                localDirectPhotons, localIndirectPhotons, localCausticPhotons, localVolumePhotons, localRadiancePhotons,
                causticDone, indirectDone, volumeDone, arena, rng, rpReflectances, rpTransmittances);
        }
}
}


void PhotonShootingTask::Run() {
    // Declare local variables for _PhotonShootingTask_
    MemoryArena arena;
    RNG rng(31 * taskNum);
    vector<Photon> localDirectPhotons, localIndirectPhotons, localCausticPhotons, localVolumePhotons;
    vector<RadiancePhoton> localRadiancePhotons;
    uint32_t totalPaths = 0;
    bool causticDone    = (shooter->nCausticPhotonsWanted == 0);
    bool indirectDone   = (shooter->nIndirectPhotonsWanted == 0);
    bool volumeDone     = (shooter->nVolumePhotonsWanted == 0);

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
                    causticDone, indirectDone, volumeDone, &arena, &rng, rpReflectances, rpTransmittances);
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
            (unsuccessful(shooter->nCausticPhotonsWanted,
                                      causticPhotons.size(), blockSize) ||
             unsuccessful(shooter->nIndirectPhotonsWanted,
                                      indirectPhotons.size(), blockSize)||
             unsuccessful(shooter->nVolumePhotonsWanted,
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
            shooter->nIndirectPaths += blockSize;
            for (uint32_t i = 0; i < localIndirectPhotons.size(); ++i)
                indirectPhotons.push_back(localIndirectPhotons[i]);
            localIndirectPhotons.erase(localIndirectPhotons.begin(),
                                       localIndirectPhotons.end());
            if (indirectPhotons.size() >= shooter->nIndirectPhotonsWanted)
                indirectDone = true;
            nDirectPaths += blockSize;
            for (uint32_t i = 0; i < localDirectPhotons.size(); ++i)
                directPhotons.push_back(localDirectPhotons[i]);
            localDirectPhotons.erase(localDirectPhotons.begin(),
                                     localDirectPhotons.end());
        }

        // Merge direct, caustic, and radiance photons into shared array
        if (!causticDone) {
            shooter->nCausticPaths += blockSize;
            for (uint32_t i = 0; i < localCausticPhotons.size(); ++i){
                causticPhotons.push_back(localCausticPhotons[i]);
            }
            localCausticPhotons.erase(localCausticPhotons.begin(), localCausticPhotons.end());
            if (causticPhotons.size() >= shooter->nCausticPhotonsWanted)
                causticDone = true;
        }

        if(!volumeDone){
            shooter->nVolumePaths += blockSize;
            for(uint32_t i = 0; i < localVolumePhotons.size(); ++i){
                localVolumePhotons[i].alpha /= float(nshot);
                // localVolumePhotons[i].alpha /= float(localVolumePhotons.size());
                volumePhotons.push_back(localVolumePhotons[i]);
            }
            localVolumePhotons.erase(localVolumePhotons.begin(), localVolumePhotons.end());
            if(volumePhotons.size() >= shooter->nVolumePhotonsWanted)
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



PhotonShooter::PhotonShooter(int ncaus, int nind, int nvol, float steps,
        int nl, int mdepth, int mphodepth, float mdist, bool fg,
        int gs, float ga) {

   
    stepSize = steps;
    nCausticPhotonsWanted = ncaus;
    nIndirectPhotonsWanted = nind;
    nVolumePhotonsWanted = nvol;
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


PhotonShooter::~PhotonShooter() {
    delete[] lightSampleOffsets;
    delete[] bsdfSampleOffsets;
    delete causticMap;
    delete indirectMap;
    delete radianceMap;
    delete volumeMap;
}


// void PhotonShooter::RequestSamples(Sampler *sampler, Sample *sample,
//         const Scene *scene) {
//     // Allocate and request samples for sampling all lights
//     uint32_t nLights = scene->lights.size();
//     lightSampleOffsets = new LightSampleOffsets[nLights];
//     bsdfSampleOffsets = new BSDFSampleOffsets[nLights];
//     for (uint32_t i = 0; i < nLights; ++i) {
//         const Light *light = scene->lights[i];
//         int nSamples = light->nSamples;
//         if (sampler) nSamples = sampler->RoundSize(nSamples);
//         lightSampleOffsets[i] = LightSampleOffsets(nSamples, sample);
//         bsdfSampleOffsets[i] = BSDFSampleOffsets(nSamples, sample);
//     }

//     // Request samples for final gathering
//     if (finalGather) {
//         gatherSamples = max(1, gatherSamples/2);
//         if (sampler) gatherSamples = sampler->RoundSize(gatherSamples);
//         bsdfGatherSampleOffsets = BSDFSampleOffsets(gatherSamples, sample);
//         indirGatherSampleOffsets = BSDFSampleOffsets(gatherSamples, sample);
//     }
// }


void PhotonShooter::Preprocess(const Scene *scene,
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


PhotonShooter *CreatePhotonShooter(const ParamSet &surfparams, const ParamSet &volparams){
    int nCaustic = surfparams.FindOneInt("causticphotons", 20000);
    int nIndirect = surfparams.FindOneInt("indirectphotons", 100000);
    int nVolume = volparams.FindOneInt("volumephotons", 0);
    float stepSize = surfparams.FindOneFloat("stepsize", 0.1f);
    int nUsed = surfparams.FindOneInt("nused", 50);
    if (PbrtOptions.quickRender) nCaustic = nCaustic / 10;
    if (PbrtOptions.quickRender) nIndirect = nIndirect / 10;
    if (PbrtOptions.quickRender) nUsed = max(1, nUsed / 10);
    int maxSpecularDepth = surfparams.FindOneInt("maxspeculardepth", 5);
    int maxPhotonDepth = surfparams.FindOneInt("maxphotondepth", 5);
    bool finalGather = surfparams.FindOneBool("finalgather", true);
    int gatherSamples = surfparams.FindOneInt("finalgathersamples", 32);
    if (PbrtOptions.quickRender) gatherSamples = max(1, gatherSamples / 4);
    float maxDist = surfparams.FindOneFloat("maxdist", .1f);
    float gatherAngle = surfparams.FindOneFloat("gatherangle", 10.f);
    return new PhotonShooter(nCaustic, nIndirect, nVolume, stepSize,
        nUsed, maxSpecularDepth, maxPhotonDepth, maxDist, finalGather, gatherSamples,
        gatherAngle);
}