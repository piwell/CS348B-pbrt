#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_PHOTONSHOOTER_H
#define PBRT_PHOTONSHOOTER_H

// core/photonshooter.h*
#include "pbrt.h"
#include "integrator.h"
#include "kdtree.h"
#include "volume.h"
#include "scene.h"
#include "montecarlo.h"
// #include "integrators/photonmap.h"

class PhotonIntegrator;
class PhotonVolumeIntegrator;

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


class PhotonShooter{
      // PhotonIntegrator Public Methods
public: 
    PhotonShooter(int ncaus, int nindir, int nvol, float steps, int nLookup, int maxspecdepth,
        int maxphotondepth, float maxdist, bool finalGather, int gatherSamples,
        float ga);
    ~PhotonShooter();
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    void Preprocess(const Scene *scene, const Camera *camera, const Renderer *renderer);

private:
    // PhotonIntegrator Private Methods
    friend class PhotonShootingTask;
    friend class PhotonIntegrator;
    friend class PhotonVolumeIntegrator;

    // PhotonIntegrator Private Data
    uint32_t nCausticPhotonsWanted, nIndirectPhotonsWanted, nVolumePhotonsWanted;
    uint32_t nLookup;
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
    KdTree<Photon> *volumeMap;
    KdTree<RadiancePhoton> *radianceMap;  
    float stepSize;
};


class PhotonShootingTask : public Task {
public:
    PhotonShootingTask(int tn, float ti, Mutex &m, PhotonShooter *in,
        ProgressReporter &prog, bool &at, int &ndp,
        vector<Photon> &direct, vector<Photon> &indir, vector<Photon> &caustic, vector<Photon> &volume,
        vector<RadiancePhoton> &rps, vector<Spectrum> &rpR, vector<Spectrum> &rpT,
        uint32_t &ns, Distribution1D *distrib, const Scene *sc,
        const Renderer *sr)
    : taskNum(tn), time(ti), mutex(m), shooter(in), progress(prog),
      abortTasks(at), nDirectPaths(ndp),
      directPhotons(direct), indirectPhotons(indir), causticPhotons(caustic),volumePhotons(volume),
      radiancePhotons(rps), rpReflectances(rpR), rpTransmittances(rpT),
      nshot(ns), lightDistribution(distrib), scene(sc), renderer (sr) { }
    void Run();

    void followPhoton(RayDifferential photonRay, Intersection photonIsect, Spectrum alpha, int nIntersections, bool specularPath,
        vector<Photon>& localDirectPhotons, vector<Photon>& localIndirectPhotons, vector<Photon>& localCausticPhotons,
        vector<Photon>& localVolumePhotons, vector<RadiancePhoton>& localRadiancePhotons,
        bool& causticDone, bool& indirectDone, bool& volumeDone,
        MemoryArena* arena, RNG* rng,
        vector<Spectrum>& localRpReflectances, vector<Spectrum>& localRpTransmittances);

    int taskNum;
    float time;
    Mutex &mutex;
    PhotonShooter *shooter;
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


PhotonShooter *CreatePhotonShooter(const ParamSet &params);
#endif