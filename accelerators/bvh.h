 
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

#ifndef PBRT_ACCELERATORS_BVH_H
#define PBRT_ACCELERATORS_BVH_H

// accelerators/bvh.h*
#include "pbrt.h"
#include "primitive.h"
#include <list>

using std::list;
struct BVHBuildNode;

// BVHAccel Forward Declarations
struct BVHPrimitiveInfo;
struct LinearBVHNode; 

// BVHAccel Declarations
class BVHAccel : public Aggregate {
public:
    // BVHAccel Public Methods
    
    BVHAccel(const vector<Reference<Primitive> > &p, uint32_t maxPrims = 1,
             const string &sm = "sah");
    BBox WorldBound() const;
    bool CanIntersect() const { return true; }
    ~BVHAccel();
    bool Intersect(const Ray &ray, Intersection *isect) const;
    bool IntersectP(const Ray &ray) const;
private:
    // BVHAccel Private Methods
    BVHBuildNode* buildAAC(MemoryArena &buildArena,
        vector<BVHPrimitiveInfo> &buildData, uint32_t *totalNodes,
        vector<Reference<Primitive> > &orderedPrims);
    
    vector<BVHBuildNode*> recursiveBuildAAC(MemoryArena &buildArena,
    						vector<BVHPrimitiveInfo>& buildData,
                            uint32_t start, uint32_t end,
    						uint32_t* totalNodes, uint32_t bit,
                            vector<Reference<Primitive> >& orderedPrims); 

    vector<BVHBuildNode*> combineClusters(MemoryArena &buildArena,
    						vector<BVHBuildNode*> &clusters,
                            uint32_t n, uint32_t* totalNodes);

    BVHBuildNode *recursiveBuild(MemoryArena &buildArena,
        vector<BVHPrimitiveInfo> &buildData, uint32_t start, uint32_t end,
        uint32_t *totalNodes, vector<Reference<Primitive> > &orderedPrims);
    
    uint32_t flattenBVHTree(BVHBuildNode *node, uint32_t *offset);

    inline uint32_t numberOfCluster(int x){return c*pow(x,alpha);}
    void initAAC(){
        delta = 4;
        epsilon = 0.2;
        c = 0.5*pow(delta,0.5+epsilon);
        alpha = 0.5-epsilon;
    };
    // BVHAccel Private Data
    uint32_t maxPrimsInNode;
    enum SplitMethod { SPLIT_MIDDLE, SPLIT_EQUAL_COUNTS, SPLIT_SAH, SPLIT_AAC};
    SplitMethod splitMethod;
    vector<Reference<Primitive> > primitives;
    LinearBVHNode *nodes;

    //AAC variables
    uint32_t delta;
    float epsilon;
    float c;
    float alpha;
};


BVHAccel *CreateBVHAccelerator(const vector<Reference<Primitive> > &prims,
        const ParamSet &ps);

#endif // PBRT_ACCELERATORS_BVH_H
