#pragma once

#include "math/ternaryuniformgenerator.h"
#include "math/discreteuniformgenerator.h"
#include "math/discretegaussiangenerator.h"
#include "seal/seal.h"
#include <iostream>
#include <random>

using namespace std;
using namespace lbcrypto;
using namespace seal;

struct regevParam{
    int n;
    int q;
    double std_dev;
    int m;
    regevParam(){
        n = 450;
        q = 65537;
        std_dev = 1.3;
        m = 16000; 
    }
    regevParam(int n, int q, double std_dev, int m)
    : n(n), q(q), std_dev(std_dev), m(m)
    {}
};

typedef NativeVector regevSK;

struct regevCiphertext{
    NativeVector a;
    NativeInteger b;
};

typedef vector<regevCiphertext> regevPK;

regevSK regevGenerateSecretKey(const regevParam& param);
regevPK regevGeneratePublicKey(const regevParam& param, const regevSK& sk);
void regevEncSK(regevCiphertext& ct, const int& msg, const regevSK& sk, const regevParam& param, const bool& pk_gen = false);
void regevEncPK(regevCiphertext& ct, const int& msg, const regevPK& pk, const regevParam& param);

/////////////////////////////////////////////////////////////////// Below are implementation

regevSK regevGenerateSecretKey(const regevParam& param){
    int n = param.n;
    int q = param.q;
    lbcrypto::TernaryUniformGeneratorImpl<regevSK> tug;
    return tug.GenerateVector(n, q);
}

void regevEncSK(regevCiphertext& ct, const int& msg, const regevSK& sk, const regevParam& param, const bool& pk_gen, const bool& add_half = false){
    NativeInteger q = param.q;
    int n = param.n;
    DiscreteUniformGeneratorImpl<NativeVector> dug;
    dug.SetModulus(q);
    ct.a = dug.GenerateVector(n);
    NativeInteger mu = q.ComputeMu();
    for (int i = 0; i < n; ++i) {
        ct.b += ct.a[i].ModMulFast(sk[i], q, mu);
    }
    ct.b.ModEq(q);
    if (add_half) {
        ct.b.ModAddFastEq(q/2, q);
    }
    if(!pk_gen)
        msg? ct.b.ModAddFastEq(3*q/4, q) : ct.b.ModAddFastEq(q/4, q);
    DiscreteGaussianGeneratorImpl<NativeVector> m_dgg(param.std_dev);
    ct.b.ModAddFastEq(m_dgg.GenerateInteger(q), q);
}

regevPK regevGeneratePublicKey(const regevParam& param, const regevSK& sk, const bool& add_half = false){
    regevPK pk(param.m);
    for(int i = 0; i < param.m; i++){
        if (add_half && i%2 == 0) {
            regevEncSK(pk[i], 0, sk, param, true, true);
        } else {
            regevEncSK(pk[i], 0, sk, param, true, false);
        }
    }
    return pk;
}

void regevEncPK(regevCiphertext& ct, const int& msg, const regevPK& pk, const regevParam& param){
    prng_seed_type seed;
    for (auto &i : seed) {
        i = random_uint64();
    }

    auto rng = make_shared<Blake2xbPRNGFactory>(Blake2xbPRNGFactory(seed));
    RandomToStandardAdapter engine(rng->create());
    uniform_int_distribution<uint64_t> dist(0, 1);

    NativeInteger q = param.q;
    ct.a = NativeVector(param.n);
    for(size_t i = 0; i < pk.size(); i++){
        if (dist(engine)){
            for(int j = 0; j < param.n; j++){
                ct.a[j].ModAddFastEq(pk[i].a[j], q);
            }
            ct.b.ModAddFastEq(pk[i].b, q);
        }
    }
    msg? ct.b.ModAddFastEq(3*q/4, q) : ct.b.ModAddFastEq(q/4, q);
}

