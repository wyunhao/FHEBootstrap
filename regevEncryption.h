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
regevPK regevGeneratePublicKey_Mod3(const regevParam& param, const regevSK& sk);
regevPK regevGenerateSquareRootInput(const regevParam& param, const regevSK& sk);
void regevEncSK(regevCiphertext& ct, const int& msg, const regevSK& sk, const regevParam& param, const bool& pk_gen = false);
void regevEncSK_Mod3(regevCiphertext& ct, const int& msg, const regevSK& sk, const regevParam& param, const int enc_num);
void regevEncSK_Value(regevCiphertext& ct, const int msg, const regevSK& sk, const regevParam& param);
void regevEncPK(regevCiphertext& ct, const int& msg, const regevPK& pk, const regevParam& param);
void regevDec(int& msg, const regevCiphertext& ct, const regevSK& sk, const regevParam& param);
void regevDec_Mod3(int& msg, const regevCiphertext& ct, const regevSK& sk, const regevParam& param);
void regevDec_Value(int& msg, const regevCiphertext& ct, const regevSK& sk, const regevParam& param, const int errorRange);

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

void regevEncSK_Mod3(regevCiphertext& ct, const int& msg, const regevSK& sk, const regevParam& param, const int enc_num){
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
    if (enc_num) {
        ct.b.ModAddFastEq(q/3 * enc_num, q);
    }
    DiscreteGaussianGeneratorImpl<NativeVector> m_dgg(param.std_dev);
    ct.b.ModAddFastEq(m_dgg.GenerateInteger(q), q);
}


regevPK regevGeneratePublicKey(const regevParam& param, const regevSK& sk, const bool& add_half = false, const bool& enc_one = false){
    regevPK pk(param.m);
    for(int i = 0; i < param.m; i++){
        if (enc_one) { // all of the ciphertexts will encrypt one
            regevEncSK(pk[i], 0, sk, param, true, true);
        } else if (add_half && i%2 == 0) { // half of the encrypted ct will add q/4, which means it encrypts one
            regevEncSK(pk[i], 0, sk, param, true, true);
        } else {
            regevEncSK(pk[i], 0, sk, param, true, false);
        }
    }
    return pk;
}


regevPK regevGeneratePublicKey_Mod3(const regevParam& param, const regevSK& sk, const int enc_num = 0){
    regevPK pk(param.m);
    for(int i = 0; i < param.m; i++){
        regevEncSK_Mod3(pk[i], 0, sk, param, enc_num);
    }
    return pk;
}

void regevEncSK_Value(regevCiphertext& ct, const int msg, const regevSK& sk, const regevParam& param, const int errorRange = 128){
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
    if (msg) {
        ct.b.ModAddFastEq(errorRange * msg, q);
    }
    DiscreteGaussianGeneratorImpl<NativeVector> m_dgg(param.std_dev);
    ct.b.ModAddFastEq(m_dgg.GenerateInteger(q), q);
}

// map 2^9 to 2^16
regevPK regevGenerateSquareRootInput(const regevParam& param, const regevSK& sk, const int plaintextSpace = 512, const int errorRange = 128){
    regevPK pk(param.m);
    for(int i = 0; i < param.m; i++){
        int val = i % plaintextSpace; // the value to encrypt
        regevEncSK_Value(pk[i], val, sk, param, errorRange);
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

void regevDec(vector<int>& msg, const vector<regevCiphertext>& ct, const regevSK& sk, const regevParam& param){
    msg.resize(ct.size());

    int q = param.q;
    int n = param.n;
    NativeInteger inner(0);
    for (int j = 0; j < (int) ct.size(); j++) {
        int r = ct[j].b.ConvertToInt();
        for (int i = 0; i < n; ++i) {
            r = (r - ct[j].a[i].ConvertToInt() * sk[i].ConvertToInt()) % q;
        }
        msg[j] = (r >= 0 && r < 65537/4) || (r < 65537 && r > 65537-65537/4) ? 0 : 1;
    }
}

// map 2^16 to 2^9
void regevDec_Value(vector<int>& msg, const vector<regevCiphertext>& ct, const regevSK& sk, const regevParam& param, const int errorRange){
    msg.resize(ct.size());

    int q = param.q;
    int n = param.n;
    NativeInteger inner(0);
    for (int i = 0; i < (int) ct.size(); i++) {
        int temp = 0;
        for (int j = 0; j < n; j++) {
            long mul_tmp = (ct[i].a[j].ConvertToInt() * sk[j].ConvertToInt()) % q;
            mul_tmp = mul_tmp < 0 ? mul_tmp + q : mul_tmp;
            temp = (temp + (int) mul_tmp) % q;
        }
        temp = (ct[i].b.ConvertToInt() + temp) % q;

        temp = (temp + errorRange/2) % q; // 64 is error bound for 2^9
        msg[i] = temp / errorRange;
    }
    cout << endl;
}


void regevDec_Mod3(vector<int>& msg, const vector<regevCiphertext>& ct, const regevSK& sk, const regevParam& param){
    msg.resize(ct.size());

    int q = param.q;
    int n = param.n;
    NativeInteger inner(0);
    for (int i = 0; i < (int) ct.size(); i++) {
        int temp = 0;
        for (int j = 0; j < n; j++) {
            long mul_tmp = (ct[i].a[j].ConvertToInt() * sk[j].ConvertToInt()) % q;
            mul_tmp = mul_tmp < 0 ? mul_tmp + q : mul_tmp;
            temp = (temp + (int) mul_tmp) % q;
        }
        temp = (temp + ct[i].b.ConvertToInt()) % q;
        if ((temp >= 0 && temp < q/6) || (temp < q && temp > q-q/6)) {
            msg[i] = 0;
        } else if ((temp >= q/6 && temp < q/2)) {
            msg[i] = 1;
        } else {
            msg[i] = 2;
        }
    }
}

bool shouldNegate(const int i, const int dim = 32768) {
    return!((i < dim/6) || (i < 3*dim/6 && i >= 2*dim/6) || (i < 5*dim/6 && i >= 4*dim/6));
}

void regevDec_Mod3_Mixed(vector<int>& msg, const vector<regevCiphertext>& ct, const regevSK& sk, const regevParam& param){
    msg.resize(ct.size());

    int q = param.q;
    int n = param.n;
    NativeInteger inner(0);
    for (int i = 0; i < (int) ct.size(); i++) {
        int temp = 0;
        for (int j = 0; j < n; j++) {
            long mul_tmp = (ct[i].a[j].ConvertToInt() * sk[j].ConvertToInt()) % q;
            mul_tmp = mul_tmp < 0 ? mul_tmp + q : mul_tmp;
            temp = (temp + (int) mul_tmp) % q;
        }
        temp = (temp + ct[i].b.ConvertToInt()) % q;
        if ((temp >= 0 && temp < q/6) || (temp < q && temp > q-q/6)) {
            msg[i] = shouldNegate(i) ? 1 : 0 ;
        } else if ((temp >= q/6 && temp < q/2)) {
            msg[i] = shouldNegate(i) ? 0 : 1 ;
        } else {
            msg[i] = 2;
        }
    }
}
