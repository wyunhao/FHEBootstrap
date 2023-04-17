#include "regevEncryption.h"
#include "global.h"
#include "util.h"
#include "seal/seal.h"
#include "seal/util/iterator.h"
#include <numeric>
#include <stdio.h>

using namespace seal;
using namespace std;


int main() {


    ////////////////////////////////////////////// PREPARE (R)LWE PARAMS ///////////////////////////////////////////////
    int ring_dim = poly_modulus_degree_glb;
    int n = 1024;
    int p = 65537;

    EncryptionParameters bfv_params(scheme_type::bfv);
    bfv_params.set_poly_modulus_degree(ring_dim);
    auto coeff_modulus = CoeffModulus::Create(ring_dim, { 28, 60, 55, 60, 60,
                                                          60, 60, 60, 60, 60,
                                                          50, 60 });
    bfv_params.set_coeff_modulus(coeff_modulus);
    bfv_params.set_plain_modulus(p);

    prng_seed_type seed;
    for (auto &i : seed) {
        i = random_uint64();
    }
    auto rng = make_shared<Blake2xbPRNGFactory>(Blake2xbPRNGFactory(seed));
    bfv_params.set_random_generator(rng);

    SEALContext seal_context(bfv_params, true, sec_level_type::none);
    cout << "primitive root: " << seal_context.first_context_data()->plain_ntt_tables()->get_root() << endl;

    KeyGenerator new_key_keygen(seal_context, n);
    SecretKey new_key = new_key_keygen.secret_key();
    inverse_ntt_negacyclic_harvey(new_key.data().data(), seal_context.key_context_data()->small_ntt_tables()[0]);

    KeyGenerator keygen(seal_context);
    SecretKey bfv_secret_key = keygen.secret_key();

    MemoryPoolHandle my_pool = MemoryPoolHandle::New();

    // generate a key switching key based on key_before and secret_key
    KSwitchKeys ksk;
    seal::util::ConstPolyIter secret_key_before(bfv_secret_key.data().data(), ring_dim, coeff_modulus.size());

    new_key_keygen.generate_kswitch_keys(secret_key_before, 1, static_cast<KSwitchKeys &>(ksk), false);
    ksk.parms_id() = seal_context.key_parms_id();

    PublicKey bfv_public_key;
    keygen.create_public_key(bfv_public_key);

    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(seal_context, bfv_public_key);
    Evaluator evaluator(seal_context);
    BatchEncoder batch_encoder(seal_context);
    Decryptor decryptor(seal_context, bfv_secret_key);

    GaloisKeys gal_keys;
    vector<int> rot_steps = {1};
    for (int i = 0; i < n;) {
        rot_steps.push_back(i);
        i += sqrt(n);
    }
    for (int i = 0; i < ring_dim/2;) {
        if (find(rot_steps.begin(), rot_steps.end(), i) == rot_steps.end()) {
            rot_steps.push_back(i);
        }
        i += sqrt(ring_dim/2);
    }
    keygen.create_galois_keys(rot_steps, gal_keys);


    ////////////////////////////////////////////// PREPARE LWE CIPHERTEXT //////////////////////////////////////////////
    
    // 32*32 = 1024
    // convert the BFV ternary sk into the new LWE key we should eventually switch to
    auto lwe_params = regevParam(n, p, 1.3, ring_dim); 
    auto lwe_sk = regevGenerateSecretKey(lwe_params);
    for (int i = 0; i < n; i++) {
        lwe_sk[i] = (int) new_key.data()[i] > p ? p-1 : new_key.data()[i];
    }

    seal::util::RNSIter new_key_rns(new_key.data().data(), ring_dim);
    ntt_negacyclic_harvey(new_key_rns, coeff_modulus.size(), seal_context.key_context_data()->small_ntt_tables());

    vector<int> msg(ring_dim);

    vector<regevCiphertext> lwe_ct_list_1 = regevGeneratePublicKey_Mod3(lwe_params, lwe_sk, 0); // enc 0
    vector<regevCiphertext> lwe_ct_list_2 = regevGeneratePublicKey_Mod3(lwe_params, lwe_sk, 1); // enc 1

    // regevDec_Mod3(msg, lwe_ct_list_1, lwe_sk, lwe_params);
    // cout << "Input Ciphertex 1: \n" << msg << endl;
    // regevDec_Mod3(msg, lwe_ct_list_2, lwe_sk, lwe_params);
    // cout << "Input Ciphertex 2: \n" << msg << endl;

    vector<regevCiphertext> lwe_ct_list = preprocess_NAND(lwe_ct_list_1, lwe_ct_list_2, lwe_params);

    // regevDec_Mod3(msg, lwe_ct_list, lwe_sk, lwe_params);
    // cout << "Expected NAND result: \n" << msg << endl;

    ////////////////////////////////////////////// ENCRYPT SK UNDER BFV ////////////////////////////////////////////////

    // one switching key for one lwe_sk
    Ciphertext lwe_sk_encrypted = encryptLWEskUnderBFV(seal_context, ring_dim, bfv_public_key, bfv_secret_key, lwe_sk, lwe_params);


    /////////////////////////////////////////////////// BOOTSTRAP //////////////////////////////////////////////////////
    bool gateEval = true;
    int f_zero = p/3;


    vector<uint64_t> q_shift_constant_NAND(ring_dim, -p/6);
    vector<regevCiphertext> lwe_ct_results = bootstrap(lwe_ct_list, lwe_sk_encrypted, seal_context, relin_keys, gal_keys,
                                                       ring_dim, n, p, ksk, rangeCheckIndices_gateEvaluation, my_pool, bfv_secret_key, q_shift_constant_NAND, f_zero, gateEval);

    regevDec_Mod3(msg, lwe_ct_results, lwe_sk, lwe_params);
    cout << "Actual NAND result: \n" << msg << endl;





    vector<uint64_t> q_shift_constant_OR(ring_dim, -p/6-p/3);
    lwe_ct_results = bootstrap(lwe_ct_list, lwe_sk_encrypted, seal_context, relin_keys, gal_keys,
                                                       ring_dim, n, p, ksk, rangeCheckIndices_gateEvaluation, my_pool, bfv_secret_key, q_shift_constant_OR, f_zero, gateEval);

    regevDec_Mod3(msg, lwe_ct_results, lwe_sk, lwe_params);
    cout << "Actual OR result: \n" << msg << endl;




    vector<uint64_t> q_shift_constant_XNOR(ring_dim, -p/6+p/3);
    lwe_ct_results = bootstrap(lwe_ct_list, lwe_sk_encrypted, seal_context, relin_keys, gal_keys,
                                                       ring_dim, n, p, ksk, rangeCheckIndices_gateEvaluation, my_pool, bfv_secret_key, q_shift_constant_XNOR, f_zero, gateEval);

    regevDec_Mod3(msg, lwe_ct_results, lwe_sk, lwe_params);
    cout << "Actual XNOR result: \n" << msg << endl;
}