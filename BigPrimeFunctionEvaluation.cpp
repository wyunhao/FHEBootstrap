#include "regevEncryption.h"
#include "global.h"
#include "util.h"
#include "seal/seal.h"
#include "seal/util/iterator.h"
#include <numeric>
#include <stdio.h>

using namespace std;
using namespace seal;


// for square root
int main() {


    ////////////////////////////////////////////// PREPARE (R)LWE PARAMS ///////////////////////////////////////////////
    int ring_dim = poly_modulus_degree_glb;
    int n = 1024;
    BootstrapParam bootstrap_param = BootstrapParam(786433, 192, 4096, 256*3, 1024);
    int p = bootstrap_param.ciphertextSpacePrime;

    EncryptionParameters bfv_params(scheme_type::bfv);
    bfv_params.set_poly_modulus_degree(ring_dim);

    auto coeff_modulus = CoeffModulus::Create(ring_dim, { 60, 50, 50, 60, 60,
                                                          60, 60, 60, 60, 60,
                                                          60, 60, 60, 60, 60 });
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

    GaloisKeys gal_keys, gal_keys_coeff;
    vector<int> rot_steps = {1};
    for (int i = 0; i < n;) {
        rot_steps.push_back(i);
        i += sqrt(n);
    }
    keygen.create_galois_keys(rot_steps, gal_keys);
    
    vector<Modulus> coeff_modulus_last = coeff_modulus;
    coeff_modulus_last.erase(coeff_modulus_last.begin() + 2, coeff_modulus_last.end()-1);
    EncryptionParameters parms_last = bfv_params;
    parms_last.set_coeff_modulus(coeff_modulus_last);
    SEALContext seal_context_last = SEALContext(parms_last, true, sec_level_type::none);

    SecretKey sk_last;
    sk_last.data().resize(coeff_modulus_last.size() * ring_dim);
    sk_last.parms_id() = seal_context_last.key_parms_id();
    util::set_poly(bfv_secret_key.data().data(), ring_dim, coeff_modulus_last.size() - 1, sk_last.data().data());
    util::set_poly(
        bfv_secret_key.data().data() + ring_dim * (coeff_modulus.size() - 1), ring_dim, 1,
        sk_last.data().data() + ring_dim * (coeff_modulus_last.size() - 1));

    vector<int> rot_steps_coeff = {1};
    for (int i = 0; i < ring_dim/2;) {
        if (find(rot_steps_coeff.begin(), rot_steps_coeff.end(), i) == rot_steps_coeff.end()) {
            rot_steps_coeff.push_back(i);
        }
        i += sqrt(ring_dim/2);
    }
    KeyGenerator keygen_last(seal_context_last, sk_last);
    keygen_last.create_galois_keys(rot_steps_coeff, gal_keys_coeff);


    ////////////////////////////////////////////// PREPARE LWE CIPHERTEXT //////////////////////////////////////////////
    
    // 32*32 = 1024
    // convert the BFV ternary sk into the new LWE key we should eventually switch to
    auto lwe_params = regevParam(n, p, 1.3, ring_dim); 
    auto lwe_sk = regevGenerateSecretKey(lwe_params);
    for (int i = 0; i < n; i++) {
        lwe_sk[i] = (int) new_key.data()[i] > p ? p-1 : new_key.data()[i];
    }
    cout << endl;

    seal::util::RNSIter new_key_rns(new_key.data().data(), ring_dim);
    ntt_negacyclic_harvey(new_key_rns, coeff_modulus.size(), seal_context.key_context_data()->small_ntt_tables());

    vector<int> msg(ring_dim);
    vector<regevCiphertext> lwe_ct_list = regevGenerateSquareRootInput(lwe_params, lwe_sk, bootstrap_param.plaintextSpace, bootstrap_param.errorRange); // enc 0


    ////////////////////////////////////////////// ENCRYPT SK UNDER BFV ////////////////////////////////////////////////

    // one switching key for one lwe_sk
    Ciphertext lwe_sk_encrypted = encryptLWEskUnderBFV(seal_context, ring_dim, bfv_public_key, bfv_secret_key, lwe_sk, lwe_params);


    /////////////////////////////////////////////////// BOOTSTRAP //////////////////////////////////////////////////////
    vector<uint64_t> q_shift_constant(ring_dim, 0);
    vector<regevCiphertext> lwe_ct_results = bootstrap_bigPrime(lwe_ct_list, lwe_sk_encrypted, seal_context, seal_context_last, relin_keys, gal_keys,
                                                                gal_keys_coeff, ring_dim, n, p, ksk, rangeCheckIndices_squareRoot_20, my_pool, bfv_secret_key,
                                                                q_shift_constant, 0, false, false, bootstrap_param.firstLevelDegree, bootstrap_param.secondLevelDegree,
                                                                1152921504581419009);


    regevDec_Value(msg, lwe_ct_results, lwe_sk, lwe_params, bootstrap_param.errorRange);

    cout << "Actual result: \n" << msg << endl;
}