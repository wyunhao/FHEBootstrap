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















    auto lwe_params = regevParam(n, p, 1.3, ring_dim); 
    auto lwe_sk = regevGenerateSecretKey(lwe_params);
    Ciphertext lwe_sk_encrypted = encryptLWEskUnderBFV(seal_context, ring_dim, bfv_public_key, bfv_secret_key, lwe_sk, lwe_params);

    vector<int> msg;

    vector<regevCiphertext> lwe_ct_list = regevGeneratePublicKey_Mod3(lwe_params, lwe_sk, 2); // enc 1


    int sq_sk = sqrt(n), sq_ct = sqrt(ring_dim/2);
    vector<Ciphertext> lwe_sk_sqrt_list(sq_sk), ct_sqrt_list(2*sq_ct);

    Ciphertext lwe_sk_column;

    evaluator.rotate_columns(lwe_sk_encrypted, gal_keys, lwe_sk_column);
    for (int i = 0; i < sq_sk; i++) {
        evaluator.rotate_rows(lwe_sk_encrypted, sq_sk * i, gal_keys, lwe_sk_sqrt_list[i]);
        evaluator.transform_to_ntt_inplace(lwe_sk_sqrt_list[i]);
    }






    Ciphertext result = evaluatePackedLWECiphertext(seal_context, lwe_ct_list, lwe_sk_sqrt_list, gal_keys, n, ring_dim);

    Ciphertext range_check_res;
    // f(0) = 1
    Bootstrap_RangeCheck_PatersonStockmeyer(range_check_res, result, rangeCheckIndices_gateEvaluation, p, ring_dim, relin_keys, seal_context, 65537/3);




    Plaintext pl;
    vector<uint64_t> v(ring_dim);

    decryptor.decrypt(range_check_res, pl);
    batch_encoder.decode(pl, v);

    cout << "Decrypted after rangeCheck:\n";
    cout << v << endl;
    




}