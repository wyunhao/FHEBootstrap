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
    int ring_dim = 8;
    int n = 4;
    int p = 65537;

    EncryptionParameters bfv_params(scheme_type::bfv);
    bfv_params.set_poly_modulus_degree(ring_dim);
    auto coeff_modulus = CoeffModulus::Create(ring_dim, { 28, 60, 60 });
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

    KeyGenerator keygen(seal_context);
    SecretKey bfv_secret_key = keygen.secret_key();

    MemoryPoolHandle my_pool = MemoryPoolHandle::New();

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

    // vector<uint64_t> msg = {21845, 21845, 32768, 32768, 0, 0, 43691, 43691};
    vector<uint64_t> msg = {0, 21845, 32768, 43490, 10922, 30000, 50000, 20000};

    Plaintext pl;
    Ciphertext c;
    batch_encoder.encode(msg, pl);

    for (int i = 0; i < ring_dim; i++) {
        cout << pl[i] << " ";
    }
    cout << endl;
    encryptor.encrypt(pl, c);


    Ciphertext c_copy(c);

    int sq_ct = sqrt(ring_dim/2);
    vector<Ciphertext> ct_sqrt_list(2*sq_ct);

    evaluator.rotate_columns_inplace(c_copy, gal_keys);
    for (int i = 0; i < sq_ct; i++) {
        evaluator.rotate_rows(c, sq_ct * i, gal_keys, ct_sqrt_list[i]);
        evaluator.transform_to_ntt_inplace(ct_sqrt_list[i]);
        evaluator.rotate_rows(c_copy, sq_ct * i, gal_keys, ct_sqrt_list[i+sq_ct]);
        evaluator.transform_to_ntt_inplace(ct_sqrt_list[i+sq_ct]);
    }

    // evaluator.rotate_rows_inplace(c, 1, gal_keys);

    // decryptor.decrypt(c, pl);
    // batch_encoder.decode(pl, msg);
    // cout << "Decode: " << msg << endl;


    // for (int i = 0; i < ct_sqrt_list.size(); i++) {
    //     Plaintext pp;
    //     vector<uint64_t> v(ring_dim);

    //     evaluator.transform_from_ntt_inplace(ct_sqrt_list[i]);

    //     decryptor.decrypt(ct_sqrt_list[i], pp);
    //     batch_encoder.decode(pp, v);
    //     cout << v << endl;

    // }


    Ciphertext coeff = slotToCoeff_WOPrepreocess(seal_context, ct_sqrt_list, gal_keys, ring_dim);

    // evaluator.rotate_columns_inplace(c, gal_keys);

    decryptor.decrypt(coeff, pl);
    for (int i = 0; i < ring_dim; i++) {
        cout << pl[i] << " ";
    }
    cout << endl;

}