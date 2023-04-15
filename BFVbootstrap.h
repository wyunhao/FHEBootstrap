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

    /////////////////////////////////////////////// PREPARE BFV PARAMS /////////////////////////////////////////////////
    int ring_dim = poly_modulus_degree_glb;
    int n = 8;
    int p = 8388449;
    // int p = 65537;
    int sparse_count = 8;

    EncryptionParameters bfv_params(scheme_type::bfv);
    bfv_params.set_poly_modulus_degree(ring_dim);
    auto coeff_modulus = CoeffModulus::Create(ring_dim, { 28, 60, 60});
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


    KeyGenerator short_key_keygen(seal_context, n, sparse_count);
    SecretKey short_key = short_key_keygen.secret_key();

    KeyGenerator keygen(seal_context, ring_dim, sparse_count);
    SecretKey bfv_secret_key = keygen.secret_key();

    KSwitchKeys ksk;
    seal::util::ConstPolyIter secret_key_before(bfv_secret_key.data().data(), ring_dim, coeff_modulus.size());

    short_key_keygen.generate_kswitch_keys(secret_key_before, 1, static_cast<KSwitchKeys &>(ksk), false);
    ksk.parms_id() = seal_context.key_parms_id();

    PublicKey bfv_public_key;
    keygen.create_public_key(bfv_public_key);

    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(seal_context, bfv_public_key);
    encryptor.set_secret_key(bfv_secret_key);
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

    MemoryPoolHandle my_pool = MemoryPoolHandle::New();


    ///////////////////////////////////////// PREPARE ORIGINAL BFV CIPHERTEXT //////////////////////////////////////////

    vector<uint64_t> values = {0,524270,1048540,1572810,2097080,2621350,3145620,3669890,4194161,4718431,5242701,5766971,6291241,6815511,7339781,7864051};
    // for (int i = 0; i < ring_dim; i++) {
    //     // values[i] = (i % 2 == 0);
    //     values[i] = i * 65537 / ring_dim;
    // }
    // cout << "values: " << values << endl;
    Plaintext pl;
    Ciphertext c;
    batch_encoder.encode(values, pl);
    for (int i = 0; i < ring_dim; i++) {
        // pl[i] = i * (8388449/65537) * (65537/ring_dim); 
        // pl[i] = (int) ((float(i) * 65537.0 * (8388449.0/65537.0) / float(ring_dim)) - 7.95*(float)i);
        pl[i] =i*128;
        // pl[i] = values[i];
        cout << pl[i] << " ";
    }
    encryptor.encrypt(pl, c);
    cout << endl << endl;


    inverse_ntt_negacyclic_harvey(bfv_secret_key.data().data(), seal_context.key_context_data()->small_ntt_tables()[0]);
    cout << ">>>>>>>>>>>>>>>>>>> BFV SECRET KEY\n";
    for (int i = 0; i < ring_dim; i++) {
        cout << bfv_secret_key.data()[i] << " ";
    }
    cout << endl << endl;
    seal::util::RNSIter bfv_key_rns(bfv_secret_key.data().data(), ring_dim);
    ntt_negacyclic_harvey(bfv_key_rns, coeff_modulus.size(), seal_context.key_context_data()->small_ntt_tables());

    while(seal_context.last_parms_id() != c.parms_id()){
        evaluator.mod_switch_to_next_inplace(c);
    }

    cout << ">>>>>>>>>>>>>>>>>>> CIPHERTEXT a\n";
    for (int i = 0; i < ring_dim; i++) {
        cout << c.data(1)[i] << " ";
    }
    cout << endl << endl;


    cout << ">>>>>>>>>>>>>>>>>>> CIPHERTEXT b\n";
    for (int i = 0; i < ring_dim; i++) {
        cout << c.data(0)[i] << " ";
    }
    cout << endl << endl;

    decryptor.decrypt(c, pl);

    cout << "Decrypted plaintext: " << endl;
    for (int i = 0; i < ring_dim; i++) {
        cout << pl[i] << " ";
    }
    cout << endl << endl;

    




















    Ciphertext input_ct;
    encryptor.encrypt(pl, input_ct);
    Ciphertext input_ct_copy(input_ct);

    int sq_rt = sqrt(ring_dim/2);
    vector<Ciphertext> input_ct_sqrt_list(2*sq_rt);
    evaluator.rotate_columns_inplace(input_ct_copy, gal_keys);
    for (int i = 0; i < sq_rt; i++) {
        evaluator.rotate_rows(input_ct, sq_rt * i, gal_keys, input_ct_sqrt_list[i]);
        evaluator.transform_to_ntt_inplace(input_ct_sqrt_list[i]);
        evaluator.rotate_rows(input_ct_copy, sq_rt * i, gal_keys, input_ct_sqrt_list[i+sq_rt]);
        evaluator.transform_to_ntt_inplace(input_ct_sqrt_list[i+sq_rt]);
    }


    /////////////////////////////////////////////////// SLOT TO COEFF //////////////////////////////////////////////////
    vector<Plaintext> U_plain_list(ring_dim);
    for (int iter = 0; iter < sq_rt; iter++) {
        for (int j = 0; j < (int) input_ct_sqrt_list.size(); j++) {
            vector<uint64_t> U_tmp = readUtemp(j*sq_rt + iter);
            batch_encoder.encode(U_tmp, U_plain_list[iter * input_ct_sqrt_list.size() + j]);
            evaluator.transform_to_ntt_inplace(U_plain_list[iter * input_ct_sqrt_list.size() + j], input_ct_sqrt_list[j].parms_id());
        }
    }

    Ciphertext coeff = slotToCoeff(seal_context, input_ct_sqrt_list, U_plain_list, gal_keys, ring_dim);

    while(seal_context.last_parms_id() != coeff.parms_id()){
        evaluator.mod_switch_to_next_inplace(coeff);
    }

    /////////////////////////////////////////////////// KEY SWITCHING //////////////////////////////////////////////////
    Ciphertext copy_coeff = coeff;
    auto ct_in_iter = util::iter(copy_coeff);
    ct_in_iter += coeff.size() - 1;
    // notice that the coeff_mod.size() is hardcoded to 1, thus this needs to be performed on the last level
    seal::util::set_zero_poly(ring_dim, 1, coeff.data(1));

    evaluator.switch_key_inplace(coeff, *ct_in_iter, static_cast<const KSwitchKeys &>(ksk), 0, my_pool);

    ///////////////////////////////////////////////////// EXTRACTION ///////////////////////////////////////////////////
    inverse_ntt_negacyclic_harvey(short_key.data().data(), seal_context.key_context_data()->small_ntt_tables()[0]);
    auto lwe_params = regevParam(n, p, 1.3, ring_dim); 
    auto lwe_sk = regevGenerateSecretKey(lwe_params);
    cout << "***************************************** Ternary SK *****************************************\n";
    for (int i = 0; i < n; i++) {
        lwe_sk[i] = short_key.data()[i] > 65537 ? 268369920 : short_key.data()[i];
        cout << lwe_sk[i] << " ";
    }
    cout << "\n**********************************************************************************************\n";

    
    vector<regevCiphertext> lwe_ct_results = extractRLWECiphertextToLWECiphertext(coeff);
    for (int i = 0; i < ring_dim; i++) {
        int temp = 0;
        for (int j = 0; j < n; j++) {
            long mul_tmp = (lwe_ct_results[i].a[j].ConvertToInt() * lwe_sk[j].ConvertToInt()) % 268369921;
            mul_tmp = mul_tmp < 0 ? mul_tmp + 268369921 : mul_tmp;
            temp = (temp + (int) mul_tmp) % 268369921;
        }
        temp = (temp + lwe_ct_results[i].b.ConvertToInt()) % 268369921;
        cout << temp << " ";

        // if (i % 2 == 0) { // around 16384
        //     int diff = abs(temp - 16384);
        //     stats[diff] += 1;
        // } else { // around 65537
        //     int diff = min(abs(temp - 65537), abs(temp - 0));
        //     stats[diff] += 1;
        // }
    }
    cout << endl;


}