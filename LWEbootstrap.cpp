#include "regevEncryption.h"
#include "global.h"
#include "util.h"
#include "seal/seal.h"
#include "seal/util/iterator.h"
#include <numeric>
#include <stdio.h>

using namespace seal;
using namespace std;


vector<regevCiphertext> bootstrap(vector<regevCiphertext>& lwe_ct_list, Ciphertext& lwe_sk_encrypted, const SEALContext& seal_context,
                                  const RelinKeys& relin_keys, const GaloisKeys& gal_keys, const int ring_dim, const int n,
                                  const int p, const KSwitchKeys& ksk, const MemoryPoolHandle& my_pool) {
    chrono::high_resolution_clock::time_point time_start, time_end;
    int total_preprocess = 0, total_online = 0;

    Evaluator evaluator(seal_context);
    BatchEncoder batch_encoder(seal_context);

    int sq_sk = sqrt(n), sq_ct = sqrt(ring_dim/2);
    vector<Ciphertext> lwe_sk_sqrt_list(sq_sk), ct_sqrt_list(2*sq_ct);

    Ciphertext lwe_sk_column;

    time_start = chrono::high_resolution_clock::now();
    evaluator.rotate_columns(lwe_sk_encrypted, gal_keys, lwe_sk_column);
    for (int i = 0; i < sq_sk; i++) {
        evaluator.rotate_rows(lwe_sk_encrypted, sq_sk * i, gal_keys, lwe_sk_sqrt_list[i]);
        evaluator.transform_to_ntt_inplace(lwe_sk_sqrt_list[i]);
    }
    time_end = chrono::high_resolution_clock::now();
    total_preprocess += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();


    time_start = chrono::high_resolution_clock::now();
    Ciphertext result = evaluatePackedLWECiphertext(seal_context, lwe_ct_list, lwe_sk_sqrt_list, gal_keys, n, ring_dim);
    time_end = chrono::high_resolution_clock::now();
    total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    cout << "TOTAL TIME for evaluation: " << total_online << endl;
    // cout << "Noise: " << decryptor.invariant_noise_budget(result) << " bits\n";

    Ciphertext range_check_res;
    time_start = chrono::high_resolution_clock::now();
    Bootstrap_RangeCheck_PatersonStockmeyer(range_check_res, result, rangeCheckIndices_bootstrap, p, ring_dim, relin_keys, seal_context);
    time_end = chrono::high_resolution_clock::now();
    total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    cout << "TOTAL TIME for rangecheck: " << total_online << endl;

    ////////////////////////////////////////// SLOT TO COEFFICIENT /////////////////////////////////////////////////////

    time_start = chrono::high_resolution_clock::now();
    evaluator.mod_switch_to_next_inplace(range_check_res);
    time_end = chrono::high_resolution_clock::now();
    total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    // cout << "Noise after range check: " << decryptor.invariant_noise_budget(range_check_res) << " bits\n";

    time_start = chrono::high_resolution_clock::now();
    Ciphertext range_check_res_copy(range_check_res);

    evaluator.rotate_columns_inplace(range_check_res_copy, gal_keys);
    for (int i = 0; i < sq_ct; i++) {
        evaluator.rotate_rows(range_check_res, sq_ct * i, gal_keys, ct_sqrt_list[i]);
        evaluator.transform_to_ntt_inplace(ct_sqrt_list[i]);
        evaluator.rotate_rows(range_check_res_copy, sq_ct * i, gal_keys, ct_sqrt_list[i+sq_ct]);
        evaluator.transform_to_ntt_inplace(ct_sqrt_list[i+sq_ct]);
    }

    vector<Plaintext> U_plain_list(ring_dim);
    for (int iter = 0; iter < sq_ct; iter++) {
        for (int j = 0; j < (int) ct_sqrt_list.size(); j++) {
            vector<uint64_t> U_tmp = readUtemp(j*sq_ct + iter);
            batch_encoder.encode(U_tmp, U_plain_list[iter * ct_sqrt_list.size() + j]);
            evaluator.transform_to_ntt_inplace(U_plain_list[iter * ct_sqrt_list.size() + j], ct_sqrt_list[j].parms_id());
        }
    }
    time_end = chrono::high_resolution_clock::now();
    total_preprocess += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

    time_start = chrono::high_resolution_clock::now();
    Ciphertext coeff = slotToCoeff(seal_context, ct_sqrt_list, U_plain_list, gal_keys, ring_dim);
    time_end = chrono::high_resolution_clock::now();
    total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    cout << "TOTAL TIME for slotToCoeff: " << total_online << endl;

    ////////////////////////////////////////////////// KEY SWITCHING ///////////////////////////////////////////////////

    time_start = chrono::high_resolution_clock::now();
    while(seal_context.last_parms_id() != coeff.parms_id()){
        evaluator.mod_switch_to_next_inplace(coeff);
    }
    // cout << "Noise before key switch: " << decryptor.invariant_noise_budget(coeff) << " bits\n";

    Ciphertext copy_coeff = coeff;
    auto ct_in_iter = util::iter(copy_coeff);
    ct_in_iter += coeff.size() - 1;
    seal::util::set_zero_poly(ring_dim, 1, coeff.data(1)); // notice that the coeff_mod.size() is hardcoded to 1, thus this needs to be performed on the last level

    evaluator.switch_key_inplace(coeff, *ct_in_iter, static_cast<const KSwitchKeys &>(ksk), 0, my_pool);

    // cout << "Noise before extraction: " << decryptor.invariant_noise_budget(coeff) << " bits\n";
    
    vector<regevCiphertext> lwe_ct_results = extractRLWECiphertextToLWECiphertext(coeff);
    time_end = chrono::high_resolution_clock::now();
    total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    cout << "TOTAL TIME for Extraction: " << total_online << endl;
    cout << "TOTAL PREPROCESS TIME: " << total_preprocess << endl;
    cout << "TOTAL ONLINE TIME: " << total_online << endl;

    return lwe_ct_results;
}



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
    cout << "***************************************** Ternary SK *****************************************\n";
    for (int i = 0; i < n; i++) {
        lwe_sk[i] = new_key.data()[i] > 65537 ? 65536 : new_key.data()[i];
        cout << lwe_sk[i] << " ";
    }
    cout << "\n**********************************************************************************************\n";

    seal::util::RNSIter new_key_rns(new_key.data().data(), ring_dim);
    ntt_negacyclic_harvey(new_key_rns, coeff_modulus.size(), seal_context.key_context_data()->small_ntt_tables());

    vector<regevCiphertext> lwe_ct_list = regevGeneratePublicKey(lwe_params, lwe_sk, true);


    ////////////////////////////////////////////// ENCRYPT SK UNDER BFV ////////////////////////////////////////////////

    // one switching key for one lwe_sk
    Ciphertext lwe_sk_encrypted = encryptLWEskUnderBFV(seal_context, ring_dim, bfv_public_key, bfv_secret_key, lwe_sk, lwe_params);
    // cout << "Noise SK: " << decryptor.invariant_noise_budget(lwe_sk_encrypted) << " bits\n";



















    vector<regevCiphertext> lwe_ct_results = bootstrap(lwe_ct_list, lwe_sk_encrypted, seal_context, relin_keys, gal_keys,
                                                       ring_dim, n, p, ksk, my_pool);





    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // int sq_sk = sqrt(n), sq_ct = sqrt(ring_dim/2);
    // vector<Ciphertext> lwe_sk_sqrt_list(sq_sk), ct_sqrt_list(2*sq_ct);

    // chrono::high_resolution_clock::time_point time_start, time_end;
    // int total_preprocess = 0, total_online = 0;

    // Ciphertext lwe_sk_column;

    // time_start = chrono::high_resolution_clock::now();
    // evaluator.rotate_columns(lwe_sk_encrypted, gal_keys, lwe_sk_column);
    // for (int i = 0; i < sq_sk; i++) {
    //     evaluator.rotate_rows(lwe_sk_encrypted, sq_sk * i, gal_keys, lwe_sk_sqrt_list[i]);
    //     evaluator.transform_to_ntt_inplace(lwe_sk_sqrt_list[i]);
    // }
    // time_end = chrono::high_resolution_clock::now();
    // total_preprocess += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();



    // //////////////////////////////////////////// EVALUATION UNDER BFV ////////////////////////////////////////////////
    // time_start = chrono::high_resolution_clock::now();
    // Ciphertext result = evaluatePackedLWECiphertext(seal_context, lwe_ct_list, lwe_sk_sqrt_list, gal_keys, n, ring_dim);
    // time_end = chrono::high_resolution_clock::now();
    // total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    // cout << "TOTAL TIME for evaluation: " << total_online << endl;
    // // cout << "Noise: " << decryptor.invariant_noise_budget(result) << " bits\n";

    // Ciphertext range_check_res;
    // time_start = chrono::high_resolution_clock::now();
    // Bootstrap_RangeCheck_PatersonStockmeyer(range_check_res, result, rangeCheckIndices_bootstrap, p, ring_dim, relin_keys, seal_context);
    // time_end = chrono::high_resolution_clock::now();
    // total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    // cout << "TOTAL TIME for rangecheck: " << total_online << endl;

    // ////////////////////////////////////////// SLOT TO COEFFICIENT /////////////////////////////////////////////////////

    // time_start = chrono::high_resolution_clock::now();
    // evaluator.mod_switch_to_next_inplace(range_check_res);
    // time_end = chrono::high_resolution_clock::now();
    // total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    // // cout << "Noise after range check: " << decryptor.invariant_noise_budget(range_check_res) << " bits\n";

    // time_start = chrono::high_resolution_clock::now();
    // Ciphertext range_check_res_copy(range_check_res);

    // evaluator.rotate_columns_inplace(range_check_res_copy, gal_keys);
    // for (int i = 0; i < sq_ct; i++) {
    //     evaluator.rotate_rows(range_check_res, sq_ct * i, gal_keys, ct_sqrt_list[i]);
    //     evaluator.transform_to_ntt_inplace(ct_sqrt_list[i]);
    //     evaluator.rotate_rows(range_check_res_copy, sq_ct * i, gal_keys, ct_sqrt_list[i+sq_ct]);
    //     evaluator.transform_to_ntt_inplace(ct_sqrt_list[i+sq_ct]);
    // }

    // vector<Plaintext> U_plain_list(ring_dim);
    // for (int iter = 0; iter < sq_ct; iter++) {
    //     for (int j = 0; j < (int) ct_sqrt_list.size(); j++) {
    //         vector<uint64_t> U_tmp = readUtemp(j*sq_ct + iter);
    //         batch_encoder.encode(U_tmp, U_plain_list[iter * ct_sqrt_list.size() + j]);
    //         evaluator.transform_to_ntt_inplace(U_plain_list[iter * ct_sqrt_list.size() + j], ct_sqrt_list[j].parms_id());
    //     }
    // }
    // time_end = chrono::high_resolution_clock::now();
    // total_preprocess += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

    // time_start = chrono::high_resolution_clock::now();
    // Ciphertext coeff = slotToCoeff(seal_context, ct_sqrt_list, U_plain_list, gal_keys, ring_dim);
    // time_end = chrono::high_resolution_clock::now();
    // total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    // cout << "TOTAL TIME for slotToCoeff: " << total_online << endl;

    // Plaintext t;
    // cout << "Coeff decoded: ";
    // decryptor.decrypt(coeff, t);
    // for (int i = 0; i < ring_dim; i++) {
    //     cout << t[i] << " ";
    // }
    // cout << endl;

    // ////////////////////////////////////////////////// KEY SWITCHING ///////////////////////////////////////////////////

    // time_start = chrono::high_resolution_clock::now();
    // while(seal_context.last_parms_id() != coeff.parms_id()){
    //     evaluator.mod_switch_to_next_inplace(coeff);
    // }
    // // cout << "Noise before key switch: " << decryptor.invariant_noise_budget(coeff) << " bits\n";

    // Ciphertext copy_coeff = coeff;
    // auto ct_in_iter = util::iter(copy_coeff);
    // ct_in_iter += coeff.size() - 1;
    // seal::util::set_zero_poly(ring_dim, 1, coeff.data(1)); // notice that the coeff_mod.size() is hardcoded to 1, thus this needs to be performed on the last level

    // evaluator.switch_key_inplace(coeff, *ct_in_iter, static_cast<const KSwitchKeys &>(ksk), 0, my_pool);

    // // cout << "Noise before extraction: " << decryptor.invariant_noise_budget(coeff) << " bits\n";
    
    // vector<regevCiphertext> lwe_ct_results = extractRLWECiphertextToLWECiphertext(coeff);
    // time_end = chrono::high_resolution_clock::now();
    // total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    // cout << "TOTAL TIME for Extraction: " << total_online << endl;
    // cout << "TOTAL PREPROCESS TIME: " << total_preprocess << endl;
    // cout << "TOTAL ONLINE TIME: " << total_online << endl;





















    vector<uint64_t> stats(100, 0);
    for (int i = 0; i < ring_dim; i++) {
        int temp = 0;
        for (int j = 0; j < n; j++) {
            long mul_tmp = (lwe_ct_results[i].a[j].ConvertToInt() * lwe_sk[j].ConvertToInt()) % p;
            mul_tmp = mul_tmp < 0 ? mul_tmp + p : mul_tmp;
            temp = (temp + (int) mul_tmp) % p;
        }
        temp = (temp + lwe_ct_results[i].b.ConvertToInt()) % p;
        cout << temp << " ";
        if (i % 2 == 0) { // around 16384
            int diff = abs(temp - 16384);
            stats[diff] += 1;
        } else { // around 65537
            int diff = min(abs(temp - 65537), abs(temp - 0));
            stats[diff] += 1;
        }
    }

    cout << endl << "stats: " << stats << endl;
}