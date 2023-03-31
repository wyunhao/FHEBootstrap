#include "regevEncryption.h"
#include "global.h"
#include "util.h"
#include "seal/seal.h"
#include "seal/util/iterator.h"
#include <numeric>
#include <stdio.h>

using namespace seal;
using namespace std;

// assume lwe_sk_len is a power of 2, and has a square root
Ciphertext evaluatePackedLWECiphertext(const SEALContext& seal_context, vector<regevCiphertext>& lwe_ct_list, const vector<Ciphertext>& lwe_sk_sqrt_list,
                                       const GaloisKeys& gal_keys, const int lwe_sk_len, const int degree = poly_modulus_degree_glb) {
    Evaluator evaluator(seal_context);
    BatchEncoder batch_encoder(seal_context);

    // rotate sqrt(degree), and get sqrt(degree)'s lwe_sk_encrypted
    int sq_rt = sqrt(lwe_sk_len);
    vector<Ciphertext> result(sq_rt);
        
    for (int iter = 0; iter < sq_rt; iter++) {
        for (int j = 0; j < (int) lwe_sk_sqrt_list.size(); j++) {
            vector<uint64_t> lwe_ct_tmp(degree);
            for (int i = 0; i < degree; i++) {
                int ct_index = (i-iter) % (degree/2) < 0 ? (i-iter) % (degree/2) + degree/2 : (i-iter) % (degree/2);
                ct_index = i < degree/2 ? ct_index : ct_index + degree/2;
                int col_index = (i + j*sq_rt) % lwe_sk_len;
                lwe_ct_tmp[i] = lwe_ct_list[ct_index].a[col_index].ConvertToInt();
            }

            Plaintext lwe_ct_pl;
            batch_encoder.encode(lwe_ct_tmp, lwe_ct_pl);
            evaluator.transform_to_ntt_inplace(lwe_ct_pl, lwe_sk_sqrt_list[j].parms_id());

            if (j == 0) {
                evaluator.multiply_plain(lwe_sk_sqrt_list[j], lwe_ct_pl, result[iter]);
            } else {
                Ciphertext temp;
                evaluator.multiply_plain(lwe_sk_sqrt_list[j], lwe_ct_pl, temp);
                evaluator.add_inplace(result[iter], temp);
            }

        }
    }

    for (int i = 0; i < (int) result.size(); i++) {
        evaluator.transform_from_ntt_inplace(result[i]);
    }

    // sum up all sq_rt tmp results to the first one, each first rotate left one and add to the previous tmp result
    for (int iter = sq_rt-1; iter > 0; iter--) {
        evaluator.rotate_rows_inplace(result[iter], 1, gal_keys);
        evaluator.add_inplace(result[iter-1], result[iter]);
    }

    vector<uint64_t> b_parts(degree);
    for(int i = 0; i < degree; i++){
        b_parts[i] = lwe_ct_list[i].b.ConvertToInt();
    }

    Plaintext lwe_b_pl;
    batch_encoder.encode(b_parts, lwe_b_pl);
    evaluator.negate_inplace(result[0]);
    evaluator.add_plain_inplace(result[0], lwe_b_pl);

    return result[0];
}

vector<vector<int>> generateMatrixU_transpose(int n) {
    vector<vector<int>> U(n,  vector<int>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == 0) {
                U[i][j] = (int) power(primitive_root, j, 65537);
            } else if (i == n/2) {
                U[i][j] = (int) modInverse(U[0][j], 65537);
            } else {
                U[i][j] = (int) power(U[i-1][j], 3, 65537);
            }
        }
    }
    return U;
}

void writeUtemp(const vector<uint64_t> U_temp, const int i) {
    ofstream datafile;
    datafile.open ("../data/U_temp/" + to_string(i) + ".txt");

    for(size_t i = 0; i < U_temp.size(); i++){
        datafile << U_temp[i] << "\n";
    }
    datafile.close();
}

vector<uint64_t> readUtemp(const int i, const int degree=poly_modulus_degree_glb) {
    ifstream datafile;
    datafile.open ("../data/U_temp/" + to_string(i) + ".txt");

    vector<uint64_t> U_temp(degree);
    
    for(size_t i = 0; i < U_temp.size(); i++){
        datafile >> U_temp[i];
    }
    datafile.close();

    return U_temp;
}

// assume that degree/2 has a square root
Ciphertext slotToCoeff(const SEALContext& context, vector<Ciphertext>& ct_sqrt_list, vector<Plaintext>& U_plain_list, const GaloisKeys& gal_keys, const int degree=poly_modulus_degree_glb) {
    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);
    // time_start = chrono::high_resolution_clock::now();
    // vector<vector<int>> U = generateMatrixU_transpose(degree);
    // time_end = chrono::high_resolution_clock::now();
    // cout << "Generate U matrix: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
    // cout << "=============================================== \n" << U << "\n=============================================== \n";

    // ofstream datafile;
    // datafile.open ("../data/U_32768.txt");

    // for(size_t i = 0; i < U.size(); i++){
    //     for (size_t j = 0; j < U[0].size(); j++) {
    //         datafile << U[i][j] << "\n";
    //     }
    // }
    // datafile.close();

    int sq_rt = sqrt(degree/2);

    vector<Ciphertext> result(sq_rt);
    for (int iter = 0; iter < sq_rt; iter++) {
        for (int j = 0; j < (int) ct_sqrt_list.size(); j++) {
            // time_start = chrono::high_resolution_clock::now();
            // for (int i = 0; i < degree; i++) {
            //     int row_index = (i-iter) % (degree/2) < 0 ? (i-iter) % (degree/2) + degree/2 : (i-iter) % (degree/2);
            //     row_index = i < degree/2 ? row_index : row_index + degree/2;
            //     int col_index = (i + j*sq_rt) % degree;
            //     U_tmp[i] = U[row_index][col_index];
            // }
            // time_end = chrono::high_resolution_clock::now();
            // total_tmp += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

            // time_start = chrono::high_resolution_clock::now();
            // writeUtemp(U_tmp, j*sq_rt + iter);
            // time_end = chrono::high_resolution_clock::now();
            // total_write += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
            if (j == 0) {
                evaluator.multiply_plain(ct_sqrt_list[j], U_plain_list[iter * ct_sqrt_list.size() + j], result[iter]);
            } else {
                Ciphertext temp;
                evaluator.multiply_plain(ct_sqrt_list[j], U_plain_list[iter * ct_sqrt_list.size() + j], temp);
                evaluator.add_inplace(result[iter], temp);
            }
        }
    }

    for (int i = 0; i < (int) result.size(); i++) {
        evaluator.transform_from_ntt_inplace(result[i]);
    }

    for (int iter = sq_rt-1; iter > 0; iter--) {
        evaluator.rotate_rows_inplace(result[iter], 1, gal_keys);
        evaluator.add_inplace(result[iter-1], result[iter]);
    }

    return result[0];
}


void Bootstrap_RangeCheck_PatersonStockmeyer(Ciphertext& ciphertext, const Ciphertext& input, int modulus, const size_t& degree,
                                             const RelinKeys &relin_keys, const SEALContext& context){
    MemoryPoolHandle my_pool_larger = MemoryPoolHandle::New(true);
    auto old_prof_larger = MemoryManager::SwitchProfile(std::make_unique<MMProfFixed>(std::move(my_pool_larger)));

    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);
    vector<Ciphertext> kCTs(256);
    calUptoDegreeK(kCTs, input, 256, relin_keys, context);

    vector<Ciphertext> kToMCTs(256);
    calUptoDegreeK(kToMCTs, kCTs[kCTs.size()-1], 256, relin_keys, context);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < (int) kCTs.size(); j++) {
            evaluator.mod_switch_to_next_inplace(kCTs[j]);
        }
    }

    // bool first_iter_done = false;
    // for (int i = 255; i >= 0; i--) {
    //     for (int j = 255; j >= 0; j--) {
    //         if (first_iter_done && j == 255) { // not first iteration, but highest degree, use the previous resulted ciphertext as the coeff
    //             evaluator.multiply_inplace(ciphertext, kCTs[j]);
    //         }
    //         else if (rangeCheckIndices_bootstrap[i*256+j] != 0) { // otherwise, always use indices as the coeff
    //             vector<uint64_t> intInd(degree, rangeCheckIndices_bootstrap[i*256+j]);
    //             Plaintext plainInd;
    //             batch_encoder.encode(intInd, plainInd);
    //             if (!first_iter_done) {
    //                 evaluator.multiply_plain(kCTs[j], plainInd, ciphertext);
    //                 first_iter_done = true;
    //             } else {
    //                 Ciphertext tmp;
    //                 evaluator.multiply_plain(kCTs[j], plainInd, tmp);
    //                 evaluator.add_inplace(ciphertext, tmp);
    //             }
    //         }
    //     }
    //     vector<uint64_t> intInd(degree, rangeCheckIndices_bootstrap[(i-1)*256+255]);
    //     Plaintext plainInd;
    //     batch_encoder.encode(intInd, plainInd);
    //     evaluator.add_inplace(ciphertext, tmp);
    // }
    Ciphertext temp_relin(context);
    for(int i = 0; i < 256; i++) {
        Ciphertext levelSum;
        bool flag = false;
        for(int j = 0; j < 256; j++) {
            if(rangeCheckIndices_bootstrap[i*256+j] != 0) {
                vector<uint64_t> intInd(degree, rangeCheckIndices_bootstrap[i*256+j]);
                Plaintext plainInd;
                batch_encoder.encode(intInd, plainInd);
                if (!flag) {
                    evaluator.multiply_plain(kCTs[j], plainInd, levelSum);
                    flag = true;
                } else {
                    Ciphertext tmp;
                    evaluator.multiply_plain(kCTs[j], plainInd, tmp);
                    evaluator.add_inplace(levelSum, tmp);
                }
            }
        }
        evaluator.mod_switch_to_inplace(levelSum, kToMCTs[i].parms_id()); // mod down the plaintext multiplication noise
        if(i == 0) {
            ciphertext = levelSum;
        } else if (i == 1) { // initialize for temp_relin, which is of ct size = 3
            evaluator.multiply(levelSum, kToMCTs[i - 1], temp_relin);
        } else {
            evaluator.multiply_inplace(levelSum, kToMCTs[i - 1]);
            evaluator.add_inplace(temp_relin, levelSum);
        }
    }

    evaluator.relinearize_inplace(temp_relin, relin_keys);
    evaluator.add_inplace(ciphertext, temp_relin);

    vector<uint64_t> intInd(degree, 0); 
    Plaintext plainInd;
    Ciphertext temp;
    batch_encoder.encode(intInd, plainInd);
    evaluator.negate_inplace(ciphertext);
    evaluator.add_plain_inplace(ciphertext, plainInd);
    temp.release();
    temp_relin.release();
    for(int i = 0; i < 256; i++){
        kCTs[i].release();
        kToMCTs[i].release();
    }
    MemoryManager::SwitchProfile(std::move(old_prof_larger));
}


Ciphertext encryptLWEskUnderBFV(const SEALContext& context, const size_t& degree,
                                const PublicKey& BFVpk, const SecretKey& BFVsk,
                                const regevSK& regSk, const regevParam& params) { 
    Ciphertext switchingKey(context);

    BatchEncoder batch_encoder(context);
    Encryptor encryptor(context, BFVpk);
    encryptor.set_secret_key(BFVsk);

    int tempn = 1;
    for(tempn = 1; tempn < params.n; tempn *= 2){}

    vector<uint64_t> skInt(degree);
    for(size_t i = 0; i < degree; i++){
        auto tempindex = i%uint64_t(tempn);
        if(int(tempindex) >= params.n) {
            skInt[i] = 0;
        } else {
            skInt[i] = uint64_t(regSk[tempindex].ConvertToInt() % params.q);
        }
    }
    Plaintext plaintext;
    batch_encoder.encode(skInt, plaintext);
    encryptor.encrypt_symmetric(plaintext, switchingKey);

    return switchingKey;
}


vector<regevCiphertext> extractRLWECiphertextToLWECiphertext(Ciphertext& rlwe_ct, const int ring_dim = poly_modulus_degree_glb, const int n = 1024, const int p = 65537, const int big_prime = 268369921) {
    vector<regevCiphertext> results(ring_dim);

    prng_seed_type seed;
    for (auto &i : seed) {
        i = random_uint64();
    }
    auto rng = std::make_shared<Blake2xbPRNGFactory>(Blake2xbPRNGFactory(seed));
    RandomToStandardAdapter engine(rng->create());
    uniform_int_distribution<uint32_t> dist(0, 100);

    for (int cnt = 0; cnt < ring_dim; cnt++) {
        results[cnt].a = NativeVector(n);
        int ind = 0;
        for (int i = cnt; i >= 0 && ind < n; i--) {
            float temp_f = ((float) rlwe_ct.data(1)[i]) * ((float) p) / ((float) big_prime);
            uint32_t decimal = (temp_f - ((int) temp_f)) * 100;
            float rounding = dist(engine) < decimal ? 1 : 0;

            long temp = ((int) (temp_f + rounding)) % p;
            results[cnt].a[ind] = temp < 0 ? p + temp : temp;

            ind++;
        }

        for (int i = ring_dim-1; i > ring_dim - n + cnt && ind < n; i--) {
            float temp_f = ((float) rlwe_ct.data(1)[i]) * ((float) p) / ((float) big_prime);
            uint32_t decimal = (temp_f - ((int) temp_f)) * 100;
            float rounding = dist(engine) < decimal ? 1 : 0;

            long temp = ((int) (temp_f + rounding)) % p;
            results[cnt].a[ind] = -temp < 0 ? p-temp : -temp;

            ind++;
        }

        float temp_f = ((float) rlwe_ct.data(0)[cnt]) * ((float) p) / ((float) big_prime);
        uint32_t decimal = temp_f - ((int) temp_f) * 100;
        float rounding = dist(engine) < decimal ? 1 : 0;

        long temp = ((int) (temp_f + rounding)) % p;
        results[cnt].b = temp % ((int) p);
    }

    return results;
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


    // Ciphertext coefft(seal_context);
    // vector<uint64_t> vc = {0, 0, 0, 1, 0, 0, 0, 0,
    //                        1, 0, 0, 0, 0, 0, 0, 0};

    // Plaintext tt(ring_dim, 0);
    // batch_encoder.encode(vc, tt);
    // encryptor.encrypt(tt, coefft);


    // vector<parms_id_type> id_list(12);
    // int cnt = 0;


    // auto context_data = seal_context.first_context_data();
    // while (context_data)
    // {
    //     cout << " Level (chain index): " << context_data->chain_index();
    //     if (context_data->parms_id() == seal_context.first_parms_id())
    //     {
    //         cout << " ...... first_context_data()" << endl;
    //         // id_list[cnt] = seal_context.first_parms_id();
    //         // cnt++;
    //     }
    //     else if (context_data->parms_id() == seal_context.last_parms_id())
    //     {
    //         cout << " ...... last_context_data()" << endl;
    //         // id_list[id_list.size()-1] = seal_context.last_parms_id();
    //     }
    //     else
    //     {
    //         cout << endl;
    //     }
    //     // id_list[cnt] = context_data->parms_id();
    //     // cnt++;
    //     cout << "      parms_id: " << context_data->parms_id() << endl;
    //     cout << "      coeff_modulus primes: ";
    //     cout << hex;
    //     for (const auto &prime : context_data->parms().coeff_modulus())
    //     {
    //         cout << prime.value() << " ";
    //     }
    //     cout << dec << endl;
    //     cout << "\\" << endl;
    //     cout << " \\-->";

    //     /*
    //     Step forward in the chain.
    //     */
    //     context_data = context_data->next_context_data();
    // }
    // cout << " End of chain reached" << endl << endl;

    // int deg = 256;
    // vector<Ciphertext> results(deg);
    // cout << coefft.parms_id() << endl;
    // calUptoDegreeK(results, coefft, deg, relin_keys, seal_context);

    // cout << results[0].parms_id() << endl;

    // for (int i = 0; i < deg; i++) {
    //     auto it = find(id_list.begin(), id_list.end(), results[i].parms_id());
    //     if (it != id_list.end()) {
    //         cout << i << ", " << it - id_list.begin()  <<" " << results[i].parms_id() << endl;
    //     } else {
    //         cout << "Not in list???????\n";
    //     }
    // }

    // for (int i = 0; i < id_list.size(); i++) {
    //     cout << id_list[i] << endl;
    // }
    // cout << endl << endl;


    // auto it_coefft = find(id_list.begin(), id_list.end(), coefft.parms_id());
    // if (it_coefft != id_list.end()) {
    //     cout << "coefft: " << it_coefft - id_list.begin() << endl;
    // } else {
    //     cout << "Not in list???????\n";
    // }


    // cout << endl << endl;


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










    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int sq_sk = sqrt(n), sq_ct = sqrt(ring_dim/2);
    vector<Ciphertext> lwe_sk_sqrt_list(sq_sk), ct_sqrt_list(2*sq_ct);

    chrono::high_resolution_clock::time_point time_start, time_end;
    int total_preprocess = 0, total_online = 0;

    Ciphertext lwe_sk_column;

    time_start = chrono::high_resolution_clock::now();
    evaluator.rotate_columns(lwe_sk_encrypted, gal_keys, lwe_sk_column);
    for (int i = 0; i < sq_sk; i++) {
        evaluator.rotate_rows(lwe_sk_encrypted, sq_sk * i, gal_keys, lwe_sk_sqrt_list[i]);
        evaluator.transform_to_ntt_inplace(lwe_sk_sqrt_list[i]);
    }
    time_end = chrono::high_resolution_clock::now();
    total_preprocess += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();



    //////////////////////////////////////////// EVALUATION UNDER BFV ////////////////////////////////////////////////
    time_start = chrono::high_resolution_clock::now();
    Ciphertext result = evaluatePackedLWECiphertext(seal_context, lwe_ct_list, lwe_sk_sqrt_list, gal_keys, n, ring_dim);
    time_end = chrono::high_resolution_clock::now();
    total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    cout << "TOTAL TIME for evaluation: " << total_online << endl;
    // cout << "Noise: " << decryptor.invariant_noise_budget(result) << " bits\n";

    Ciphertext range_check_res;
    time_start = chrono::high_resolution_clock::now();
    Bootstrap_RangeCheck_PatersonStockmeyer(range_check_res, result, p, ring_dim, relin_keys, seal_context);
    time_end = chrono::high_resolution_clock::now();
    total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    cout << "TOTAL TIME for rangecheck: " << total_online << endl;

    // decryptor.decrypt(range_check_res, t);
    // batch_encoder.decode(t, v);
    // cout << ">>>>>>>>>>>> after range check" << v << endl << endl;


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
    // cout << "Noise after slotToCoeff: " << decryptor.invariant_noise_budget(coeff) << " bits\n";

    // cout << "Coeff decoded: ";
    // decryptor.decrypt(coeff, t);
    // for (int i = 0; i < ring_dim; i++) {
    //     cout << t[i] << " ";
    // }
    // cout << endl;

    // batch_encoder.decode(t, v);

    // for (int i = 0; i < ring_dim; i++) {
    //     cout << v[i] << " ";
    // }
    // cout << endl;

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


    vector<uint64_t> stats(100, 0);
    for (int i = 0; i < ring_dim; i++) {
        int temp = 0;
        for (int j = 0; j < n; j++) {
            long mul_tmp = (lwe_ct_results[i].a[j].ConvertToInt() * lwe_sk[j].ConvertToInt()) % p;
            mul_tmp = mul_tmp < 0 ? mul_tmp + p : mul_tmp;
            temp = (temp + (int) mul_tmp) % p;
        }
        temp = (temp + lwe_ct_results[i].b.ConvertToInt()) % p;
        if (i % 2 == 0) { // around 16384
            int diff = abs(temp - 16384);
            stats[diff] += 1;
        } else { // around 65537
            int diff = min(abs(temp - 65537), abs(temp - 0));
            stats[diff] += 1;
        }
    }
    // }

    cout << endl << "stats: " << stats << endl;
}