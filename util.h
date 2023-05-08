
using namespace std;
using namespace seal;


struct BootstrapParam {
    int ciphertextSpacePrime;
    int errorRange;
    int plaintextSpace;

    int firstLevelDegree; // used for CallToDegreeK in range check procedure
    int secondLevelDegree;

    BootstrapParam(int ciphertextSpacePrime, int errorRange, int plaintextSpace, int firstLevelDegree, int secondLevelDegree)
    : ciphertextSpacePrime(ciphertextSpacePrime), errorRange(errorRange), plaintextSpace(plaintextSpace),
      firstLevelDegree(firstLevelDegree), secondLevelDegree(secondLevelDegree)
    {}
};


inline
long power(long x, long y, long m)
{
    if (y == 0)
        return 1;
    long p = power(x, y / 2, m) % m;
    p = (p * p) % m;
 
    return (y % 2 == 0) ? p : (x * p) % m;
}

inline
long modInverse(long a, long m)
{
    return power(a, m - 2, m);
}

inline void calUptoDegreeK(vector<Ciphertext>& output, const Ciphertext& input, const int DegreeK, const RelinKeys &relin_keys,
                           const SEALContext& context, const bool skip_odd=false) {
    vector<int> calculated(DegreeK, 0);
    Evaluator evaluator(context);
    output[0] = input;
    calculated[0] = 1; // degree 1, x
    Ciphertext res, base;
    vector<int> numMod(DegreeK, 0);

    for(int i = DegreeK; i > 0; i--){
        // if (i == 512 && skip) {
        //     cout << "Calculated: " << calculated << endl;
        //     break;
        // }
        if (skip_odd && i % 2 == 1) { // 0 is for degree 1, 1 is for degree 2, skip all 2k+1 degree
            calculated[i-1] = 1;
            output[i-1] = input;
        } else if(calculated[i-1] == 0){
            auto toCalculate = i;
            int resdeg = 0;
            int basedeg = 1;
            base = input;
            while(toCalculate > 0){
                if(toCalculate & 1){
                    toCalculate -= 1;
                    resdeg += basedeg;
                    if(calculated[resdeg-1] != 0){
                        res = output[resdeg - 1];
                    } else {
                        if(resdeg == basedeg){
                            res = base; // should've never be used as base should have made calculated[basedeg-1]
                        } else {
                            numMod[resdeg-1] = numMod[basedeg-1];

                            evaluator.mod_switch_to_inplace(res, base.parms_id()); // match modulus
                            evaluator.multiply_inplace(res, base);
                            evaluator.relinearize_inplace(res, relin_keys);
                            while(numMod[resdeg-1] < (ceil(log2(resdeg))/2)){
                                evaluator.mod_switch_to_next_inplace(res);
                                numMod[resdeg-1]+=1;
                            }
                        }
                        output[resdeg-1] = res;
                        calculated[resdeg-1] += 1;
                    }
                } else {
                    toCalculate /= 2;
                    basedeg *= 2;
                    if(calculated[basedeg-1] != 0){
                        base = output[basedeg - 1];
                    } else {
                        numMod[basedeg-1] = numMod[basedeg/2-1];
                        evaluator.square_inplace(base);
                        evaluator.relinearize_inplace(base, relin_keys);
                        while(numMod[basedeg-1] < (ceil(log2(basedeg))/2)){
                            evaluator.mod_switch_to_next_inplace(base);
                            numMod[basedeg-1]+=1;
                        }
                        output[basedeg-1] = base;
                        calculated[basedeg-1] += 1;
                    }
                }
            }
        }
    }

    // if (!skip) {
    for(size_t i = 0; i < output.size()-1; i++){
        evaluator.mod_switch_to_inplace(output[i], output[output.size()-1].parms_id()); // match modulus
    }
    // }
    return;
}


inline void calUptoDegreeK_bigPrime(vector<Ciphertext>& output, const Ciphertext& input, const int DegreeK, const RelinKeys &relin_keys,
                                    const SEALContext& context, map<int, bool>& modDownIndices, const bool skip_odd=false) {
    vector<int> calculated(DegreeK, 0);
    Evaluator evaluator(context);
    output[0] = input;
    calculated[0] = 1; // degree 1, x
    Ciphertext res, base;
    vector<int> numMod(DegreeK, 0);

    for(int i = DegreeK; i > 0; i--){
        if (skip_odd && i % 2 == 1) { // 0 is for degree 1, 1 is for degree 2, skip all 2k+1 degree
            calculated[i-1] = 1;
            output[i-1] = input;
        } else if(calculated[i-1] == 0){
            auto toCalculate = i;
            int resdeg = 0;
            int basedeg = 1;
            base = input;
            while(toCalculate > 0){
                if(toCalculate & 1){
                    toCalculate -= 1;
                    resdeg += basedeg;
                    if(calculated[resdeg-1] != 0){
                        res = output[resdeg - 1];
                    } else {
                        if(resdeg == basedeg){
                            res = base; // should've never be used as base should have made calculated[basedeg-1]
                        } else {
                            numMod[resdeg-1] = numMod[basedeg-1];

                            evaluator.mod_switch_to_inplace(res, base.parms_id()); // match modulus
                            evaluator.multiply_inplace(res, base);
                            evaluator.relinearize_inplace(res, relin_keys);
                            if(modDownIndices.count(resdeg) && !modDownIndices[resdeg]) {
                                modDownIndices[resdeg] = true;
                                cout << "mod down on " << resdeg << endl;
                                evaluator.mod_switch_to_next_inplace(res);
                                numMod[resdeg-1]+=1;
                            }
                        }
                        output[resdeg-1] = res;
                        calculated[resdeg-1] += 1;
                    }
                } else {
                    toCalculate /= 2;
                    basedeg *= 2;
                    if(calculated[basedeg-1] != 0){
                        base = output[basedeg - 1];
                    } else {
                        numMod[basedeg-1] = numMod[basedeg/2-1];
                        evaluator.square_inplace(base);
                        evaluator.relinearize_inplace(base, relin_keys);
                        while(modDownIndices.count(basedeg) && !modDownIndices[basedeg]) {
                            modDownIndices[basedeg] = true;
                            cout << "mod down on " << basedeg << endl;
                            evaluator.mod_switch_to_next_inplace(base);
                            numMod[basedeg-1]+=1;
                        }
                        output[basedeg-1] = base;
                        calculated[basedeg-1] += 1;
                    }
                }
            }
        }
    }

    cout << "Mod level: " << numMod << endl;

    // if (!skip) {
    for(size_t i = 0; i < output.size()-1; i++){
        evaluator.mod_switch_to_inplace(output[i], output[output.size()-1].parms_id()); // match modulus
    }
    // }
    return;
}


// assume lwe_sk_len is a power of 2, and has a square root
Ciphertext evaluatePackedLWECiphertext(const SEALContext& seal_context, vector<regevCiphertext>& lwe_ct_list,
                                       const vector<Ciphertext>& lwe_sk_sqrt_list, const GaloisKeys& gal_keys, const int lwe_sk_len,
                                       const vector<uint64_t>& q_shift_constant, const int degree = poly_modulus_degree_glb,
                                       const bool gateEval = false, const int q = 65537) {
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

    if (gateEval) {
        Plaintext q_shift_pl;
        batch_encoder.encode(q_shift_constant, q_shift_pl);
        evaluator.add_plain_inplace(result[0], q_shift_pl);
    }

    return result[0];
}

vector<vector<int>> generateMatrixU_transpose(int n, const int q = 65537) {
    vector<vector<int>> U(n,  vector<int>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == 0) {
                U[i][j] = (int) power(primitive_root, j, q);
            } else if (i == n/2) {
                U[i][j] = (int) modInverse(U[0][j], q);
            } else {
                U[i][j] = (int) power(U[i-1][j], 3, q);
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
/**
 * @brief Pre-processed Version.
 * 
 * @param context 
 * @param ct_sqrt_list 
 * @param U_plain_list_c 
 * @param gal_keys 
 * @param degree 
 * @return Ciphertext 
 */
Ciphertext slotToCoeff(const SEALContext& context, vector<Ciphertext>& ct_sqrt_list, vector<Plaintext>& U_plain_list_c, const GaloisKeys& gal_keys, const int degree=poly_modulus_degree_glb) {
    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);
    int sq_rt = sqrt(degree/2);
    vector<Plaintext> U_plain_list(U_plain_list_c);

    chrono::high_resolution_clock::time_point time_start, time_end;
    int total_mul = 0;

    vector<Ciphertext> result(sq_rt);
    for (int iter = 0; iter < sq_rt; iter++) {
        for (int j = 0; j < (int) ct_sqrt_list.size(); j++) {
            time_start = chrono::high_resolution_clock::now();
            if (j == 0) {
                evaluator.multiply_plain(ct_sqrt_list[j], U_plain_list[iter * ct_sqrt_list.size() + j], result[iter]);
            } else {
                Ciphertext temp;
                evaluator.multiply_plain(ct_sqrt_list[j], U_plain_list[iter * ct_sqrt_list.size() + j], temp);
                evaluator.add_inplace(result[iter], temp);
            }
            time_end = chrono::high_resolution_clock::now();
            total_mul += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
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


/**
 * @brief Non Pre-processed Version.
 * 
 * @param context 
 * @param ct_sqrt_list 
 * @param U_plain_list_c 
 * @param gal_keys 
 * @param degree 
 * @return Ciphertext 
 */
Ciphertext slotToCoeff_WOPrepreocess(const SEALContext& context, vector<Ciphertext>& ct_sqrt_list, const GaloisKeys& gal_keys, const int degree=poly_modulus_degree_glb) {
    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);
    vector<vector<int>> U = generateMatrixU_transpose(degree);
    int sq_rt = sqrt(degree/2);

    chrono::high_resolution_clock::time_point time_start, time_end;
    int total_mul = 0;

    vector<Ciphertext> result(sq_rt);
    for (int iter = 0; iter < sq_rt; iter++) {
        for (int j = 0; j < (int) ct_sqrt_list.size(); j++) {
            vector<uint64_t> U_tmp(degree);
            for (int i = 0; i < degree; i++) {
                int row_index = (i-iter) % (degree/2) < 0 ? (i-iter) % (degree/2) + degree/2 : (i-iter) % (degree/2);
                row_index = i < degree/2 ? row_index : row_index + degree/2;
                int col_index = (i + j*sq_rt) % (degree/2);
                if (j < (int) ct_sqrt_list.size() / 2) { // first half
                    col_index = i < degree/2 ? col_index : col_index + degree/2;
                } else {
                    col_index = i < degree/2 ? col_index + degree/2 : col_index;
                }
                U_tmp[i] = U[row_index][col_index];
            }
            writeUtemp(U_tmp, j*sq_rt + iter);

            Plaintext U_plain;
            batch_encoder.encode(U_tmp, U_plain);
            evaluator.transform_to_ntt_inplace(U_plain, ct_sqrt_list[j].parms_id());

            if (j == 0) {
                evaluator.multiply_plain(ct_sqrt_list[j], U_plain, result[iter]);
            } else {
                Ciphertext temp;
                evaluator.multiply_plain(ct_sqrt_list[j], U_plain, temp);
                evaluator.add_inplace(result[iter], temp);
            }
            time_end = chrono::high_resolution_clock::now();
            total_mul += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
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




void Bootstrap_RangeCheck_PatersonStockmeyer_bigPrime(Ciphertext& ciphertext, const Ciphertext& input, const vector<uint64_t>& rangeCheckIndices,
                                             const int modulus, const size_t& degree, const RelinKeys &relin_keys, const SEALContext& context, const SecretKey& bfv_secret_key,
                                             const int f_zero = 0, const bool gateEval = false, const bool skip_first_odd = true,
                                             const int firstDegree = 256, const int secondDegree = 256) {
    MemoryPoolHandle my_pool_larger = MemoryPoolHandle::New(true);
    auto old_prof_larger = MemoryManager::SwitchProfile(std::make_unique<MMProfFixed>(std::move(my_pool_larger)));

    Decryptor decryptor(context, bfv_secret_key);
    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);
    vector<Ciphertext> kCTs(firstDegree);

    map<int, bool> modDownIndices_firstLevel = {{4, false}, {16, false}, {64, false}, {256, false}};
    map<int, bool> modDownIndices_secondLevel = {{2, false}, {8, false}, {32, false}, {64, false}, {128, false}, {512, false}};

    calUptoDegreeK_bigPrime(kCTs, input, firstDegree, relin_keys, context, modDownIndices_firstLevel, skip_first_odd);
    cout << "KCT :\n";
    Plaintext pl;
    vector<uint64_t> v(degree);
    for (int i = 0; i < kCTs.size(); i++) {
        decryptor.decrypt(kCTs[i], pl);
        batch_encoder.decode(pl, v);
        cout << v[1] << ",";
    }
    cout << endl;


    vector<Ciphertext> kToMCTs(secondDegree);
    calUptoDegreeK_bigPrime(kToMCTs, kCTs[kCTs.size()-1], secondDegree, relin_keys, context, modDownIndices_secondLevel);

    cout << "KToMCT :\n";
    for (int i = 0; i < kToMCTs.size(); i++) {
        decryptor.decrypt(kToMCTs[i], pl);
        batch_encoder.decode(pl, v);
        cout << v[1] << ",";
    }
    cout << endl;
    cout << "Noise of last: " << decryptor.invariant_noise_budget(kToMCTs[kToMCTs.size()-1]) << " bits\n";

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < (int) kCTs.size(); j++) {
            evaluator.mod_switch_to_next_inplace(kCTs[j]);
        }
    }

    Ciphertext temp_relin(context);

    for(int i = 0; i < secondDegree; i++) {
        Ciphertext levelSum;
        bool flag = false;
        for(int j = 0; j < firstDegree; j++) {
            if(rangeCheckIndices[i*firstDegree+j] != 0) {
                vector<uint64_t> intInd(degree, rangeCheckIndices[i*firstDegree+j]);
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

    for(int i = 0; i < firstDegree; i++){
        kCTs[i].release();
    }
    for(int i = 0; i < secondDegree; i++){
        kToMCTs[i].release();
    }

    evaluator.relinearize_inplace(temp_relin, relin_keys);
    evaluator.add_inplace(ciphertext, temp_relin);
    temp_relin.release();

    vector<uint64_t> intInd(degree, f_zero);
    Plaintext plainInd;
    batch_encoder.encode(intInd, plainInd);
    evaluator.negate_inplace(ciphertext);
    evaluator.add_plain_inplace(ciphertext, plainInd);

    decryptor.decrypt(ciphertext, pl);
    batch_encoder.decode(pl, v);
    cout << "Decrypt of function should be all ZERO: \n" << v << endl;
    cout << "Noise after function: " << decryptor.invariant_noise_budget(ciphertext) << " bits\n";

    MemoryManager::SwitchProfile(std::move(old_prof_larger));
}


void Bootstrap_RangeCheck_PatersonStockmeyer(Ciphertext& ciphertext, const Ciphertext& input, const vector<uint64_t>& rangeCheckIndices,
                                             const int modulus, const size_t& degree, const RelinKeys &relin_keys, const SEALContext& context, const SecretKey& bfv_secret_key, 
                                             const int f_zero = 0, const bool gateEval = false, const bool skip_first_odd = true,
                                             const int firstDegree = 256, const int secondDegree = 256) {
    MemoryPoolHandle my_pool_larger = MemoryPoolHandle::New(true);
    auto old_prof_larger = MemoryManager::SwitchProfile(std::make_unique<MMProfFixed>(std::move(my_pool_larger)));

    Decryptor decryptor(context, bfv_secret_key);
    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);
    vector<Ciphertext> kCTs(firstDegree);
    // chrono::high_resolution_clock::time_point time_start, time_end;
    // time_start = chrono::high_resolution_clock::now();
    calUptoDegreeK(kCTs, input, firstDegree, relin_keys, context, skip_first_odd);
    // time_end = chrono::high_resolution_clock::now();
    // cout << "first level: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;

    // cout << "KCT :\n";
    // Plaintext pl;
    // vector<uint64_t> v(degree);
    // for (int i = 0; i < kCTs.size(); i++) {
    //     decryptor.decrypt(kCTs[i], pl);
    //     batch_encoder.decode(pl, v);
    //     cout << v[1] << ",";
    // }
    // cout << endl;

    vector<Ciphertext> kToMCTs(secondDegree);
    // time_start = chrono::high_resolution_clock::now();
    calUptoDegreeK(kToMCTs, kCTs[kCTs.size()-1], secondDegree, relin_keys, context);
    // time_end = chrono::high_resolution_clock::now();
    // cout << "second level: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;

    // cout << "KToMCT :\n";
    // for (int i = 0; i < kToMCTs.size(); i++) {
    //     decryptor.decrypt(kToMCTs[i], pl);
    //     batch_encoder.decode(pl, v);
    //     cout << v[1] << ",";
    // }
    // cout << endl;
    // decryptor.decrypt(kToMCTs[kToMCTs.size()-1], pl);
    // batch_encoder.decode(pl, v);
    // cout << v[1];

    cout << "Noise of last: " << decryptor.invariant_noise_budget(kToMCTs[kToMCTs.size()-1]) << " bits\n";
    cout << endl;


    // decryptor.decrypt(kToMCTs[kToMCTs.size()-1], pl);
    // batch_encoder.decode(pl, v);
    // cout << "Decrypt SHOULD BE ALL 1: \n" << v << endl;


    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < (int) kCTs.size(); j++) {
            evaluator.mod_switch_to_next_inplace(kCTs[j]);
        }
    }

    Ciphertext temp_relin(context);
    // int third = 0;
    
    for(int i = 0; i < secondDegree; i++) {
        Ciphertext levelSum;
        bool flag = false;
        for(int j = 0; j < firstDegree; j++) {
            if(rangeCheckIndices[i*firstDegree+j] != 0) {
                vector<uint64_t> intInd(degree, rangeCheckIndices[i*firstDegree+j]);
            // if (true) {
            //     vector<uint64_t> intInd(degree, 1);
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
        // time_start = chrono::high_resolution_clock::now();
        if(i == 0) {
            ciphertext = levelSum;
        } else if (i == 1) { // initialize for temp_relin, which is of ct size = 3
            evaluator.multiply(levelSum, kToMCTs[i - 1], temp_relin);
        } else {
            evaluator.multiply_inplace(levelSum, kToMCTs[i - 1]);
            evaluator.add_inplace(temp_relin, levelSum);
        }
        // time_end = chrono::high_resolution_clock::now();
        // third += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    }

    // cout << "third level: " << third << endl;
    
    for(int i = 0; i < firstDegree; i++){
        kCTs[i].release();
    }
    for(int i = 0; i < secondDegree; i++){
        kToMCTs[i].release();
    }

    evaluator.relinearize_inplace(temp_relin, relin_keys);
    evaluator.add_inplace(ciphertext, temp_relin);
    temp_relin.release();

    vector<uint64_t> intInd(degree, f_zero); 
    Plaintext plainInd;
    batch_encoder.encode(intInd, plainInd);
    evaluator.negate_inplace(ciphertext);
    evaluator.add_plain_inplace(ciphertext, plainInd);

    // decryptor.decrypt(ciphertext, pl);
    // batch_encoder.decode(pl, v);
    // cout << "Decrypt of function should be all ZERO: \n" << v << endl;
    cout << "Noise after function: " << decryptor.invariant_noise_budget(ciphertext) << " bits\n";

    if (gateEval) { // flip 0 to q/3, q/3 to 0
        vector<uint64_t> flip_constant(degree, modulus/3); 
        Plaintext flip_pl;
        batch_encoder.encode(flip_constant, flip_pl);
        evaluator.negate_inplace(ciphertext);
        evaluator.add_plain_inplace(ciphertext, flip_pl);
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


vector<regevCiphertext> extractRLWECiphertextToLWECiphertext(Ciphertext& rlwe_ct, const int ring_dim = poly_modulus_degree_glb,
                                                             const int n = 1024, const int p = 65537, const int big_prime = 268369921) {
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


vector<regevCiphertext> extractRLWECiphertextToLWECiphertext_NoModSwitch(Ciphertext& rlwe_ct, const int ring_dim = poly_modulus_degree_glb,
                                                             const int n = 1024, const int big_prime = 268369921) {
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
            results[cnt].a[ind] = rlwe_ct.data(1)[i];
            ind++;
        }

        for (int i = ring_dim-1; i > ring_dim - n + cnt && ind < n; i--) {
            results[cnt].a[ind] = -rlwe_ct.data(1)[i];
            ind++;
        }
        results[cnt].b = rlwe_ct.data(0)[cnt];
    }

    return results;
}


vector<regevCiphertext> preprocess_NAND(const vector<regevCiphertext>& ct_list_1, const vector<regevCiphertext>& ct_list_2,
                                        const regevParam& params) {

    chrono::high_resolution_clock::time_point time_start, time_end;
    time_start = chrono::high_resolution_clock::now();

    vector<regevCiphertext> result(ct_list_1.size());
    for (int i = 0; i < (int) ct_list_1.size(); i++) {
        result[i].a = NativeVector(params.n);
        for (int j = 0; j < params.n; j++) {
            result[i].a[j] = (ct_list_1[i].a[j].ConvertToInt() + ct_list_2[i].a[j].ConvertToInt()) % params.q;
        }
        result[i].b = (ct_list_1[i].b.ConvertToInt() + ct_list_2[i].b.ConvertToInt()) % params.q;
    }

    time_end = chrono::high_resolution_clock::now();
    cout << "TOTAL prepare NAND input time: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
    
    return result;
}


vector<regevCiphertext> bootstrap_bigPrime(vector<regevCiphertext>& lwe_ct_list, Ciphertext& lwe_sk_encrypted, const SEALContext& seal_context,
                                           const RelinKeys& relin_keys, const GaloisKeys& gal_keys, const int ring_dim, const int n,
                                           const int p, const KSwitchKeys& ksk, const vector<uint64_t>& rangeCheckIndices,
                                           const MemoryPoolHandle& my_pool, const SecretKey& bfv_secret_key, const vector<uint64_t>& q_shift_constant,
                                           const int f_zero = 0, const bool gateEval = false, const bool skip_first_odd = true,
                                           const int firstDegree = 256, const int secondDegree = 256) {
    chrono::high_resolution_clock::time_point time_start, time_end;
    int total_preprocess = 0, total_online = 0;

    Evaluator evaluator(seal_context);
    BatchEncoder batch_encoder(seal_context);
    Decryptor decryptor(seal_context, bfv_secret_key);

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
    Ciphertext result = evaluatePackedLWECiphertext(seal_context, lwe_ct_list, lwe_sk_sqrt_list, gal_keys, n, q_shift_constant, ring_dim, gateEval);
    time_end = chrono::high_resolution_clock::now();
    total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    cout << "TOTAL TIME for evaluation: " << total_online << endl;
    cout << "Noise: " << decryptor.invariant_noise_budget(result) << " bits\n";

    evaluator.mod_switch_to_next_inplace(result);
    cout << "Noise after mod down: " << decryptor.invariant_noise_budget(result) << " bits\n";

    Plaintext pl;
    vector<uint64_t> v(ring_dim);
    // decryptor.decrypt(result, pl);
    // batch_encoder.decode(pl, v);
    // cout << "Decrypt after evaluation: \n" << v << endl;

    Ciphertext range_check_res;
    time_start = chrono::high_resolution_clock::now();
    Bootstrap_RangeCheck_PatersonStockmeyer_bigPrime(range_check_res, result, rangeCheckIndices, p, ring_dim,
                                                     relin_keys, seal_context, bfv_secret_key, f_zero, gateEval, skip_first_odd,
                                                     firstDegree, secondDegree);
    time_end = chrono::high_resolution_clock::now();
    total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    cout << "TOTAL TIME for rangecheck: " << total_online << endl;

    decryptor.decrypt(range_check_res, pl);
    batch_encoder.decode(pl, v);
    cout << "Decrypt after rangeCheck: \n" << v << endl;

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
    // Ciphertext coeff = slotToCoeff(seal_context, ct_sqrt_list, U_plain_list, gal_keys, ring_dim);
    Ciphertext coeff = slotToCoeff_WOPrepreocess(seal_context, ct_sqrt_list, gal_keys, ring_dim);
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

vector<regevCiphertext> bootstrap(vector<regevCiphertext>& lwe_ct_list, Ciphertext& lwe_sk_encrypted, const SEALContext& seal_context,
                                  const RelinKeys& relin_keys, const GaloisKeys& gal_keys, const int ring_dim, const int n,
                                  const int p, const KSwitchKeys& ksk, const vector<uint64_t>& rangeCheckIndices,
                                  const MemoryPoolHandle& my_pool, const SecretKey& bfv_secret_key, const vector<uint64_t>& q_shift_constant,
                                  const int f_zero = 0, const bool gateEval = false, const bool skip_first_odd = true,
                                  const int firstDegree = 256, const int secondDegree = 256) {
    chrono::high_resolution_clock::time_point time_start, time_end;
    int total_preprocess = 0, total_online = 0;

    Evaluator evaluator(seal_context);
    BatchEncoder batch_encoder(seal_context);
    Decryptor decryptor(seal_context, bfv_secret_key);

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
    Ciphertext result = evaluatePackedLWECiphertext(seal_context, lwe_ct_list, lwe_sk_sqrt_list, gal_keys, n, q_shift_constant, ring_dim, gateEval);
    time_end = chrono::high_resolution_clock::now();
    total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    cout << "TOTAL TIME for evaluation: " << total_online << endl;
    // cout << "Noise: " << decryptor.invariant_noise_budget(result) << " bits\n";

    Plaintext pl;
    vector<uint64_t> v(ring_dim);
    // decryptor.decrypt(result, pl);
    // batch_encoder.decode(pl, v);
    // cout << "Decrypt after evaluation: \n" << v << endl;

    Ciphertext range_check_res;
    time_start = chrono::high_resolution_clock::now();
    Bootstrap_RangeCheck_PatersonStockmeyer(range_check_res, result, rangeCheckIndices, p, ring_dim,
                                            relin_keys, seal_context, bfv_secret_key, f_zero, gateEval, skip_first_odd,
                                            firstDegree, secondDegree);
    time_end = chrono::high_resolution_clock::now();
    total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    cout << "TOTAL TIME for rangecheck: " << total_online << endl;

    decryptor.decrypt(range_check_res, pl);
    batch_encoder.decode(pl, v);
    cout << "Decrypt after rangeCheck: \n" << v << endl;

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
    // Ciphertext coeff = slotToCoeff_WOPrepreocess(seal_context, ct_sqrt_list, gal_keys, ring_dim);
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