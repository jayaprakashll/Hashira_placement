#include <bits/stdc++.h>
using namespace std;
// ---------- small helpers ----------
static inline bool is_ws(char c){ return c==' '||c=='\n'||c=='\t'||c=='\r'; }
static inline bool is_digit(char c){ return c>='0' && c<='9'; }
static inline int  digval(char c){
    if (c>='0' && c<='9') return c-'0';
    if (c>='a' && c<='z') return 10 + (c-'a');
    if (c>='A' && c<='Z') return 10 + (c-'A');
    return -1;
}

// ---------- simple JSON-ish parser (just enough for the provided format) ----------
struct Share { long long x; int base; string val; }; // val may be very long
struct Input {
    int n=0, k=0;
    vector<Share> shares;
};

// Very minimal: extract "n", "k", then each numeric string key -> { "base": "...", "value": "..." }
Input parse_input(const string& s) {
    Input R;
    // Find "n": <int>, "k": <int>
    auto find_num_after_key = [&](const string& key)->long long{
        size_t p = s.find("\""+key+"\"");
        if (p==string::npos) return 0;
        p = s.find(':', p);
        if (p==string::npos) return 0;
        ++p;
        while (p<s.size() && is_ws(s[p])) ++p;
        // read number
        long long val=0; bool neg=false;
        if (p<s.size() && s[p]=='-'){ neg=true; ++p; }
        while (p<s.size() && is_digit(s[p])){ val = val*10 + (s[p]-'0'); ++p; }
        return neg ? -val : val;
    };
    R.n = (int)find_num_after_key("n");
    R.k = (int)find_num_after_key("k");

    // Collect all top-level keys that look like:  "123": { "base": "...", "value":"..." }
    // We'll scan for quotes followed by digits.
    size_t i=0;
    while (true) {
        size_t q = s.find('"', i);
        if (q==string::npos) break;
        size_t q2 = s.find('"', q+1);
        if (q2==string::npos) break;
        string key = s.substr(q+1, q2-q-1);
        i = q2+1;

        // numeric key?
        bool allnum = !key.empty();
        for (char c: key) if (!is_digit(c)) { allnum=false; break; }
        if (!allnum) continue;

        // find object after :
        size_t col = s.find(':', i);
        if (col==string::npos) break;
        // move to '{'
        size_t br = s.find('{', col);
        if (br==string::npos) continue;

        // find matching '}' (shallow)
        int depth=0; size_t j=br;
        for (; j<s.size(); ++j) {
            if (s[j]=='{') depth++;
            else if (s[j]=='}') { depth--; if (depth==0){ ++j; break; } }
        }
        if (j>=s.size()) break;
        string obj = s.substr(br, j-br);

        // inside obj, find "base": "..." and "value": "..."
        auto find_str_after_key = [&](const string& key2)->string{
            size_t p = obj.find("\""+key2+"\"");
            if (p==string::npos) return "";
            p = obj.find(':', p);
            if (p==string::npos) return "";
            p++;
            while (p<obj.size() && is_ws(obj[p])) ++p;
            if (p<obj.size() && obj[p]=='"') {
                size_t e = obj.find('"', p+1);
                if (e!=string::npos) return obj.substr(p+1, e-(p+1));
            }
            // also allow unquoted numbers (not expected here)
            size_t st=p;
            while (st<obj.size() && (is_ws(obj[st])||obj[st]==':')) ++st;
            size_t en=st;
            while (en<obj.size() && (is_digit(obj[en])||obj[en]=='-')) ++en;
            return obj.substr(st, en-st);
        };

        string bstr = find_str_after_key("base");
        string vstr = find_str_after_key("value");
        if (!bstr.empty() && !vstr.empty()) {
            Share sh;
            sh.x = stoll(key);
            sh.base = stoi(bstr);
            sh.val = vstr;
            R.shares.push_back(sh);
        }
        i = j;
    }
    sort(R.shares.begin(), R.shares.end(), [](auto& a, auto& b){ return a.x<b.x; });
    return R;
}

// ---------- modular arithmetic ----------
using u64 = unsigned long long;
using i64 = long long;
using u128 = __uint128_t;

u64 mod_add(u64 a, u64 b, u64 mod){ a%=mod; b%=mod; u64 r=a+b; if (r>=mod) r-=mod; return r; }
u64 mod_sub(u64 a, u64 b, u64 mod){ a%=mod; b%=mod; return (a>=b)?(a-b):(a + (mod-b)); }
u64 mod_mul(u64 a, u64 b, u64 mod){ return (u128)a*b % mod; }

u64 mod_pow(u64 a, u64 e, u64 mod){
    u64 r=1%mod; a%=mod;
    while (e){
        if (e&1) r = mod_mul(r,a,mod);
        a = mod_mul(a,a,mod);
        e >>= 1;
    }
    return r;
}
u64 mod_inv(u64 a, u64 mod){ // mod is prime
    return mod_pow(a, mod-2, mod);
}

// ---------- convert huge base-B string to residues modulo several primes ----------
static const u64 PRIMES[] = {
    1000000007ULL, 1000000009ULL, 1000000033ULL, 1000000087ULL,
    1000000093ULL, 1000000097ULL, 1000000103ULL, 1000000123ULL
};
static const int NP = sizeof(PRIMES)/sizeof(PRIMES[0]);

struct ModY {
    u64 r[8];
};

ModY value_mod_primes(const string& s, int base) {
    ModY out{};
    for (int i=0;i<NP;i++) out.r[i]=0;
    for (char c: s){
        int v = digval(c);
        if (v<0 || v>=base){
            // ignore (robustness) – or treat as error in production
            continue;
        }
        for (int i=0;i<NP;i++){
            u64 p = PRIMES[i];
            out.r[i] = ( (out.r[i]*base)%p + (u64)v ) % p;
        }
    }
    return out;
}

// ---------- Lagrange interpolation at x=0 modulo p ----------
u64 lagrange_secret_at_zero(const vector<pair<u64,u64>>& pts, u64 mod){
    // pts.size() == k; x distinct modulo mod
    int k = (int)pts.size();
    u64 S = 0;
    for (int i=0;i<k;i++){
        u64 xi = pts[i].first % mod;
        u64 yi = pts[i].second % mod;
        u64 num = 1, den = 1;
        for (int j=0;j<k;j++){
            if (i==j) continue;
            u64 xj = pts[j].first % mod;
            // l_i(0) = prod_{j!=i} (-xj)/(xi-xj)
            num = mod_mul(num, (mod - xj) % mod, mod);
            u64 denom = (xi + mod - xj) % mod;
            den = mod_mul(den, denom, mod);
        }
        u64 li0 = mod_mul(num, mod_inv(den, mod), mod);
        u64 term = mod_mul(yi, li0, mod);
        S = mod_add(S, term, mod);
    }
    return S; // f(0)
}

// ---------- Fit polynomial coefficients (degree k-1) modulo p via Vandermonde ----------
vector<u64> fit_coeffs(const vector<pair<u64,u64>>& pts, u64 mod){
    int k = (int)pts.size();
    // A[i][j] = x_i^j, solve A a = y
    vector<vector<u64>> A(k, vector<u64>(k,0));
    vector<u64> y(k,0);
    for (int i=0;i<k;i++){
        u64 xi = pts[i].first % mod;
        u64 pwr = 1;
        for (int j=0;j<k;j++){
            A[i][j] = pwr;
            pwr = mod_mul(pwr, xi, mod);
        }
        y[i] = pts[i].second % mod;
    }
    // Gaussian elimination mod prime
    for (int col=0, row=0; col<k && row<k; col++, row++){
        int piv = row;
        for (int r=row;r<k;r++) if (A[r][col]!=0){piv=r;break;}
        if (A[piv][col]==0){ continue; } // singular in this modulus; unlikely with good input
        if (piv!=row){ swap(A[piv],A[row]); swap(y[piv],y[row]); }
        u64 inv = mod_inv(A[row][col], mod);
        for (int j=col;j<k;j++) A[row][j] = mod_mul(A[row][j], inv, mod);
        y[row] = mod_mul(y[row], inv, mod);
        for (int r=0;r<k;r++){
            if (r==row) continue;
            u64 f = A[r][col];
            if (!f) continue;
            for (int j=col;j<k;j++){
                u64 sub = mod_mul(f, A[row][j], mod);
                A[r][j] = mod_sub(A[r][j], sub, mod);
            }
            u64 ysub = mod_mul(f, y[row], mod);
            y[r] = mod_sub(y[r], ysub, mod);
        }
    }
    // y now holds the solution a0..a_{k-1}
    return y;
}

u64 eval_poly(const vector<u64>& a, u64 x, u64 mod){
    u64 s=0, p=1%mod;
    for (u64 c: a){
        s = mod_add(s, mod_mul(c,p,mod), mod);
        p = mod_mul(p, x, mod);
    }
    return s;
}

// ---------- Garner CRT to decimal string ----------
// We get residues R[i] modulo M[i] (PRIMES). Compute x uniquely modulo product.
// Garner gives x = c0 + c1*M0 + c2*M0*M1 + ...
// Then we build a decimal big-int by repeated (big * small) + (big * small).
struct Big {
    // store as base 1e9
    static const uint32_t BASE = 1000000000u;
    vector<uint32_t> d; // little-endian
    Big(){ d.clear(); }
    bool isZero() const { return d.empty(); }
    void trim(){ while(!d.empty() && d.back()==0) d.pop_back(); }

    static Big from_u64(u64 x){
        Big b;
        while (x){ b.d.push_back((uint32_t)(x % BASE)); x /= BASE; }
        return b;
    }
    static Big mul_small(const Big& a, u64 m){
        if (a.d.empty() || m==0) return Big();
        Big r; r.d.resize(a.d.size());
        unsigned long long carry=0;
        for (size_t i=0;i<a.d.size();i++){
            unsigned long long cur = (unsigned long long)a.d[i]*m + carry;
            r.d[i] = (uint32_t)(cur % BASE);
            carry = cur / BASE;
        }
        while (carry){
            r.d.push_back((uint32_t)(carry % BASE));
            carry /= BASE;
        }
        return r;
    }
    static Big add(const Big& a, const Big& b){
        Big r;
        size_t n = max(a.d.size(), b.d.size());
        r.d.resize(n);
        unsigned long long carry=0;
        for (size_t i=0;i<n;i++){
            unsigned long long cur = carry;
            if (i<a.d.size()) cur += a.d[i];
            if (i<b.d.size()) cur += b.d[i];
            r.d[i] = (uint32_t)(cur % BASE);
            carry = cur / BASE;
        }
        while (carry){
            r.d.push_back((uint32_t)(carry % BASE));
            carry /= BASE;
        }
        r.trim();
        return r;
    }
    string to_string() const {
        if (d.empty()) return "0";
        stringstream ss;
        int n = (int)d.size();
        ss << d.back();
        for (int i=n-2;i>=0;i--){
            ss << setw(9) << setfill('0') << d[i];
        }
        return ss.str();
    }
};

string garner_to_decimal(const array<u64,8>& R, int np){
    // precompute inverses invM[i][j] s.t. invM[i][j] = (M[i])^{-1} mod M[j], i<j
    // Garner core:
    vector<u64> c(np,0);
    vector<u64> m(PRIMES, PRIMES+np);
    // For speed, precompute inverses
    vector<vector<u64>> inv(np, vector<u64>(np, 0));
    for (int i=0;i<np;i++){
        for (int j=i+1;j<np;j++){
            inv[i][j] = mod_inv(m[i]%m[j], m[j]);
        }
    }
    for (int i=0;i<np;i++){
        u64 t = R[i];
        for (int j=0;j<i;j++){
            // t = (t - c[j]) * inv(prod_{<j}) mod m[i]
            t = mod_sub(t, c[j], m[i]);
            t = mod_mul(t, inv[j][i], m[i]);
        }
        c[i] = t; // < m[i]
    }
    // Build decimal: x = sum c[i] * (M0*M1*...*M_{i-1})
    Big res; // 0
    Big pref; // 1
    pref.d = {1};
    for (int i=0;i<np;i++){
        Big term = Big::mul_small(pref, c[i]);
        res = Big::add(res, term);
        pref = Big::mul_small(pref, m[i]);
    }
    return res.to_string();
}

// ---------- main solve for one JSON ----------
struct SolveResult {
    string secret_dec;             // decimal
    vector<long long> wrong_xs;    // x-values that are wrong
};

SolveResult solve_one(const string& jsonText){
    Input in = parse_input(jsonText);
    int n = (int)in.shares.size();
    int k = in.k;

    // Precompute y residues for each share
    struct PShare {
        u64 x;
        array<u64,8> ymod;
    };
    vector<PShare> pts(n);
    for (int i=0;i<n;i++){
        pts[i].x = (u64)in.shares[i].x;
        ModY my = value_mod_primes(in.shares[i].val, in.shares[i].base);
        for (int t=0;t<NP;t++) pts[i].ymod[t] = my.r[t];
    }

    // 1) Compute secret residues for ALL k-subsets; choose majority secret vector
    //    (We key by the vector<residuals> across primes)
    struct VecHash {
        size_t operator()(const vector<u64>& v) const {
            size_t h=0;
            for (u64 x: v) h = (h*1315423911u) ^ (x+0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2));
            return h;
        }
    };
    unordered_map<vector<u64>, int, VecHash> cnt;
    vector<vector<int>> subsets_idx; // store idxs for one representative per secret if needed
    vector<vector<u64>> secret_vecs;

    vector<int> idx(n);
    iota(idx.begin(), idx.end(), 0);

    // To keep time bounded, we’ll pick a few random k-subsets if n is large,
    // but here n<=10 so iterate all.
    vector<int> bitmask(n,0);
    fill(bitmask.begin(), bitmask.begin()+k, 1);
    // generate combinations
    vector<int> indices;
    do{
        indices.clear();
        for (int i=0;i<n;i++) if (bitmask[i]) indices.push_back(i);
        vector<u64> sres(NP,0);
        for (int t=0;t<NP;t++){
            vector<pair<u64,u64>> sub;
            sub.reserve(k);
            for (int id: indices) sub.emplace_back(pts[id].x, pts[id].ymod[t]);
            sres[t] = lagrange_secret_at_zero(sub, PRIMES[t]);
        }
        cnt[sres]++;
        secret_vecs.push_back(sres);
        subsets_idx.push_back(indices);
    } while (prev_permutation(bitmask.begin(), bitmask.end()));

    // pick majority secret vector
    int bestc=-1, besti=-1;
    for (int i=0;i<(int)secret_vecs.size();i++){
        int ccur = cnt[secret_vecs[i]];
        if (ccur>bestc){ bestc=ccur; besti=i; }
    }
    vector<u64> best_secret_res = secret_vecs[besti];

    // 2) Using that representative subset, fit a polynomial modulo each prime,
    //    then check all points for consistency.
    vector<int> good_subset = subsets_idx[besti];
    vector<vector<u64>> coeffs_mod(NP);
    vector<pair<u64,u64>> subpts;
    subpts.reserve(k);
    for (int t=0;t<NP;t++){
        subpts.clear();
        for (int id: good_subset) subpts.emplace_back(pts[id].x, pts[id].ymod[t]);
        coeffs_mod[t] = fit_coeffs(subpts, PRIMES[t]);
    }

    auto point_is_ok = [&](int i)->bool{
        for (int t=0;t<NP;t++){
            u64 expect = eval_poly(coeffs_mod[t], pts[i].x % PRIMES[t], PRIMES[t]);
            if (expect != pts[i].ymod[t]) return false;
        }
        return true;
    };

    vector<long long> wrong;
    for (int i=0;i<n;i++){
        if (!point_is_ok(i)) wrong.push_back((long long)pts[i].x);
    }

    // 3) Garner CRT to decimal for the secret
    array<u64,8> R{};
    for (int t=0;t<NP;t++) R[t] = best_secret_res[t];
    string secret_dec = garner_to_decimal(R, NP);

    return {secret_dec, wrong};
}

// ---------- driver ----------
int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    string jsonText, line;
    while (getline(cin, line)){
        jsonText += line;
        jsonText.push_back('\n');
    }
    if (jsonText.find("\"keys\"")==string::npos){
        cerr<<"No JSON detected on stdin.\n";
        return 0;
    }
    // You can paste the first JSON then the second (run twice), or keep them in separate files.
    // For convenience, we attempt to split if two JSON blobs are concatenated.
    // Simple split on a '}\n{' boundary.
    vector<string> blobs;
    {
        // naive split
        int depth=0;
        size_t start=0;
        for (size_t i=0;i<jsonText.size();i++){
            if (jsonText[i]=='{') depth++;
            else if (jsonText[i]=='}'){
                depth--;
                if (depth==0){
                    blobs.push_back(jsonText.substr(start, i-start+1));
                    // skip whitespace
                    size_t j=i+1;
                    while (j<jsonText.size() && is_ws(jsonText[j])) j++;
                    start=j;
                }
            }
        }
        if (blobs.empty()) blobs.push_back(jsonText);
    }

    for (int bi=0; bi<(int)blobs.size(); ++bi){
        SolveResult R = solve_one(blobs[bi]);
        cout << "TestCase-" << (bi+1) << " Secret: " << R.secret_dec << "\n";
        if (R.wrong_xs.empty()){
            cout << "Wrong points: None\n";
        } else {
            cout << "Wrong points at x =";
            for (auto x: R.wrong_xs) cout << " " << x;
            cout << "\n";
        }
        if (bi+1 < (int)blobs.size()) cout << "-----\n";
    }
    return 0;
}
