#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <thread>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <unistd.h> // POSIX 命令行解析

using namespace std;

// ============================================================================
// 全局配置与数据结构
// ============================================================================
struct Config {
    int k = 21;
    int max_freq = 4;
    int threshold = 15;      // LIS 阈值
    int threads = 4;
    string output_mode = "filter"; // "tag" 或 "filter"
    string aln_format = "";  // "sam" 或 "paf"
    string out_file = "";    // 输出文件

    string ref_fasta;
    string aln_file;         // "-" 代表 stdin 管道
    string query_fasta;      // PAF 模式下必须提供
};

unordered_set<uint64_t> rare_kmers;
unordered_map<string, string> ref_genomes;

struct CompressedSeq {
    int length;
    vector<uint8_t> data;
};
unordered_map<string, CompressedSeq> query_genomes;

// ============================================================================
// 底层算法：2-bit 序列编码与内存压缩
// ============================================================================
inline uint64_t encode_base(char c) {
    switch (c) {
        case 'A': case 'a': return 0ULL;
        case 'C': case 'c': return 1ULL;
        case 'G': case 'g': return 2ULL;
        case 'T': case 't': return 3ULL;
        default: return 4ULL; // N 碱基处理
    }
}

inline char decode_base(uint8_t val) {
    switch (val) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'N';
    }
}

CompressedSeq compress_sequence(const string& seq) {
    CompressedSeq cs;
    cs.length = seq.length();
    cs.data.reserve((seq.length() + 3) / 4);
    uint8_t current_byte = 0;
    for (int i = 0; i < seq.length(); ++i) {
        uint8_t val = encode_base(seq[i]);
        if (val > 3) val = 0; // 将 N 简化处理为 A，避免打断整体压缩
        current_byte = (current_byte << 2) | val;
        if ((i + 1) % 4 == 0) {
            cs.data.push_back(current_byte);
            current_byte = 0;
        }
    }
    if (seq.length() % 4 != 0) {
        current_byte <<= (2 * (4 - (seq.length() % 4))); 
        cs.data.push_back(current_byte);
    }
    return cs;
}

string decompress_sequence(const CompressedSeq& cs) {
    string seq(cs.length, 'A');
    int idx = 0;
    for (uint8_t byte : cs.data) {
        for (int i = 3; i >= 0 && idx < cs.length; --i) {
            seq[idx++] = decode_base((byte >> (i * 2)) & 0x03);
        }
    }
    return seq;
}

// ============================================================================
// 核心数学引擎：O(n log n) 求解最长递增子序列 (LIS)
// ============================================================================
int calculate_LIS(const vector<int>& arr) {
    if (arr.empty()) return 0;
    vector<int> tails;
    for (int x : arr) {
        auto it = lower_bound(tails.begin(), tails.end(), x);
        if (it == tails.end()) tails.push_back(x);
        else *it = x;
    }
    return tails.size();
}

// ============================================================================
// 模块 1：构建参考基因组的全局稀有 k-mer 字典
// ============================================================================
void build_rare_kmer_dictionary(const Config& cfg) {
    cerr << "[Info] Building Rare k-mer Dictionary from Reference..." << endl;
    ifstream fas(cfg.ref_fasta);
    if(!fas.is_open()) { cerr << "[Error] Cannot open ref fasta: " << cfg.ref_fasta << endl; exit(1); }
    
    string line, name, seq;
    unordered_map<uint64_t, int> kmer_counts;
    uint64_t mask = (1ULL << (2 * cfg.k)) - 1;

    auto process_ref = [&](const string& n, const string& s) {
        ref_genomes[n] = s;
        if (s.length() < cfg.k) return;
        uint64_t kmer = 0; int valid = 0;
        for (char base : s) {
            uint64_t val = encode_base(base);
            if (val > 3) { valid = 0; kmer = 0; continue; }
            kmer = ((kmer << 2) | val) & mask;
            valid++;
            if (valid >= cfg.k) kmer_counts[kmer]++;
        }
    };

    while (getline(fas, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!name.empty()) process_ref(name, seq);
            name = line.substr(1, line.find_first_of(" \t") - 1);
            seq.clear();
        } else seq += line;
    }
    if (!name.empty()) process_ref(name, seq);

    for (const auto& kv : kmer_counts) {
        if (kv.second <= cfg.max_freq) rare_kmers.insert(kv.first);
    }
    cerr << "[Info] Retained " << rare_kmers.size() << " rare k-mers globally." << endl;
}

// ============================================================================
// 模块 2：预载 Query Reads (仅限 PAF 模式使用)
// ============================================================================
void load_query_reads(const Config& cfg) {
    if (cfg.aln_format != "paf") return;
    cerr << "[Info] PAF Mode: Loading and compressing query reads into memory..." << endl;
    ifstream fas(cfg.query_fasta);
    if(!fas.is_open()) { cerr << "[Error] Cannot open query fasta: " << cfg.query_fasta << endl; exit(1); }
    
    string line, name, seq;
    auto store_query = [&](const string& n, const string& s) {
        query_genomes[n] = compress_sequence(s);
    };

    while (getline(fas, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!name.empty()) store_query(name, seq);
            name = line.substr(1, line.find_first_of(" \t") - 1);
            seq.clear();
        } else seq += line;
    }
    if (!name.empty()) store_query(name, seq);
    cerr << "[Info] Query reads loaded and compressed successfully." << endl;
}

// ============================================================================
// 模块 3：核心校验逻辑 (供多线程现场提取调用)
// ============================================================================
int validate_alignment(const string& ref_seq, const string& query_seq, int ref_start, const Config& cfg) {
    uint64_t mask = (1ULL << (2 * cfg.k)) - 1;

    // 1. 提取 Reference 区间的稀有 k-mer 并计数
    unordered_map<uint64_t, vector<int>> ref_kmers_pos;
    uint64_t kmer = 0; int valid = 0;
    int ref_rare_count = 0; // 记录 Ref 中的稀有 k-mer 总数

    for (int i = 0; i < ref_seq.length(); ++i) {
        uint64_t val = encode_base(ref_seq[i]);
        if (val > 3) { valid = 0; continue; }
        kmer = ((kmer << 2) | val) & mask;
        valid++;
        if (valid >= cfg.k && rare_kmers.count(kmer)) {
            ref_kmers_pos[kmer].push_back(i);
            ref_rare_count++;
        }
    }

    // 2. 提取 Query 区间的稀有 k-mer 并映射坐标
    vector<pair<int, int>> shared_pairs; // <ref_pos, query_pos>
    kmer = 0; valid = 0;
    int query_rare_count = 0; // 记录 Query 中的稀有 k-mer 总数

    for (int i = 0; i < query_seq.length(); ++i) {
        uint64_t val = encode_base(query_seq[i]);
        if (val > 3) { valid = 0; continue; }
        kmer = ((kmer << 2) | val) & mask;
        valid++;
        if (valid >= cfg.k && rare_kmers.count(kmer)) {
            query_rare_count++;
            if (ref_kmers_pos.count(kmer)) {
                shared_pairs.push_back({ref_kmers_pos[kmer][0], i});
            }
        }
    }

    // 如果该区域没有任何稀有 k-mer 作为锚点，返回极低分 0
    int min_kmers = min(ref_rare_count, query_rare_count);
    if (min_kmers == 0) return 255;

    // 3. 计算 LIS (最长公共递增序列)
    sort(shared_pairs.begin(), shared_pairs.end());
    vector<int> query_coords;
    for (auto& sp : shared_pairs) query_coords.push_back(sp.second);

    int lis_length = calculate_LIS(query_coords);

    // 4. 计算 Phred-scaled 共线性得分 (模拟 MAPQ)
    double ratio = (double)lis_length / min_kmers;

    // 保护机制：如果出现异常（如超过 1.0），强制归一化
    if (ratio >= 1.0) ratio = 1.0;

    double error_rate = 1.0 - ratio;

    // 限制最小错误率为 1e-6，防止 log10(0) 崩溃，同时设定最高分为 60 分
    if (error_rate < 1e-6) {
        error_rate = 1e-6;
    }

    double phred_score = -10.0 * log10(error_rate);

    // 返回四舍五入后的整数分数 (范围 0 ~ 60)
    return static_cast<int>(round(phred_score));
}
// ============================================================================
// 模块 4：处理单行对齐结果
// ============================================================================
string process_line(const string& line, const Config& cfg) {
    if (line.empty() || (cfg.aln_format == "sam" && line[0] == '@')) return line; // 保留头文件

    istringstream iss(line);
    string qname, rname, seq;
    int lis_score = 0;

    if (cfg.aln_format == "sam") {
        string flag, pos_str, mapq, cigar, rnext, pnext, tlen, qual;
        iss >> qname >> flag >> rname >> pos_str >> mapq >> cigar >> rnext >> pnext >> tlen >> seq >> qual;
        if (rname == "*" || ref_genomes.find(rname) == ref_genomes.end() || seq == "*") {
            return (cfg.output_mode == "tag") ? line + "\tkm:i:0" : "";
        }
        
        int pos = stoi(pos_str) - 1;
        if(pos < 0) pos = 0;
        // 越界保护
        int len = min((int)seq.length(), (int)ref_genomes[rname].length() - pos);
        if(len <= 0) return (cfg.output_mode == "tag") ? line + "\tkm:i:0" : "";
        
        string ref_sub = ref_genomes[rname].substr(pos, len);
        lis_score = validate_alignment(ref_sub, seq, pos, cfg);

    } else if (cfg.aln_format == "paf") {
        string qlen, qstart, qend, strand, rlen, rstart, rend;
        iss >> qname >> qlen >> qstart >> qend >> strand >> rname >> rlen >> rstart >> rend;
        if (ref_genomes.find(rname) == ref_genomes.end() || query_genomes.find(qname) == query_genomes.end()) {
            return (cfg.output_mode == "tag") ? line + "\tkm:i:0" : "";
        }
        
        int r_st = stoi(rstart);
		int r_ed = stoi(rend);
        // 越界保护
        r_ed = min(r_ed, (int)ref_genomes[rname].length());
        if(r_ed <= r_st) return (cfg.output_mode == "tag") ? line + "\tkm:i:0" : "";
        
        string ref_sub = ref_genomes[rname].substr(r_st, r_ed - r_st);
        string query_sub = decompress_sequence(query_genomes[qname]);
        lis_score = validate_alignment(ref_sub, query_sub, r_st, cfg);
    }

    if (cfg.output_mode == "filter") {
        return (lis_score >= cfg.threshold) ? line : "";
    } else {
        return line + "\tkm:i:" + to_string(lis_score);
    }
}

// ============================================================================
// 模块 5：流水线多线程调度器
// ============================================================================
void run_pipeline(const Config& cfg) {
    istream* in_stream;
    ifstream file_stream;
    if (cfg.aln_file == "-") {
        in_stream = &cin;
        cerr << "[Info] Input: Reading alignment stream from Standard Input (stdin)..." << endl;
    } else {
        file_stream.open(cfg.aln_file);
        if (!file_stream.is_open()) { cerr << "[Error] Cannot open input file.\n"; exit(1); }
        in_stream = &file_stream;
        cerr << "[Info] Input: Reading alignment file (" << cfg.aln_file << ")..." << endl;
    }

    ostream* out_stream = &cout;
    ofstream out_file_stream;
    if (!cfg.out_file.empty()) {
        out_file_stream.open(cfg.out_file);
        if (!out_file_stream.is_open()) { cerr << "[Error] Cannot open output file.\n"; exit(1); }
        out_stream = &out_file_stream;
        cerr << "[Info] Output: Writing directly to disk (" << cfg.out_file << ")..." << endl;
    }

    string line;
    const int CHUNK_SIZE = 10000;
    vector<string> lines_chunk, results_chunk;

    auto worker = [&](int start, int end) {
        for (int i = start; i < end; ++i) {
            results_chunk[i] = process_line(lines_chunk[i], cfg);
        }
    };

    while (true) {
        lines_chunk.clear();
        while (lines_chunk.size() < CHUNK_SIZE && getline(*in_stream, line)) {
            lines_chunk.push_back(line);
        }
        if (lines_chunk.empty()) break;

        results_chunk.assign(lines_chunk.size(), "");
        vector<thread> threads;
        int chunk_per_thread = lines_chunk.size() / cfg.threads;
        int start = 0;
        
        for (int i = 0; i < cfg.threads; ++i) {
            int end = (i == cfg.threads - 1) ? lines_chunk.size() : start + chunk_per_thread;
            if (start < end) threads.emplace_back(worker, start, end);
            start = end;
        }
        for (auto& t : threads) t.join();

        for (const string& res : results_chunk) {
            if (!res.empty()) (*out_stream) << res << "\n";
        }
    }
}

// ============================================================================
// 命令行帮助与入口函数
// ============================================================================
void print_usage(const char* prog_name) {
    cerr << "========================================================================\n"
         << " rafilter: Rare-kmer based Alignment Filter using O(n log n) LIS engine\n"
         << "========================================================================\n"
         << "Usage: " << prog_name << " [options] <ref.fa> <in.aln> [query.fa]\n\n"
         << "Options:\n"
         << "  -k INT     k-mer length [default: 21]\n"
         << "  -n INT     max frequency for rare k-mers [default: 4]\n"
         << "  -c INT     LIS score cutoff for filtering [default: 15]\n"
         << "  -t INT     number of threads to use [default: 4]\n"
         << "  -m STR     output mode: 'filter' or 'tag' [default: filter]\n"
         << "  -f STR     input format: 'sam' or 'paf' [default: auto-detect]\n"
         << "  -o FILE    output file name [default: stdout]\n"
         << "  -h         show this help message and exit\n\n"
         << "Positional Arguments:\n"
         << "  <ref.fa>   Reference genome in FASTA format\n"
         << "  <in.aln>   Alignment file (use '-' for reading from stdin pipe)\n"
         << "  [query.fa] Required ONLY if format is 'paf' to supply read sequences\n\n"
         << "Example:\n"
         << "  samtools view -h in.bam | " << prog_name << " -k 21 -n 4 -t 8 -o out.sam ref.fa -\n"
         << "========================================================================\n";
}

int main(int argc, char* argv[]) {
    // 开启 C++ 极速 I/O (关闭与 C stdio 的同步，极大加快管道吞吐率)
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    Config cfg;
    int opt;
    
    while ((opt = getopt(argc, argv, "k:n:c:t:m:f:o:h")) != -1) {
        switch (opt) {
            case 'k': cfg.k = stoi(optarg); break;
            case 'n': cfg.max_freq = stoi(optarg); break;
            case 'c': cfg.threshold = stoi(optarg); break;
            case 't': cfg.threads = stoi(optarg); break;
            case 'm': cfg.output_mode = optarg; break;
            case 'f': cfg.aln_format = optarg; break;
            case 'o': cfg.out_file = optarg; break;
            case 'h': print_usage(argv[0]); return 0;
            default:  print_usage(argv[0]); return 1;
        }
    }

    if (optind + 2 > argc) {
        cerr << "[Error] Missing required positional arguments.\n\n";
        print_usage(argv[0]);
        return 1;
    }

    cfg.ref_fasta = argv[optind];
    cfg.aln_file = argv[optind + 1];

    if (optind + 2 < argc) {
        cfg.query_fasta = argv[optind + 2];
    }

    // 自动检测文件格式
    if (cfg.aln_format.empty()) {
        if (cfg.aln_file.length() >= 4 && cfg.aln_file.substr(cfg.aln_file.length() - 4) == ".paf") {
            cfg.aln_format = "paf";
        } else {
            cfg.aln_format = "sam"; 
        }
    }

    if (cfg.aln_format == "paf" && cfg.query_fasta.empty()) {
        cerr << "[Error] PAF format explicitly requires [query.fa] as the last argument.\n";
        return 1;
    }

    cerr << "[Info] Starting rafilter with parameters:\n"
         << "       k-mer length: " << cfg.k << "\n"
         << "       max rare freq: " << cfg.max_freq << "\n"
         << "       threads: " << cfg.threads << "\n"
         << "       mode: " << cfg.output_mode << "\n"
         << "       format: " << cfg.aln_format << "\n"
         << "       cutoff: " << cfg.threshold << "\n"
         << "--------------------------------------------------------\n";

    build_rare_kmer_dictionary(cfg);
    load_query_reads(cfg);
    run_pipeline(cfg);

    cerr << "[Info] rafilter pipeline completed successfully." << endl;
    return 0;
}
