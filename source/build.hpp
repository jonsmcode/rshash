

inline uint64_t mark_sequences(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &input, const size_t k,
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> &endpoints)
{
    size_t length = 0;
    uint64_t kmers = 0;
    uint64_t no_sequences = 0;
    for(auto & record : input) {
        length += record.size();
        kmers += record.size() - k + 1;
        no_sequences++;
    }

    std::cout << "text length: " << length << "\n";
    std::cout << "text kmers: " << kmers <<  '\n';
    std::cout << "no sequences: " << no_sequences << "\n";

    std::cout << "mark endpoints BV...\n";
    bit_vector sequences = bit_vector(length+33, 0);
    sequences[0] = 1;
    sequences[32] = 1;
    uint64_t j = 32;
    for(uint64_t i=0; i < no_sequences; i++) {
        j += input[i].size();
        sequences[j] = 1;
    }
    endpoints = sux::bits::EliasFano(reinterpret_cast<uint64_t*>(sequences.data()), length+33);
    sequences = bit_vector();

    return kmers;
}


inline uint64_t get_unfrequent_minimizers(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    const size_t minimiser_length, const size_t threshold, const size_t k, const uint64_t seed,
    std::vector<uint64_t> &unfreq_minimizers, std::vector<uint8_t> &counts)
{
    std::vector<uint64_t> minimizers;

    std::cout << "computing minimizers...\n";
    auto minimiserview = rshash::views::xor_minimiser_and_positions({.minimiser_size = minimiser_length, .window_size = k, .seed=seed});
    for(auto & sequence : sequences) {
        for(auto && minimiser : sequence | minimiserview)
            minimizers.emplace_back(minimiser.minimiser_value);
    }

    std::cout << "sorting minimizers...\n";
    std::sort(minimizers.begin(), minimizers.end());

    std::cout << "get unfrequent minimizers...\n";
    uint64_t current_minimizer = minimizers[0];
    uint64_t occurences = 1;
    uint64_t no_skmers = 0;
    for(size_t i = 1; i < minimizers.size(); i++) {
        if(minimizers[i] != current_minimizer) {
            if(occurences <= threshold) {
                unfreq_minimizers.emplace_back(current_minimizer);
                counts.emplace_back(static_cast<uint8_t>(occurences));
                no_skmers += occurences;
            }
            current_minimizer = minimizers[i];
            occurences = 1;
        }
        else
            occurences++;
    }
    if(minimizers.back() == current_minimizer && occurences <= threshold) {
        unfreq_minimizers.emplace_back(current_minimizer);
        counts.emplace_back(static_cast<uint8_t>(occurences));
        no_skmers += occurences;
    }

    minimizers.clear();

    return no_skmers;
}



inline std::vector<uint64_t> pack_dna4_to_uint64(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &input)
{
    auto ranks = input | std::views::join | seqan3::views::to_rank;

    std::vector<uint64_t> packed;
    packed.push_back(UINT64_MAX); // padding

    uint64_t word = 0;
    size_t shift = 0;

    for (uint8_t r : ranks)
    {
        word |= uint64_t(r) << shift; // pack 2 bits per base
        shift += 2;

        if (shift == 64)
        {
            packed.push_back(word);
            word = 0;
            shift = 0;
        }
    }

    if (shift != 0)
        packed.push_back(word);

    packed.push_back(UINT64_MAX); // padding

    return packed;
}

