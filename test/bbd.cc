#include <sndfile.h>
#include <getopt.h>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "../bbd_filter.cc"
#include "../bbd_line.cc"

static constexpr unsigned default_bucket_count = 256;
static constexpr float default_clock = 50000;

///
static void print_usage();
static std::vector<float> read_mono_samples(SNDFILE *sf, const SF_INFO &info);

int main(int argc, char *argv[])
{
    const char *output_path = nullptr;
    unsigned bucket_count = default_bucket_count;
    float fbbd = default_clock;

    if (argc <= 1) {
        print_usage();
        return 0;
    }

    for (int c; (c = getopt(argc, argv, "o:b:c:h")) != -1;) {
        switch (c) {
        case 'o':
            output_path = optarg;
            break;
        case 'b':
            bucket_count = atoi(optarg);
            if ((int)bucket_count < 16 || bucket_count > 8192) {
                fprintf(stderr, "Invalid value for bucket count.\n");
                return 1;
            }
            break;
        case 'c':
            fbbd = atof(optarg);
            if (fbbd <= 0) {
                fprintf(stderr, "Invalid value for BBD clock.\n");
                return 1;
            }
            break;
        case 'h':
            print_usage();
            return 0;
        default:
            print_usage();
            return 1;
        }
    }

    if (argc - optind != 1) {
        fprintf(stderr, "No file has been specified for input.\n");
        return 1;
    }

    const char *input_path = argv[optind];

    //
    SF_INFO info_in;
    SNDFILE *sf_in = sf_open(input_path, SFM_READ, &info_in);

    if (!sf_in) {
        fprintf(stderr, "Cannot open the sound file for reading.\n");
        return 1;
    }

    //
    std::vector<float> input = read_mono_samples(sf_in, info_in);
    sf_close(sf_in);

    //
    float fs = info_in.samplerate;
    unsigned ns = (unsigned)input.size();

    fprintf(stderr, "* Sample rate: %f\n", fs);
    fprintf(stderr, "* Frame count: %u\n", ns);
    fprintf(stderr, "* File channels: %u\n", info_in.channels);

    BBD_Line line(fs, ns, bbd_fin_j60, bbd_fout_j60);
    line.set_delay_size(256);

    std::unique_ptr<float[]> output(new float[ns]);
    std::unique_ptr<float[]> clock(new float[ns]);

    for (unsigned i = 0; i < ns; ++i) {
        output[i] = input[i];
        clock[i] = fbbd / fs;
    }

    line.process(ns, output.get(), clock.get());

    //
    if (!output_path) {
        fprintf(stderr, "! Output file not specified, not saving result.\n");
    }
    else {
        SF_INFO info_out = {};
        info_out.samplerate = fs;
        info_out.frames = ns;
        info_out.channels = 1;
        info_out.format = SF_FORMAT_AIFF|SF_FORMAT_FLOAT;

        SNDFILE *sf_out = sf_open(output_path, SFM_WRITE, &info_out);

        if (!sf_out) {
            fprintf(stderr, "Cannot open the sound file for writing.\n");
            return 1;
        }

        fprintf(stderr, "* Writing result to: %s\n", output_path);

        sf_writef_float(sf_out, output.get(), ns);
        sf_close(sf_out);
    }

    return 0;
}

static void print_usage()
{
    fprintf(stderr,
            "bbd [OPTION]... <SOUND-FILE>\n"
            "    -o <SOUND-FILE>:  output result to AIFF file\n"
            "    -b <BUCKETS>:     set the number of buckets (default: %u)\n"
            "    -c <CLOCK>:       set the clock (default: %f)\n"
            "    -h:               display help\n",
            default_bucket_count,
            default_clock);
}

static std::vector<float> read_mono_samples(SNDFILE *sf, const SF_INFO &info)
{
    std::vector<float> samples;
    samples.reserve(info.frames);

    unsigned channels = info.channels;
    float mixdown_gain = M_SQRT1_2 / std::sqrt((M_SQRT1_2 * M_SQRT1_2) * channels);

    unsigned bufsize = 256;
    std::unique_ptr<float[]> buffer(new float[bufsize * channels]);

    for (sf_count_t count; (count = sf_readf_float(sf, buffer.get(), bufsize)) > 0;) {
        for (unsigned i = 0; i < (unsigned)count; ++i) {
            float sample = 0;
            for (unsigned j = 0; j < channels; ++j)
                sample += buffer[i * channels + j];
            samples.push_back(sample * mixdown_gain);
        }
    }

    return samples;
}
