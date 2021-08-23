
#include <savvy/reader.hpp>
#include <savvy/writer.hpp>

#include <cstdint>
#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <algorithm>
#include <getopt.h>

std::vector<std::string> split_string_to_vector(const char* in, char delim)
{
  std::vector<std::string> ret;
  const char* d = nullptr;
  std::string token;
  const char* s = in;
  const char*const e = in + strlen(in);
  while ((d = std::find(s, e,  delim)) != e)
  {
    ret.emplace_back(std::string(s, d));
    s = d ? d + 1 : d;
  }
  ret.emplace_back(std::string(s,d));
  return ret;
}

savvy::genomic_region string_to_region(const std::string& s)
{
  const std::size_t colon_pos = s.find(':');
  if (colon_pos == std::string::npos)
  {
    return savvy::genomic_region(s);
  }
  else
  {
    std::string chr = s.substr(0, colon_pos);
    const std::size_t hyphen_pos = s.find('-', colon_pos + 1);
    if (hyphen_pos == std::string::npos)
    {
      std::string slocus = s.substr(colon_pos + 1);
      std::uint64_t ilocus = std::uint64_t(std::atoll(slocus.c_str()));
      return savvy::genomic_region(chr, ilocus, ilocus);
    }
    else
    {
      std::string sbeg = s.substr(colon_pos + 1, hyphen_pos - chr.size() - 1);
      std::string send = s.substr(hyphen_pos + 1);
      if (send.empty())
      {
        return savvy::genomic_region(chr, std::uint64_t(std::atoll(sbeg.c_str())));
      }
      else
      {
        return savvy::genomic_region(chr, std::uint64_t(std::atoll(sbeg.c_str())), std::uint64_t(std::atoll(send.c_str())));
      }
    }
  }

}

class prog_args
{
private:
  static const int default_compression_level = 3;
  static const int default_block_size = 2048;

  std::vector<option> long_options_;
  std::string input_path_;
  std::string patch_input_path_;
  std::vector<savvy::region> patch_regions_;
  std::string output_path_ = "/dev/stdout";
  savvy::file::format output_mode_ = savvy::file::format::vcf;
  int compression_level_ = 0;
  bool help_ = false;

public:
  prog_args() :
    long_options_(
      {
        {"help", no_argument, 0, 'h'},
        {"output", required_argument, 0, 'o'},
        {"output-format", required_argument, 0, 'O'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::string& input_path() const { return input_path_; }
  const std::string& patch_input_path() const { return patch_input_path_; }
  const std::string& output_path() const { return output_path_; }
  savvy::file::format output_mode() const { return output_mode_; }
  const std::vector<savvy::region>& patch_regions() const { return patch_regions_; }
  int compression_level() const { return compression_level_; }
  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "Usage: vcf-patch [opts ...] <in.bcf> <patch.bcf> <regions_string>\n";
    os << "\n";
    os << " -h, --help           Print usage\n";
    os << " -o, --output         Output path (default: /dev/stdout)\n";
    os << " -O, --output-format  Output format (default: vcf)\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "ho:O:", long_options_.data(), &long_index )) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
      case 'h':
        help_ = true;
        return true;
      case 'o':
        output_path_ = optarg ? optarg : "";
        break;
      case 'O':
      {
        using fmt = savvy::file::format;
        std::string ot = optarg ? optarg : "";
        if (ot == "vcf")
        {
          output_mode_ = fmt::vcf;
        }
        else if (ot == "vcf.gz")
        {
          output_mode_ = fmt::vcf;
          compression_level_ = 6;
        }
        else if (ot == "bcf")
        {
          output_mode_ = fmt::bcf;
          compression_level_ = 6;
        }
        else if (ot == "ubcf")
        {
          output_mode_ = fmt::bcf;
        }
        else
        {
          std::cerr << "Invalid --output-format: " << ot << std::endl;
          return false;
        }
        break;
      }
      default:
        return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count == 3)
    {
      input_path_ = argv[optind];
      patch_input_path_ = argv[optind + 1];
      auto reg_vec = split_string_to_vector(argv[optind + 2], ',');
      patch_regions_.reserve(reg_vec.size());
      for (const auto& r : reg_vec)
        patch_regions_.emplace_back(string_to_region(r));
    }
    else if (remaining_arg_count < 3)
    {
      std::cerr << "Too few arguments\n";
      return false;
    }
    else
    {
      std::cerr << "Too many arguments\n";
      return false;
    }

    if (input_path_ == "/dev/stdin" || input_path_ == "/dev/fd/0")
    {
//      std::cerr << "Input file cannot be stdin\n";
//      return false;
    }

    return true;
  }
};

int main(int argc, char** argv)
{
  prog_args args;
  if (!args.parse(argc, argv))
  {
    args.print_usage(std::cerr);
    return EXIT_FAILURE;
  }

  if (args.help_is_set())
  {
    args.print_usage(std::cout);
    return EXIT_SUCCESS;
  }

  savvy::reader in(args.input_path().c_str());
  if (!in)
    return std::fprintf(stderr, "Failed to open \"%s\" : %s", args.input_path().c_str(), strerror(errno)), EXIT_FAILURE;

  auto hdrs = in.headers();
  hdrs.emplace_back("INFO","<ID=PATCH_VARIANT,Number=0,Type=Flag,Description=\"Patched variant\">");

  std::unordered_set<std::string> out_info_keys {"PATCH_VARIANT"};
  for (const auto& h : in.info_headers())
    out_info_keys.insert(h.id);

  savvy::writer out(args.output_path(), args.output_mode(), hdrs, in.samples(), args.compression_level());
  if (!out)
    return fprintf(stderr, "Couldn't open \"%s\" for writing: %s\n", args.output_path().c_str(), strerror(errno)), EXIT_FAILURE;

  savvy::variant rec;

  in >> rec;

  for (const auto& reg : args.patch_regions())
  {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    // Copy records before region and discard those in region
    while (in && rec.position() <= reg.to())
    {
      if (rec.position() < reg.from())
      {
        if (!(out << rec))
          return std::fprintf(stderr, "Failed to write to %s\n", args.output_path().c_str()), EXIT_FAILURE;
      }
      in >> rec;
    }
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    // Copy region from patch file
    savvy::reader in_patch(args.patch_input_path());
    in_patch.reset_bounds(reg);
    if (!in_patch)
      return std::fprintf(stderr, "Failed to query %s from %s : %s", (reg.chromosome() + ":" + std::to_string(reg.from()) + "-" + std::to_string(reg.to())).c_str(), args.input_path().c_str(), strerror(errno)), EXIT_FAILURE;


    savvy::variant rec_patch;
    while (in_patch >> rec_patch && rec_patch.position() >= reg.from() && rec_patch.position() < reg.to())
    {
      if (in_patch.bad())
        return std::fprintf(stderr, "Failed to read from %s\n", args.patch_input_path().c_str()), EXIT_FAILURE;

      std::vector<std::string> info_fields_to_remove;
      for (const auto& f : rec_patch.info_fields())
      {
        if (out_info_keys.find(f.first) == out_info_keys.end())
          info_fields_to_remove.emplace_back(f.first);
      }

      for (const auto& k : info_fields_to_remove)
        rec_patch.remove_info(k);

      rec_patch.set_info("PATCH_VARIANT", std::vector<std::int8_t>());

      if (!(out << rec_patch))
        return std::fprintf(stderr, "Failed to write to %s\n", args.output_path().c_str()), EXIT_FAILURE;
    }
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  // Copy records after last region
  while (in)
  {
    if (!(out << rec))
      return std::fprintf(stderr, "Failed to write to %s\n", args.output_path().c_str()), EXIT_FAILURE;
    in >> rec;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

  return out.good() && !in.bad() ? EXIT_SUCCESS : EXIT_FAILURE;
}
