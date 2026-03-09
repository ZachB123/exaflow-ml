#include <fstream>
#include <iomanip>

#include "metadata_writer.h"


void MetadataWriter::write(const std::filesystem::path& path) const {
    nlohmann::json metadata;
    for (const std::reference_wrapper<const IMetadataProvider>& provider : providers) {
        provider.get().appendMetadata(metadata);
    }
    std::ofstream out(path);
    out << std::setw(4) << metadata << "\n";
}