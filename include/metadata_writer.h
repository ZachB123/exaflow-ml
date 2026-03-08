#pragma once

#include <filesystem>
#include <functional>
#include <vector>

#include <nlohmann/json.hpp>


class IMetadataProvider {
public:
    virtual void appendMetadata(nlohmann::json& metadata) const = 0;
    virtual ~IMetadataProvider() = default;
};


class MetadataWriter {
public:
    template<typename... Ps>
    MetadataWriter(const Ps&... ps)
        : providers{std::cref(static_cast<const IMetadataProvider&>(ps))...}
    {}

    void write(const std::filesystem::path& path) const;

private:
    std::vector<std::reference_wrapper<const IMetadataProvider>> providers;
};
