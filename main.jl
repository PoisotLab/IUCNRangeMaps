# We will use the polygon reading functions in SDT
using SpeciesDistributionToolkit
const SDT = SpeciesDistributionToolkit

# Needed for some sources
using STAC

# Save the data
using DataFrames
using Statistics
import CSV

# Function to get the rangemaps
function load_shapefile(shapefile_dir)
    Pol = SpeciesDistributionToolkit.SimpleSDMPolygons
    sf_table, _ = Pol._read(Pol._Shapefile, shapefile_dir)

    mps = Pol.GI.coordinates.(sf_table.geometry)
    taxa_labels = getproperty(sf_table, :sci_name)

    # We ONLY get the bats here - this needs to be abstracted away eventually
    chiros = findall(isequal("CHIROPTERA"), sf_table.order_)
    taxa_labels = taxa_labels[chiros]

    feats = []
    for sp in unique(taxa_labels)
        idx = findall(isequal(sp), taxa_labels)

        # We need to set the CRS because reasons?
        _crs = Pol.AG.importEPSG(4326)
        mp = MultiPolygon(Pol._add_crs(Pol.AG.createmultipolygon(mps[idx][1]), _crs))

        push!(feats, Feature(mp, Dict("Name" => sp)))
    end

    return FeatureCollection(feats)
end

# This loads the shapefile
iucn = load_shapefile("data")

# Get STAC data - these are examples
biab = STAC.Catalog("https://stac.geobon.org/")

# This is the source of layers we want -- they are in the STAC catalogue as collection
# / item / asset, so we just drop the names in a named tuple, which will make it much easier
# later on
layers_spec = [
    (collection="accessibility_to_cities", item="accessibility", asset="data"),
    (collection="bii_nhm", item="bii_nhm_10km_2020", asset="bii_nhm_10km_2020"),
]

# This will have the actual objects that we can pass through SDMLayer to get the data
biab_sources = [biab[s.collection].items[s.item].assets[s.data] for s in layers_spec]

# To create the dataframe, we need to have the type of the values as well as the column
# name. To avoid making a super WIDE dataframe, we will simply use a key-value system to
# have a LONG dataframe that we can then reshape
df = DataFrame(sp=AbstractString[], collection=AbstractString[], item=AbstractString[], asset=AbstractString[], statistic=AbstractString[], value=Real[])

stats = [mean, median, maximum, minimum, std]

# And now we can get started with the retrieval
for sp in uniqueproperties(iucn)["Name"]
    @info "🦇  $(sp)"
    _range = iucn[sp]
    _bbox = SDT.boundingbox(_range)
    for layer in layers_spec
        @info "\t 🗺️  $(layer.collection)"
        stac_asset = biab[layer.collection].items[layer.item].assets[layer.asset]
        L = SDMLayer(stac_asset; _bbox...)
        mask!(L, _range)
        if !iszero(count(L))
            for s in stats
                push!(df, (sp, layer.collection, layer.item, layer.asset, "$(s)", s(L)))
            end
        end
    end
end

CSV.write("new_variables.csv", df)
