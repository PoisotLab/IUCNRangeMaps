# We will use the polygon reading functions in SDT
using SpeciesDistributionToolkit
const SDT = SpeciesDistributionToolkit

# Needed for some sources
using STAC

# Save the data
using DataFrames
using Statistics
import Random
import CSV

# Make the data folder
_out_dir = "outputs"
if ~isdir(_out_dir)
    mkdir(_out_dir)
end

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
    # (collection="accessibility_to_cities", item="accessibility", asset="data"),
    (collection="bii_nhm", item="bii_nhm_10km_2020", asset="bii_nhm_10km_2020"),
    # (collection="ghmts", item="GHMTS", asset="GHMTS"),
]

stats = [mean, median, maximum, minimum, std]

function prepare_summary(ranges, species, catalogue, layers, measures)
    # This is the name of the expected output file
    fname = joinpath(_out_dir, replace(species, " " => "_") * ".csv")

    # If the output file is already here, we skip this species -- this is useful if we want
    # to re-start the process later without re-doing species for which we have the results
    if isfile(fname)
        return nothing
    end

    # We extract the range from the object containing the spatial features
    _range = ranges[species]

    # This is the bounding box of the range, so we store it once for all sources of
    # information on this species
    _bbox = SDT.boundingbox(_range)

    # We start an empty collection in which we will put named tuples that will become the
    # rows of a data frame
    outputs = []

    # Now for each layer -- actually defined as a named tuple for access to the STAC
    # catalogue -- we start extracting information
    for layer in layers
        # The first step is to get the data itself
        stac_asset = catalogue[layer.collection].items[layer.item].assets[layer.asset]
        L = SDMLayer(stac_asset; _bbox...)

        # And then we clip it with the bounding box of the species
        mask!(L, _range)

        # If there are no cells in the layer, which can happen with species with super
        # small ranges and layers with poor resolution, we don't do anything
        if !iszero(count(L))

            # For every measure we want to apply to the layer...
            for m in measures

                # We add the output of applying this measure, as well as information about
                # the species and data, to the outputs collection
                push!(outputs, (
                    species=species,
                    collection=layer.collection,
                    item=layer.item,
                    asset=layer.asset,
                    measure=string(m),
                    value=m(L)
                ))
            end
        end
    end

    # If we have at least one entry in outputs...
    if ~isempty(outputs)
        # We change the outputs into a data frame
        df = DataFrame(outputs)
        # We save it
        CSV.write(fname, df)
        # And we return the filename
        return fname
    end
    return nothing
end

# Demo
#prepare_summary(iucn, "Artibeus jamaicensis", biab, layers_spec, stats)

# We do the species in random order
sp_names = Random.shuffle(uniqueproperties(iucn)["Name"])

# Now we run these chunks in parallel
tasks = fetch.([Threads.@spawn prepare_summary(iucn, sp, biab, layers_spec, stats) for sp in sp_names])
