# Ksat data derived by
# Gupta, S., Lehmann, P., Bonetti, S., Papritz, A., and Or, D., (2020)
# Global prediction of soil saturated hydraulic conductivity using random
# forest in a Covariate-based Geo Transfer Functions (CoGTF) framework.
# Journal of Advances in Modeling Earth Systems, 13(4), e2020MS002242.
# https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020MS002242

# Water retention curve parameters derived by
# Gupta, S., Papritz, A., Lehmann, P., Hengl, T., Bonetti, S., & Or, D. (2022).
# Global Mapping of Soil Water Characteristics Parameters—Fusing Curated Data with
# Machine Learning and Environmental Covariates. Remote Sensing, 14(8), 1947.

# First, download the data in geotiff form using the script "download_soilgrids.sh"

# Then, convert each tif file to nc using the script "transform_geotiff_to_netcdf.sh" found here.
# Note that you must supply the SRC_DIR and DEST_DIR to match your local paths.

# You will also need the executable "gdal_translate".

# The following `filedir` is missing intentionally - it must be replaced with paths to the nc files
# on your local machine.

using NCDatasets
using Statistics
using ClimaArtifactsHelper

filedir = "soilgrids_nc/"

outputdir = "soilgrids_artifacts"
if isdir(outputdir)
    @warn "$outputdir already exists. Content will end up in the artifact and may be overwritten."
    @warn "Abort this calculation, unless you know what you are doing."
else
    mkdir(outputdir)
end

include("utils.jl")
# Parameters specific to this data
z = [-1.5, -0.8, -0.45, -0.225, -0.1, -0.025] # depth of soil layer
level_names = ["100-200cm","60-100cm","30-60cm","15-30cm","5-15cm","0-5cm"]
vars = ["bdod","silt","sand","clay","cfvo","soc"] # varnames
nvars = length(vars)
attrib_bdod = (;
    vartitle = "Dry bulk density of fine earth fraction",
    varunits = "kg/m^3",
    varname = "bdod",
)
transform_bdod(x) = x*1e-5*1e6 # how to convert to from cg/cm^3 to kg/m^3

attrib_silt = (;
    vartitle = "Mass fraction of silt in the mineral part of the fine earth fraction of soil",
    varunits = "kg/kg",
    varname = "Q_silt")
attrib_sand = (;
    vartitle = "Mass fraction of sand in the mineral part of the fine earth fraction of soil",
    varunits = "kg/kg",
    varname = "Q_sand")
attrib_clay = (;
    vartitle = "Mass fraction of clay in the mineral part of the fine earth fraction of soil",
    varunits = "kg/kg",
    varname = "Q_clay")
transform_comp(x) = x*1e-3 # how to convert to from g/kg to kg/kg

attrib_cfvo = (;
    vartitle = "Volumetric fraction of coarse fragments",
    varunits = "m^3/m^3",
    varname = "cfvo")
transform_cfvo(x) = x*1e-3 # how to convert to from (cm/dm)^3 to (m/m)^3

attrib_soc = (;
    vartitle = "Mass fraction of soil organic carbon in the fine earth fraction of soil",
    varunits = "kg/kg",
    varname = "q_soc")
transform_soc(x) = x * 1e-4 # how to convert to from dg/kg to kg/kg
attribs = [attrib_bdod, attrib_silt, attrib_sand, attrib_clay, attrib_cfvo, attrib_soc]
transforms = [transform_bdod, transform_comp, transform_comp, transform_comp, transform_cfvo, transform_soc]
nlayers = length(z)

# just pick one of the files to get lat and lon values
file= joinpath(filedir, "bdod_0-5cm_mean_5000.nc")
nc_data = NCDatasets.NCDataset(file)
lat = nc_data["lat"][:];
lon = nc_data["lon"][:];
lat_ct = length(lat)
lon_ct = length(lon)
data = Array{Union{Missing, Float32}}(missing, lon_ct, lat_ct, nlayers);

# Function which reads in the data by layer and writes the file with all layers to the correct output location.
function create_combined_data(data, files, attrib, transform, outfilepath)
    # get parameter values at each layer
    read_nc_data!(data, files, filedir)
    # Replace missing with NaN
    data[typeof.(data) .== Missing] .= NaN
    data .= transform(data)
    write_nc_out(data, lat, lon, z, attrib, outfilepath)
    Base.mv(outfilepath, joinpath(outputdir, outfilepath))
end

for i in 1:nvars
    var = vars[i]
    attrib = attribs[i]
    transform = transforms[i]
    files = ["$(var)_$(ln)_mean_5000.nc" for ln in level_names]
    outfilepath = "$(var)_soilgrids_combined.nc"
    create_combined_data(data, files, attrib, transform, outfilepath)
end

# Particle density
ρp_min = 2.65*1e3 # convert from g/cm^3 to kg/m^3
ρp_om = 1.3*1e3 # convert from g/cm^3 to kg/m^3
ρp_cf = ρp_min # Guess
# Check - lower case q_i = mass of i/mass of fine earth
#       - upper case Q_i = mass of i/mass of minerals in fine earth
#        - fine earth - particles < 2mm, this includes SOC
#        - coarse fragments - particles > 2mm
# q_soc + q_min = 1 (fine earth mass fractions)
# Q_i * (1-q_soc) = q_i (fine earth mass fraction of mineral i)
# ∑Q = Q_silt .+ Q_sand .+ Q_clay; # sums to 1 
soc = NCDataset("soilgrids_artifacts/soc_soilgrids_combined.nc");
q_soc = soc["q_soc"][:,:,:]
silt = NCDataset("soilgrids_artifacts/silt_soilgrids_combined.nc");
q_silt = silt["Q_silt"][:,:,:] .* (1 .- q_soc)
clay = NCDataset("soilgrids_artifacts/clay_soilgrids_combined.nc");
q_clay = clay["Q_clay"][:,:,:].* (1 .- q_soc)
sand = NCDataset("soilgrids_artifacts/sand_soilgrids_combined.nc");
q_sand = sand["Q_sand"][:,:,:].* (1 .- q_soc)
# Fine earth bulk density
density = NCDataset("soilgrids_artifacts/bdod_soilgrids_combined.nc");
ρ_bulk_fine_earth = density["bdod"][:,:,:] # Mass of fine earth/ Volume of whole soil

# Volumetric fractions relative to whole soil
cf = NCDataset("soilgrids_artifacts/cfvo_soilgrids_combined.nc");
θ_cf = cf["cfvo"][:,:,:] # volume of gravel/volume of whole soil
θ_silt = @. q_silt / ρp_min * ρ_bulk_fine_earth;
θ_clay = @. q_clay / ρp_min * ρ_bulk_fine_earth;
θ_sand = @. q_sand / ρp_min * ρ_bulk_fine_earth;
θ_om = @. q_soc / ρp_om * ρ_bulk_fine_earth;

# Sanity checks
check_extrema(x) = extrema(x[ .! isnan.(x)])
check_extrema(q_sand .+ q_soc .+ q_silt .+ q_clay)# [1,1]

θ_min = @. θ_silt + θ_clay + θ_sand
ν_est = 1 .- (θ_cf .+ θ_min .+ θ_om)
check_extrema(ν_est)

ρ_bulk_min = θ_min * ρp_min
ρ_bulk_om = θ_om * ρp_om
ρ_bulk_cf = θ_cf * ρp_cf
ρ_bulk_soil = (1-θ_cf)*ρ_bulk_fine_earth + θ_cf * ρ_bulk_cf
ν_formula = 1 .- ρ_bulk_soil * ()



create_artifact_guided(outputdir; artifact_name = basename(@__DIR__))
