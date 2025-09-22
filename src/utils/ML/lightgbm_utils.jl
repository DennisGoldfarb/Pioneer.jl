# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""
    PioneerLightGBMModel

LightGBM booster wrapper storing feature metadata and missing value fill values.
"""
struct PioneerLightGBMModel
    booster::LightGBM.Booster
    feature_names::Vector{Symbol}
    fill_values::Vector{Float32}
end

"""
    _convert_label(value) -> Float32

Convert labels to `Float32`, treating booleans as 1/0 and missing values as 0.
"""
@inline function _convert_label(value)
    if value isa Missing
        return 0f0
    elseif value isa Bool
        return value ? 1f0 : 0f0
    else
        return Float32(value)
    end
end

function _fill_value(column)
    collected = Float32[]
    for value in skipmissing(column)
        push!(collected, Float32(value))
    end
    return isempty(collected) ? 0f0 : Float32(median(collected))
end

"""
    prepare_feature_matrix(data, features) -> Matrix{Float32}, Vector{Float32}

Convert feature columns to a dense matrix while imputing missing values with
their column medians. Returns the matrix and the fill values used.
"""
function prepare_feature_matrix(data::AbstractDataFrame, features::Vector{Symbol})
    n_rows = nrow(data)
    n_features = length(features)
    matrix = Matrix{Float32}(undef, n_rows, n_features)
    fill_values = Vector{Float32}(undef, n_features)

    for (j, feature) in enumerate(features)
        column = data[!, feature]
        fill = _fill_value(column)
        fill_values[j] = fill

        @inbounds for i in 1:n_rows
            value = column[i]
            matrix[i, j] = ismissing(value) ? fill : Float32(value)
        end
    end

    return matrix, fill_values
end

"""
    prepare_feature_matrix(data, features, fill_values) -> Matrix{Float32}

Convert feature columns to a matrix using previously computed fill values.
"""
function prepare_feature_matrix(
    data::AbstractDataFrame,
    features::Vector{Symbol},
    fill_values::Vector{Float32},
)
    n_rows = nrow(data)
    n_features = length(features)
    matrix = Matrix{Float32}(undef, n_rows, n_features)

    @inbounds for (j, feature) in enumerate(features)
        column = data[!, feature]
        fill = fill_values[j]

        for i in 1:n_rows
            value = column[i]
            matrix[i, j] = ismissing(value) ? fill : Float32(value)
        end
    end

    return matrix
end

"""
    train_lightgbm_booster(psms, features, num_round; kwargs...) -> PioneerLightGBMModel

Train a LightGBM booster on the provided PSM table using the selected features.
Hyperparameters mirror the previous EvoTrees interface for compatibility.
"""
function train_lightgbm_booster(
    psms::AbstractDataFrame,
    features::Vector{Symbol},
    num_round::Integer;
    target_name::Symbol = :target,
    colsample::Float64 = 1.0,
    eta::Float64 = 0.1,
    min_child_weight::Integer = 1,
    subsample::Float64 = 1.0,
    gamma::Float64 = 0.0,
    max_depth::Integer = -1,
    verbosity::Integer = 0,
)
    feature_matrix, fill_values = prepare_feature_matrix(psms, features)

    labels_raw = psms[!, target_name]
    labels = Vector{Float32}(undef, length(labels_raw))
    @inbounds for (i, value) in pairs(labels_raw)
        labels[i] = _convert_label(value)
    end

    dataset = _create_lightgbm_dataset(feature_matrix, labels, features)

    params = Dict{String, Any}(
        "objective" => "binary",
        "metric" => "binary_logloss",
        "learning_rate" => eta,
        "feature_fraction" => colsample,
        "bagging_fraction" => subsample,
        "bagging_freq" => 1,
        "max_depth" => max_depth,
        "min_child_weight" => min_child_weight,
        "min_split_gain" => gamma,
        "verbosity" => verbosity > 0 ? verbosity : -1,
    )

    booster = LightGBM.train(params, dataset; num_boost_round = Int(num_round))
    return PioneerLightGBMModel(booster, Symbol.(features), fill_values)
end

"""
    _create_lightgbm_dataset(feature_matrix, labels, features)

Construct a LightGBM dataset using the keyword-based API introduced in
LightGBM.jl 2.0.0, falling back to the positional constructor used by
earlier releases. This keeps Pioneer compatible across LightGBM versions
without requiring callers to care about the specific signature.
"""
function _create_lightgbm_dataset(
    feature_matrix::AbstractMatrix{Float32},
    labels::AbstractVector{Float32},
    features::Vector{Symbol},
)
    feature_names = String.(features)
    attempts = (
        ("keyword-only constructor", () -> LightGBM.Dataset(; data = feature_matrix, label = labels, feature_name = feature_names)),
        ("positional data with feature_name keyword", () -> LightGBM.Dataset(feature_matrix; label = labels, feature_name = feature_names)),
        ("positional data with feature_names keyword", () -> LightGBM.Dataset(feature_matrix; label = labels, feature_names = feature_names)),
        ("positional data, label keyword only", () -> LightGBM.Dataset(feature_matrix; label = labels)),
        ("positional data and label", () -> LightGBM.Dataset(feature_matrix, labels)),
        ("positional data, label, and feature_name", () -> LightGBM.Dataset(feature_matrix, labels, feature_names)),
    )

    collected_errors = Vector{Pair{String, Any}}()

    for (label, attempt) in attempts
        try
            return attempt()
        catch err
            if err isa MethodError || err isa ArgumentError
                push!(collected_errors, label => err)
            else
                rethrow(err)
            end
        end
    end

    if _has_lightgbm_low_level_api()
        try
            return _create_dataset_via_low_level(feature_matrix, labels, feature_names)
        catch err
            push!(collected_errors, "low-level LightGBM C API" => err)
        end
    end

    if !isempty(collected_errors)
        messages = map(collected_errors) do (label, err)
            sprint() do io
                print(io, label, " failed with ")
                showerror(io, err)
            end
        end
        error(
            "Failed to construct LightGBM dataset using available constructors.\n" *
            join(messages, "\n"),
        )
    else
        error("Failed to construct LightGBM dataset from the provided feature matrix.")
    end
end

function _has_lightgbm_low_level_api()
    return all(
        name -> isdefined(LightGBM, name),
        (:LGBM_DatasetCreateFromMat, :LGBM_DatasetSetField, :LGBM_DatasetSetFeatureNames),
    )
end

const _LIGHTGBM_FLOAT32 = isdefined(LightGBM, :C_API_DTYPE_FLOAT32) ? LightGBM.C_API_DTYPE_FLOAT32 : Cint(0)

function _create_dataset_via_low_level(
    feature_matrix::AbstractMatrix{Float32},
    labels::AbstractVector{Float32},
    feature_names::Vector{String},
)
    handle_ref = Ref{Ptr{Nothing}}()
    status = _call_lgbm_dataset_create_from_mat(feature_matrix, handle_ref)
    _check_lightgbm_status(status, "creating dataset from matrix")

    dataset = LightGBM.Dataset(handle_ref[])

    status = _call_lgbm_dataset_set_field(dataset, "label", labels, _LIGHTGBM_FLOAT32)
    _check_lightgbm_status(status, "attaching labels")

    if !isempty(feature_names)
        status = _call_lgbm_dataset_set_feature_names(dataset, feature_names)
        _check_lightgbm_status(status, "setting feature names")
    end

    return dataset
end

function _call_lgbm_dataset_create_from_mat(
    feature_matrix::AbstractMatrix{Float32},
    handle_ref::Ref{Ptr{Nothing}},
)
    try
        return LightGBM.LGBM_DatasetCreateFromMat(
            feature_matrix,
            _LIGHTGBM_FLOAT32,
            size(feature_matrix, 1),
            size(feature_matrix, 2),
            0,
            "",
            nothing,
            handle_ref,
        )
    catch err
        if err isa MethodError
            return LightGBM.LGBM_DatasetCreateFromMat(
                pointer(feature_matrix),
                _LIGHTGBM_FLOAT32,
                size(feature_matrix, 1),
                size(feature_matrix, 2),
                0,
                "",
                Ptr{Nothing}(C_NULL),
                handle_ref,
            )
        else
            rethrow(err)
        end
    end
end

function _call_lgbm_dataset_set_field(
    dataset,
    field_name::AbstractString,
    values::AbstractVector{Float32},
    dtype,
)
    try
        return LightGBM.LGBM_DatasetSetField(
            dataset,
            field_name,
            values,
            length(values),
            dtype,
        )
    catch err
        if err isa MethodError
            handle = getfield(dataset, :handle)
            return LightGBM.LGBM_DatasetSetField(
                handle,
                field_name,
                pointer(values),
                length(values),
                dtype,
            )
        else
            rethrow(err)
        end
    end
end

function _call_lgbm_dataset_set_feature_names(dataset, feature_names::Vector{String})
    try
        return LightGBM.LGBM_DatasetSetFeatureNames(dataset, feature_names)
    catch err
        if err isa MethodError
            handle = getfield(dataset, :handle)
            c_feature_names = Base.cconvert(Vector{Cstring}, feature_names)
            return LightGBM.LGBM_DatasetSetFeatureNames(
                handle,
                c_feature_names,
                length(c_feature_names),
            )
        else
            rethrow(err)
        end
    end
end

function _check_lightgbm_status(status::Integer, context::AbstractString)
    status == 0 && return
    last_error = _lightgbm_last_error()
    if isempty(last_error)
        error("LightGBM $context failed with status $status")
    else
        error("LightGBM $context failed with status $status: $last_error")
    end
end

function _lightgbm_last_error()
    if isdefined(LightGBM, :LGBM_GetLastError)
        ptr = LightGBM.LGBM_GetLastError()
        if ptr !== C_NULL
            return unsafe_string(ptr)
        end
    end
    return ""
end

"""
    predict(model::PioneerLightGBMModel, data) -> Vector{Float32}

Generate probabilities for the provided data using a trained LightGBM model.
"""
function predict(model::PioneerLightGBMModel, data::AbstractDataFrame)
    feature_matrix = prepare_feature_matrix(data, model.feature_names, model.fill_values)
    preds = LightGBM.predict(model.booster, feature_matrix)
    return Float32.(preds)
end

"""
    feature_importances(model; importance_type="gain")

Return feature importances paired with their feature names.
"""
function feature_importances(
    model::PioneerLightGBMModel;
    importance_type::AbstractString = "gain",
)
    importances = LightGBM.feature_importance(
        model.booster;
        importance_type = importance_type,
    )
    return collect(zip(model.feature_names, importances))
end
