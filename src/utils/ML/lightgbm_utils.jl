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

    dataset = LightGBM.Dataset(
        feature_matrix;
        label = labels,
        feature_name = String.(features),
    )

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
