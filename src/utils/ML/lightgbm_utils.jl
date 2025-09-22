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
    LightGBMModel

LightGBM booster wrapper storing the trained booster and feature metadata.
"""
struct LightGBMModel
    booster::Any
    features::Vector{Symbol}
end

"""
Vector alias for collections of `LightGBMModel` instances.
"""
const LightGBMModelVector = Vector{LightGBMModel}

"""
    lightgbm_feature_matrix(psms, features) -> Matrix{Float32}

Convert the requested feature columns from a `DataFrame` into a dense
`Float32` matrix. Missing values are imputed with sensible defaults based
on the underlying column type.
"""
function lightgbm_feature_matrix(psms::AbstractDataFrame, features::Vector{Symbol})
    n = nrow(psms)
    m = length(features)
    matrix = Matrix{Float32}(undef, n, m)

    for (j, feat) in enumerate(features)
        column = psms[!, feat]
        T = nonmissingtype(eltype(column))

        if T <: AbstractFloat
            if eltype(column) <: Union{Missing, T}
                matrix[:, j] = Float32.(coalesce.(column, zero(T)))
            else
                matrix[:, j] = Float32.(column)
            end
        elseif T <: Integer
            if eltype(column) <: Union{Missing, T}
                matrix[:, j] = Float32.(coalesce.(column, zero(T)))
            else
                matrix[:, j] = Float32.(column)
            end
        elseif T <: Bool
            if eltype(column) <: Union{Missing, Bool}
                matrix[:, j] = Float32.(coalesce.(column, false))
            else
                matrix[:, j] = Float32.(column)
            end
        else
            throw(ArgumentError("Unsupported feature type $(eltype(column)) for LightGBM"))
        end
    end

    return matrix
end

function _prepare_labels(labels)
    result = Vector{Float64}(undef, length(labels))
    @inbounds for (idx, value) in pairs(labels)
        if value isa Missing
            result[idx] = 0.0
        elseif value isa Bool
            result[idx] = value ? 1.0 : 0.0
        else
            result[idx] = Float64(value)
        end
    end
    return result
end

"""
    train_lightgbm_booster(psms, features, num_round; kwargs...) -> LightGBMModel

Train a LightGBM booster on the provided PSM table using the selected
features. Hyperparameters mirror the prior EvoTrees interface for
backwards compatibility.
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
    available_features = Symbol[]
    for feat in features
        if hasproperty(psms, feat)
            push!(available_features, feat)
        end
    end

    if isempty(available_features)
        throw(ArgumentError("No available features for LightGBM training"))
    end

    matrix = lightgbm_feature_matrix(psms, available_features)
    labels = _prepare_labels(psms[!, target_name])

    dataset = _create_lightgbm_dataset(matrix, labels, available_features)

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

    if isdefined(LightGBM, :free_dataset!)
        try
            LightGBM.free_dataset!(dataset)
        catch
        end
    end

    return LightGBMModel(booster, available_features)
end

function _create_lightgbm_dataset(
    matrix::AbstractMatrix{Float32},
    labels::AbstractVector{Float64},
    features::Vector{Symbol},
)
    return LightGBM.Dataset(matrix; label = labels, feature_name = String.(features))
end

"""
    predict(model::LightGBMModel, psms) -> Vector{Float32}

Generate probability predictions for the provided PSMs.
"""
function predict(model::LightGBMModel, psms::AbstractDataFrame)
    matrix = lightgbm_feature_matrix(psms, model.features)
    return Float32.(LightGBM.predict(model.booster, matrix))
end

"""
    feature_importances(model; importance_type = :gain)

Return feature importances paired with their feature names.
"""
function feature_importances(
    model::LightGBMModel;
    importance_type::Symbol = :gain,
)
    scores = LightGBM.feature_importance(model.booster; importance_type = importance_type)
    return collect(zip(model.features, scores))
end

"""
    importance(model)

Alias for `feature_importances` to match the historical wrapper API.
"""
importance(model::LightGBMModel) = feature_importances(model)
