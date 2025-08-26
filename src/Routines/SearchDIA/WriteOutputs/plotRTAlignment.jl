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

function plotRTAlign(RT::Vector{T},
                    iRT::Vector{T},
                    rt_map::Any;
                    out_fdir::String = "./",
                    out_fname::String = "rt_align_plot") where {T<:AbstractFloat}
    n = length(RT)

    plot_title = ""
    n = 0
    for i in range(1, length(out_fname))
        n += 1
        if n > 24
            n = 1
            plot_title *= "\n"
        end
        plot_title *= out_fname[i]
    end

    # Original library iRTs
    p_orig = Plots.plot(RT, iRT, seriestype = :scatter,
                        title = plot_title*"\n n = $n",
                        xlabel = "Retention Time RT (min)",
                        ylabel = "Indexed Retention Time iRT (min)",
                        label = nothing,
                        size = 100*[13.3, 7.5],
                        fontsize = 24,
                        titlefontsize = 24,
                        legendfontsize = 24,
                        tickfontsize = 24,
                        guidefontsize = 24,
                        margin = 10Plots.mm,
                        alpha = 0.1,
                        dpi = 300)
    savefig(p_orig, joinpath(out_fdir, out_fname)*"_original.pdf")

    # Run-specific predicted iRTs using the alignment model
    aligned_iRT = rt_map.(RT)
    p_aligned = Plots.plot(RT, aligned_iRT, seriestype = :scatter,
                        title = plot_title*"\n n = $n",
                        xlabel = "Retention Time RT (min)",
                        ylabel = "Indexed Retention Time iRT (min)",
                        label = nothing,
                        size = 100*[13.3, 7.5],
                        fontsize = 24,
                        titlefontsize = 24,
                        legendfontsize = 24,
                        tickfontsize = 24,
                        guidefontsize = 24,
                        margin = 10Plots.mm,
                        alpha = 0.1,
                        dpi = 300)
    savefig(p_aligned, joinpath(out_fdir, out_fname)*"_aligned.pdf")

    return p_orig, p_aligned
end
